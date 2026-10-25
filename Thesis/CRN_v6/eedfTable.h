// EEDF table shared by ChemPlasKin (masterCRN.cpp) and EEDFTable (eedfTableBuilder.cpp).
// One point = one converged CppBOLOS solution for (E/N, T_gas, mixture of the Boltzmann species)
// with the outputs the reactors need: Te, power/N and the rate constant of every Boltzmann process.
// The file is YAML: a header with the solver settings and a "points" list, one flow-style line per
// point, so that new points are simply appended at the end of the file.
#ifndef EEDF_TABLE_H
#define EEDF_TABLE_H

#include "cantera/base/AnyMap.h"
#include "cantera/base/ct_defs.h"
#include "cantera/ext/bolos/solver.h" // also includes parser.h (no include guard)

#include <algorithm>
#include <cctype>
#include <deque>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

namespace EEDF {

struct Point {
    double EN = 0.0;                   // reduced electric field [Td]
    double T = 0.0;                    // gas temperature [K]
    std::map<std::string, double> X;   // mole fractions of the Boltzmann species
    double Te = 0.0;                   // electron temperature from CppBOLOS [K] (not clipped to T)
    double power = 0.0;                // power/N absorbed by electrons [eV m^3/s]
    std::map<std::string, double> k;   // rate constant of each Boltzmann process [m^3/kmol/s]
};

// Settings that must be identical for all the points of a table
struct Settings {
    std::string csDataFile;
    std::string gridType;
    double gridStart = 0.0;
    double gridEnd = 0.0;
    int gridPoints = 0;
    int maxItr = 200;
    double rtol = 1e-5;
};

// Matching of stored points. A request is served by linear interpolation in T between two
// points with the same E/N and mixture at most maxGap apart, or, when T is not between two such
// points, by the closest point within T (TGAS_TOLERANCE, the error CRN_v4 accepts between two
// EEDF updates). T is also the distance under which the builder considers a point already present.
struct Tolerance {
    double EN = 0.5;    // [Td], EN_TOLERANCE
    double T = 0.5;     // [K], TGAS_TOLERANCE
    double maxGap = 5.0; // [K], largest T interval used for interpolation (0 = exact T only, < 0 = no reuse)
    double Xrel = 0.0;  // mixture: |x - x_point| <= max(Xabs, Xrel*x_point) for every species
    double Xabs = 0.0;  // (0 and 0 = identical mixture only)
};

// Grid and cross sections, as in the CRN_v4 set-up of BoltzmannRate::bsolver
inline void setupSolver(CppBOLOS::BoltzmannSolver& bs, const Settings& set)
{
    auto lines = CppBOLOS::clean_file(set.csDataFile);
    std::vector<CppBOLOS::Collision> collisions = CppBOLOS::parse(lines);
    bs.set_grid(set.gridType, set.gridStart, set.gridEnd, set.gridPoints);
    bs.load_collisions(collisions);
}

// Converge the EEDF for (EN, T, X) and evaluate the rate constants of the given processes.
// Same sequence as CRN_v4: F0 is the initial guess (2 eV Maxwellian if empty, Maxwellian at T_gas
// if ENguess < 1 Td) and returns the converged EEDF; delta0 is the first pseudo-time step of converge().
// X must contain all the Boltzmann species: the solver keeps the density of a species not set here.
inline Point solvePoint(CppBOLOS::BoltzmannSolver& bs, const Settings& set, double EN, double T,
                        const std::map<std::string, double>& X,
                        const std::vector<std::string>& processes, Eigen::VectorXd& F0,
                        double ENguess, double delta0)
{
    bs.set_EN(EN);
    bs.set_kT(T);
    bs.set_density(X);
    bs.init();

    if (F0.size() == 0) {
        F0 = bs.maxwell(2.0);
    } else if (ENguess < 1) {
        F0 = bs.maxwell(T * 8.6173E-5);
    }
    F0 = bs.converge(F0, set.maxItr, set.rtol, delta0);

    Point p;
    p.EN = bs.get_EN(); // values as seen by the solver, as CRN_v4
    p.T = bs.get_kT();
    p.X = X;
    p.Te = bs.get_Te();
    p.power = bs.elec_power();
    for (const auto& process : processes) {
        if (p.k.count(process) == 0) {
            p.k[process] = Cantera::Avogadro * bs.rate(F0, process);
        }
    }
    return p;
}

class Table
{
public:
    // path empty: in-memory table only. Otherwise the file is loaded if it exists (its settings
    // must match set) or created with the header, and every new point is appended to it.
    Table(const std::string& path, const Settings& set, const Tolerance& tolerance)
        : tol(tolerance), m_path(path), m_set(set)
    {
        if (m_path.empty()) return;
        if (std::filesystem::exists(m_path) && std::filesystem::file_size(m_path) > 0) {
            load();
        } else {
            const auto dir = std::filesystem::path(m_path).parent_path();
            if (!dir.empty()) std::filesystem::create_directories(dir);
            writeHeader();
        }
        m_nLoaded = m_points.size();
    }

    // Closest stored point within tolerance that contains all the processes, nullptr if none
    const Point* find(double EN, double T, const std::map<std::string, double>& X,
                      const std::vector<std::string>& processes) const
    {
        const Point* best = nullptr;
        double bestScore = std::numeric_limits<double>::max();
        for (const Group* gp : groupsNear(EN)) {
            const Group& g = *gp;
            const double dX = groupDistance(g, EN, X, processes);
            if (dX < 0) continue;
            for (const Point* p : neighbours(g, T)) {
                const double dT = std::abs(p->T - T);
                if (dT > tol.T || !hasProcesses(*p, processes)) continue;
                const double score = (tol.T > 0 ? dT / tol.T : 0.0) + (tol.EN > 0 ? std::abs(g.EN - EN) / tol.EN : 0.0) + dX;
                if (score < bestScore) {
                    bestScore = score;
                    best = p;
                }
            }
        }
        return best;
    }

    // Point for (EN, T, X), from the points of one E/N and mixture (within tolerance of the request):
    // a point at exactly T, else the linear interpolation in T between the two points around T if
    // at most maxGap apart, else the closest point within T. Empty if none.
    std::optional<Point> interpolate(double EN, double T, const std::map<std::string, double>& X,
                                     const std::vector<std::string>& processes) const
    {
        if (tol.maxGap < 0) return std::nullopt; // every EEDF computed, as CRN_v4
        const Point* lo = nullptr;
        const Point* hi = nullptr;
        const Point* nearest = nullptr;
        double bestBracket = std::numeric_limits<double>::max();
        double bestNearest = std::numeric_limits<double>::max();
        for (const Group* gp : groupsNear(EN)) {
            const Group& g = *gp;
            const double dX = groupDistance(g, EN, X, processes);
            if (dX < 0) continue;
            const double dEN = tol.EN > 0 ? std::abs(g.EN - EN) / tol.EN : 0.0;
            auto [a, b] = bracket(g, T);
            if (a && !hasProcesses(*a, processes)) a = nullptr;
            if (b && !hasProcesses(*b, processes)) b = nullptr;
            if (a && a->T == T) return *a;
            if (a && b && b->T - a->T <= tol.maxGap) {
                const double score = (b->T - a->T) / std::max(tol.maxGap, 1e-300) + dX + dEN;
                if (score < bestBracket) {
                    bestBracket = score;
                    lo = a;
                    hi = b;
                }
            }
            for (const Point* p : {a, b}) {
                if (!p || std::abs(p->T - T) > tol.T) continue;
                const double score = (tol.T > 0 ? std::abs(p->T - T) / tol.T : 0.0) + dX + dEN;
                if (score < bestNearest) {
                    bestNearest = score;
                    nearest = p;
                }
            }
        }
        if (!lo) {
            if (!nearest) return std::nullopt;
            return *nearest;
        }

        const double w = (T - lo->T) / (hi->T - lo->T);
        auto lin = [w](double a, double b) { return a + w * (b - a); };
        Point p;
        p.EN = lo->EN;
        p.T = T;
        p.X = lo->X;
        p.Te = lin(lo->Te, hi->Te);
        p.power = lin(lo->power, hi->power);
        for (const auto& process : processes) {
            p.k[process] = lin(lo->k.at(process), hi->k.at(process));
        }
        return p;
    }

    // Add a point to the table (and to the file). The reference is valid until the next add().
    const Point& add(const Point& p)
    {
        m_points.push_back(p);
        index(m_points.size() - 1);
        if (!m_path.empty()) {
            // single write() in append mode: lines of different processes sharing the file do not mix
            const std::string line = formatPoint(p);
            const int fd = ::open(m_path.c_str(), O_WRONLY | O_APPEND | O_CREAT, 0644);
            if (fd < 0 || ::write(fd, line.data(), line.size()) != static_cast<ssize_t>(line.size())) {
                if (fd >= 0) ::close(fd);
                throw std::runtime_error("EEDF table: cannot append to " + m_path);
            }
            ::close(fd);
        }
        return m_points.back();
    }

    const std::deque<Point>& points() const { return m_points; }
    size_t size() const { return m_points.size(); }
    size_t nLoaded() const { return m_nLoaded; }
    const std::string& path() const { return m_path; }

    Tolerance tol;

private:
    std::string m_path;
    Settings m_set;
    std::deque<Point> m_points; // deque: stable addresses for the index

    // Two-level index. Level 1: E/N -> groups with that E/N (ordered map, so the groups within
    // EN_TOLERANCE of a request are found by binary search). Level 2: a group holds the points of
    // one E/N and one mixture sorted by T (binary search of the two points around T).
    struct Group {
        double EN;
        std::map<std::string, double> X;
        std::vector<const Point*> byT;
    };
    std::deque<Group> m_groups; // deque: stable addresses for m_groupsByEN
    std::map<double, std::vector<Group*>> m_groupsByEN;
    size_t m_nLoaded = 0;

    void index(size_t i)
    {
        const Point* p = &m_points[i];
        auto& sameEN = m_groupsByEN[p->EN];
        auto it = std::find_if(sameEN.begin(), sameEN.end(), [&](const Group* gr) { return gr->X == p->X; });
        Group* g;
        if (it == sameEN.end()) {
            m_groups.push_back({p->EN, p->X, {}});
            g = &m_groups.back();
            sameEN.push_back(g);
        } else {
            g = *it;
        }
        auto pos = std::upper_bound(g->byT.begin(), g->byT.end(), p->T,
                                    [](double T, const Point* q) { return T < q->T; });
        g->byT.insert(pos, p);
    }

    // Groups whose E/N is within EN_TOLERANCE of EN
    std::vector<const Group*> groupsNear(double EN) const
    {
        std::vector<const Group*> v;
        for (auto it = m_groupsByEN.lower_bound(EN - tol.EN); it != m_groupsByEN.end() && it->first <= EN + tol.EN; ++it) {
            v.insert(v.end(), it->second.begin(), it->second.end());
        }
        return v;
    }

    // Mixture distance of a group for a request, -1 if the group cannot serve it
    double groupDistance(const Group& g, double EN, const std::map<std::string, double>& X,
                         const std::vector<std::string>& processes) const
    {
        if (std::abs(g.EN - EN) > tol.EN || g.byT.empty()) return -1.0;
        const double dX = mixtureDistance(X, g.X);
        if (dX > 1.0) return -1.0;
        return dX;
    }

    // Last point with T_point <= T and first point with T_point > T (nullptr if missing)
    static std::pair<const Point*, const Point*> bracket(const Group& g, double T)
    {
        auto it = std::upper_bound(g.byT.begin(), g.byT.end(), T,
                                   [](double t, const Point* q) { return t < q->T; });
        const Point* a = (it != g.byT.begin()) ? *(it - 1) : nullptr;
        const Point* b = (it != g.byT.end()) ? *it : nullptr;
        return {a, b};
    }

    static std::vector<const Point*> neighbours(const Group& g, double T)
    {
        auto [a, b] = bracket(g, T);
        std::vector<const Point*> v;
        if (a) v.push_back(a);
        if (b) v.push_back(b);
        return v;
    }

    static bool hasProcesses(const Point& p, const std::vector<std::string>& processes)
    {
        for (const auto& process : processes) {
            if (p.k.count(process) == 0) return false;
        }
        return true;
    }

    // Largest species difference relative to the tolerance: <= 1 means same mixture.
    // Zero tolerances = exact match (0 for identical mixtures, infinity otherwise).
    double mixtureDistance(const std::map<std::string, double>& X, const std::map<std::string, double>& Xp) const
    {
        double d = 0.0;
        auto add = [&](double x, double xp) {
            const double diff = std::abs(x - xp);
            const double allowed = std::max(tol.Xabs, tol.Xrel * std::abs(xp));
            if (diff == 0.0) return;
            d = std::max(d, allowed > 0 ? diff / allowed : std::numeric_limits<double>::infinity());
        };
        for (const auto& [sp, xp] : Xp) {
            auto it = X.find(sp);
            add(it != X.end() ? it->second : 0.0, xp);
        }
        for (const auto& [sp, x] : X) {
            if (Xp.count(sp) == 0) add(x, 0.0);
        }
        return d;
    }

    static std::string num(double v)
    {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.17g", v); // round-trip exact: mixtures are matched exactly
        return buf;
    }

    static std::string formatMap(const std::map<std::string, double>& m)
    {
        std::string s = "{";
        bool first = true;
        for (const auto& [key, v] : m) {
            s += (first ? "\"" : ", \"") + key + "\": " + num(v);
            first = false;
        }
        return s + "}";
    }

    static std::string formatPoint(const Point& p)
    {
        return "  - {EN: " + num(p.EN) + ", T: " + num(p.T) + ", Te: " + num(p.Te)
               + ", power: " + num(p.power) + ", X: " + formatMap(p.X) + ", k: " + formatMap(p.k) + "}\n";
    }

    static double toDouble(const Cantera::AnyValue& v)
    {
        return v.is<long int>() ? static_cast<double>(v.asInt()) : v.asDouble();
    }

    std::uintmax_t csFileSize() const
    {
        return std::filesystem::exists(m_set.csDataFile) ? std::filesystem::file_size(m_set.csDataFile) : 0;
    }

    void writeHeader() const
    {
        std::ofstream out(m_path);
        if (!out.is_open()) throw std::runtime_error("EEDF table: cannot create " + m_path);
        out << "# EEDF table (CppBOLOS) for ChemPlasKin CRN_v5\n"
            << "# Inputs: EN [Td], T [K], X = mole fractions of the Boltzmann species.\n"
            << "# Outputs: Te [K], power = power/N [eV m^3/s], k = rate constants [m^3/kmol/s].\n"
            << "# New points are appended at the end of the file. Do not edit the settings:\n"
            << "# points computed with other cross sections, grid or tolerance would be invalid.\n"
            << "format: CRN-EEDF-table-1\n"
            << "settings:\n"
            << "  csDataFile: \"" << std::filesystem::path(m_set.csDataFile).filename().string() << "\"\n"
            << "  csDataFileSize: " << csFileSize() << "\n"
            << "  grid: {type: " << m_set.gridType << ", start: " << num(m_set.gridStart)
            << ", end: " << num(m_set.gridEnd) << ", points: " << m_set.gridPoints << "}\n"
            << "  convergence: {maxIterations: " << m_set.maxItr << ", rtol: " << num(m_set.rtol) << "}\n"
            << "points:\n";
    }

    // Parser of the flow-style map of one point, i.e. what formatPoint() writes: numbers and
    // maps of numbers, keys quoted or not. Much faster than the general YAML reader on big tables.
    struct PointParser {
        const std::string& s;
        size_t i = 0;

        void ws() { while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) i++; }
        void expect(char c)
        {
            ws();
            if (i >= s.size() || s[i] != c) throw std::runtime_error(std::string("expected '") + c + "'");
            i++;
        }
        bool accept(char c)
        {
            ws();
            if (i < s.size() && s[i] == c) {
                i++;
                return true;
            }
            return false;
        }
        std::string key()
        {
            ws();
            std::string k;
            if (accept('"')) {
                const size_t end = s.find('"', i);
                if (end == std::string::npos) throw std::runtime_error("unterminated key");
                k = s.substr(i, end - i);
                i = end + 1;
            } else {
                const size_t end = s.find(':', i);
                if (end == std::string::npos) throw std::runtime_error("missing ':'");
                k = s.substr(i, end - i);
                while (!k.empty() && std::isspace(static_cast<unsigned char>(k.back()))) k.pop_back();
                i = end;
            }
            expect(':');
            return k;
        }
        double number()
        {
            ws();
            char* end = nullptr;
            const double v = std::strtod(s.c_str() + i, &end);
            if (end == s.c_str() + i) throw std::runtime_error("expected a number");
            i = end - s.c_str();
            return v;
        }
        std::map<std::string, double> numberMap()
        {
            std::map<std::string, double> m;
            expect('{');
            if (accept('}')) return m;
            do {
                const std::string k = key();
                m[k] = number();
            } while (accept(','));
            expect('}');
            return m;
        }
        Point point()
        {
            Point p;
            expect('{');
            do {
                const std::string k = key();
                if (k == "X") p.X = numberMap();
                else if (k == "k") p.k = numberMap();
                else if (k == "EN") p.EN = number();
                else if (k == "T") p.T = number();
                else if (k == "Te") p.Te = number();
                else if (k == "power") p.power = number();
                else throw std::runtime_error("unknown key '" + k + "'");
            } while (accept(','));
            expect('}');
            return p;
        }
    };

    void load()
    {
        // Header (everything before "points:") with the YAML reader, points with PointParser
        std::ifstream in(m_path);
        std::string line, header;
        size_t lineNo = 0;
        while (std::getline(in, line)) {
            lineNo++;
            header += line + "\n";
            if (line.rfind("points:", 0) == 0) break;
        }
        Cantera::AnyMap root = Cantera::AnyMap::fromYamlString(header);
        auto fail = [&](const std::string& what) {
            throw std::runtime_error("EEDF table " + m_path + ": " + what
                                     + ". Use another eedfTable file or remove this one.");
        };
        if (!root.hasKey("format") || root["format"].asString() != "CRN-EEDF-table-1") {
            fail("unknown format");
        }
        const auto& set = root["settings"].as<Cantera::AnyMap>();
        const std::string csName = std::filesystem::path(m_set.csDataFile).filename().string();
        if (set["csDataFile"].asString() != csName
            || static_cast<std::uintmax_t>(set["csDataFileSize"].asInt()) != csFileSize()) {
            fail("computed with another cross-section file (" + set["csDataFile"].asString() + ")");
        }
        const auto& grid = set["grid"].as<Cantera::AnyMap>();
        if (grid["type"].asString() != m_set.gridType || toDouble(grid["start"]) != m_set.gridStart
            || toDouble(grid["end"]) != m_set.gridEnd || grid["points"].asInt() != m_set.gridPoints) {
            fail("computed with another Boltzmann grid");
        }
        const auto& conv = set["convergence"].as<Cantera::AnyMap>();
        if (conv["maxIterations"].asInt() != m_set.maxItr || toDouble(conv["rtol"]) != m_set.rtol) {
            fail("computed with other EEDF convergence settings (eedfRtol, eedfMaxIterations)");
        }

        if (!root.hasKey("points")) fail("missing 'points:'");
        while (std::getline(in, line)) {
            lineNo++;
            const size_t start = line.find_first_not_of(" \t");
            if (start == std::string::npos || line[start] == '#') continue;
            if (line[start] != '-') fail("line " + std::to_string(lineNo) + " is not a point ('- {...}')");
            try {
                PointParser parser{line, start + 1};
                m_points.push_back(parser.point());
                index(m_points.size() - 1);
            } catch (const std::exception& e) {
                fail("line " + std::to_string(lineNo) + ": " + e.what());
            }
        }
    }
};

} // namespace EEDF

#endif // EEDF_TABLE_H
