/*---------------------------------------------------------------------------*\
|   EEDFTable: builds the EEDF table of a ChemPlasKin CRN_v5 case in advance  |
|   Reactor Network module made by: Cristian Casalanguida                     |
|   Contact:      cristian.casalanguida@studenti.polito.it                    |
|-----------------------------------------------------------------------------|
|   Usage:                                                                    |
|     ./EEDFTable -case <dir> -T <list> [-EN <list>] [-comp <source>]...      |
|                 [-table <file>] [-tolT <K>] [-tolEN <Td>]                   |
|   <list>   comma separated values and ranges start:end:step,                |
|            e.g. -T 800:2500:1   -EN 0.1,20,40:60:2.5                        |
|   -EN      default: the E/N used by the case (E/N(NRP), E/N(DC), 0.1 Td)    |
|   -comp    case  : mixtures of the reactors in crnDict (default)            |
|            table : mixtures already present in the table                    |
|            <file>: one mixture per line, e.g. "H2:0.07, O2:0.20, N2:0.73"   |
|            can be repeated, the mixtures of all the sources are used        |
|   -table   default: eedfTable of chemPlasProperties or <case>/eedfTable.yaml|
|   -tolT, -tolEN  default: TGAS_TOLERANCE and EN_TOLERANCE of the case       |
\*---------------------------------------------------------------------------*/

#include "cantera/core.h"
#include "cantera/ext/bolos/Logger.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace Cantera; // required by utilities.h, as in masterCRN.cpp

#include "../../src/utilities.h"
#include "crnDict.h"
#include "eedfTable.h"

using Mixture = std::map<std::string, double>;

// Compute the EEDF point (EN, T, X) and add it to the table, unless the table already holds a
// point within tolerance. F0 is the initial guess and returns the last converged EEDF.
// Returns true when the point has been computed.
bool addPoint(EEDF::Table& table, CppBOLOS::BoltzmannSolver& bs, const EEDF::Settings& settings,
              double EN, double T, const Mixture& X, const std::vector<std::string>& processes,
              Eigen::VectorXd& F0)
{
    if (table.find(EN, T, X, processes)) {
        return false;
    }
    table.add(EEDF::solvePoint(bs, settings, EN, T, X, processes, F0, EN, 1E20 / (std::max(EN, 0.1) * std::max(EN, 0.1))));
    return true;
}

// Call addPoint for every combination of E/N, temperature and mixture. The temperatures are
// swept in increasing order so that each EEDF starts from the one of the previous temperature.
void buildTable(EEDF::Table& table, CppBOLOS::BoltzmannSolver& bs, const EEDF::Settings& settings,
                const std::vector<double>& ENs, std::vector<double> Ts, const std::vector<Mixture>& mixtures,
                const std::vector<std::string>& processes)
{
    std::sort(Ts.begin(), Ts.end());
    const size_t nTot = ENs.size() * Ts.size() * mixtures.size();
    size_t i = 0, nComputed = 0, nFailed = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (size_t iX = 0; iX < mixtures.size(); iX++) {
        for (const double EN : ENs) {
            Eigen::VectorXd F0;
            for (const double T : Ts) {
                i++;
                bool computed = false;
                try {
                    computed = addPoint(table, bs, settings, EN, T, mixtures[iX], processes, F0);
                } catch (const std::exception& e) {
                    std::cerr << "\n[Warning] EN = " << EN << " Td, T = " << T << " K, mixture " << iX
                              << ": " << e.what() << "\n";
                    F0.resize(0);
                    nFailed++;
                    continue;
                }
                if (computed) nComputed++;
                const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                std::cout << "\r[" << i << "/" << nTot << "] mixture " << iX << ", EN = " << EN << " Td, T = " << T
                          << " K: " << (computed ? "computed" : "in table") << " | computed " << nComputed
                          << ", " << static_cast<long>(elapsed) << " s      " << std::flush;
            }
        }
    }
    std::cout << "\n" << nComputed << " points computed, " << (nTot - nComputed - nFailed)
              << " already in the table, " << nFailed << " failed. Table size: " << table.size() << std::endl;
}

// "a,b,start:end:step,..." -> values
std::vector<double> parseList(const std::string& s)
{
    std::vector<double> values;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) continue;
        std::vector<double> r;
        std::stringstream is(item);
        std::string v;
        while (std::getline(is, v, ':')) r.push_back(std::stod(v));
        if (r.size() == 1) {
            values.push_back(r[0]);
        } else if (r.size() == 3 && r[2] > 0 && r[1] >= r[0]) {
            const int n = static_cast<int>(std::floor((r[1] - r[0]) / r[2] + 1e-9));
            for (int j = 0; j <= n; j++) values.push_back(r[0] + j * r[2]);
        } else {
            throw std::runtime_error("Invalid list item '" + item + "' (use value or start:end:step)");
        }
    }
    return values;
}

// Mixture restricted to the Boltzmann species, missing species = 0 (as masterCRN.cpp)
Mixture boltzmannMixture(const Mixture& composition, const std::vector<std::string>& bSpecies)
{
    Mixture X;
    for (const auto& sp : bSpecies) {
        X[sp] = composition.count(sp) ? composition.at(sp) : 0.0;
    }
    return X;
}

// One mixture per line: "SP1:x1, SP2:x2, ..."; '#' starts a comment
std::vector<Mixture> readMixtureFile(const std::string& fileName)
{
    std::ifstream in(fileName);
    if (!in.is_open()) throw std::runtime_error("Unable to open mixture file: " + fileName);
    std::vector<Mixture> mixtures;
    std::string line;
    while (std::getline(in, line)) {
        line = line.substr(0, line.find('#'));
        Mixture m;
        std::stringstream ss(line);
        std::string pair;
        while (std::getline(ss, pair, ',')) {
            const auto colon = pair.find(':');
            if (colon == std::string::npos) continue;
            std::string sp = pair.substr(0, colon);
            sp.erase(std::remove_if(sp.begin(), sp.end(), ::isspace), sp.end());
            m[sp] = std::stod(pair.substr(colon + 1));
        }
        if (!m.empty()) mixtures.push_back(m);
    }
    return mixtures;
}

// Distinct process strings of the Boltzmann reactions of a mechanism
void addMechanismProcesses(const std::string& mechFile, std::vector<std::string>& processes)
{
    Cantera::AnyMap root = Cantera::AnyMap::fromYamlFile(mechFile);
    for (const auto& r : root["reactions"].asVector<Cantera::AnyMap>()) {
        if (r.hasKey("type") && r["type"].asString() == "Boltzmann" && r.hasKey("process")) {
            const std::string& process = r["process"].asString();
            if (std::find(processes.begin(), processes.end(), process) == processes.end()) {
                processes.push_back(process);
            }
        }
    }
}

int main(int argc, char* argv[])
{
    CppBOLOS::currentLogLevel = CppBOLOS::LOG_WARNING;

    std::string casePath = "..", tablePath, TArg, ENArg;
    std::vector<std::string> compSources;
    double tolT = -1.0, tolEN = -1.0;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value after " + arg);
            return argv[++i];
        };
        if (arg == "-case") casePath = next();
        else if (arg == "-table") tablePath = next();
        else if (arg == "-T") TArg = next();
        else if (arg == "-EN") ENArg = next();
        else if (arg == "-comp") compSources.push_back(next());
        else if (arg == "-tolT") tolT = std::stod(next());
        else if (arg == "-tolEN") tolEN = std::stod(next());
        else {
            std::cerr << "Unknown option " << arg << " (see the header of eedfTableBuilder.cpp)\n";
            return 1;
        }
    }
    if (TArg.empty()) {
        std::cerr << "Usage: EEDFTable -case <dir> -T <list> [-EN <list>] [-comp case|table|<file>]..."
                  << " [-table <file>] [-tolT <K>] [-tolEN <Td>]\n";
        return 1;
    }

    const std::string parameterPath = casePath + "/chemPlasProperties";
    const std::string crnDictPath = casePath + "/crnDict";

    // Solver settings and tolerances: the same values ChemPlasKin reads from the case
    EEDF::Settings settings;
    settings.csDataFile = readParameter<std::string>(parameterPath, "csDataFile");
    settings.gridType = readParameter<std::string>(parameterPath, "gridType");
    auto gridSize = readParameter<std::map<std::string, double>>(parameterPath, "gridSize");
    settings.gridStart = gridSize["start"];
    settings.gridEnd = gridSize["end"];
    settings.gridPoints = (int)gridSize["points"];
    settings.rtol = readParameter<double>(parameterPath, "eedfRtol", 1e-5);
    settings.maxItr = (int)readParameter<double>(parameterPath, "eedfMaxIterations", 200.0);

    EEDF::Tolerance tolerance;
    tolerance.T = (tolT >= 0) ? tolT : readParameter<double>(parameterPath, "TGAS_TOLERANCE");
    tolerance.EN = (tolEN >= 0) ? tolEN : readParameter<double>(parameterPath, "EN_TOLERANCE");
    tolerance.maxGap = readParameter<double>(parameterPath, "eedfTableMaxGap", tolerance.maxGap);
    tolerance.Xrel = readParameter<double>(parameterPath, "eedfTableXrel", tolerance.Xrel);
    tolerance.Xabs = readParameter<double>(parameterPath, "eedfTableXabs", tolerance.Xabs);

    if (tablePath.empty()) {
        tablePath = readParameter<std::string>(parameterPath, "eedfTable", std::string("<case>/eedfTable.yaml"));
        if (tablePath.rfind("<case>", 0) == 0) tablePath.replace(0, 6, casePath);
    }
    if (tablePath == "none") {
        std::cerr << "eedfTable is 'none' in chemPlasProperties: give a file with -table\n";
        return 1;
    }

    // E/N values: the ones requested by ChemPlasKin unless given
    std::vector<double> ENs;
    if (ENArg.empty()) {
        const double EN_NRP = readParameter<double>(parameterPath, "E/N(NRP)");
        const double EN_DC = readParameter<double>(parameterPath, "E/N(DC)");
        ENs = {EN_NRP, std::max(EN_DC, 0.1), 0.1};
    } else {
        ENs = parseList(ENArg);
    }
    std::sort(ENs.begin(), ENs.end());
    ENs.erase(std::unique(ENs.begin(), ENs.end()), ENs.end());
    const std::vector<double> Ts = parseList(TArg);

    // Boltzmann processes of all the mechanisms of the case
    const auto bSpecies = readParameter<std::vector<std::string>>(parameterPath, "BoltzmannSpecies");
    const std::string mechFile = readParameter<std::string>(parameterPath, "mechFile");
    auto reactorZones = readReactors(crnDictPath);
    std::vector<std::string> processes;
    addMechanismProcesses(mechFile, processes);
    for (const auto& [name, zone] : reactorZones) {
        if (zone.mechanismFile.has_value()) addMechanismProcesses(zone.mechanismFile.value(), processes);
    }

    EEDF::Table table(tablePath, settings, tolerance);

    // Mixtures
    if (compSources.empty()) compSources.push_back("case");
    std::vector<Mixture> mixtures;
    auto addMixture = [&](const Mixture& X) {
        if (std::find(mixtures.begin(), mixtures.end(), X) == mixtures.end()) mixtures.push_back(X);
    };
    for (const auto& source : compSources) {
        if (source == "case") {
            const auto globalMixture = readParameter<Mixture>(parameterPath, "Mixture");
            for (const auto& [name, zone] : reactorZones) {
                if (zone.type == "reservoir") continue;
                addMixture(boltzmannMixture(zone.composition.value_or(globalMixture), bSpecies));
            }
        } else if (source == "table") {
            for (const auto& p : table.points()) addMixture(p.X);
        } else {
            for (const auto& m : readMixtureFile(source)) addMixture(boltzmannMixture(m, bSpecies));
        }
    }

    std::cout << "EEDF table: " << tablePath << " (" << table.nLoaded() << " points loaded)\n"
              << "Processes: " << processes.size() << " | E/N [Td]:";
    for (const double EN : ENs) std::cout << " " << EN;
    std::cout << " | T: " << Ts.size() << " values in [" << *std::min_element(Ts.begin(), Ts.end()) << ", "
              << *std::max_element(Ts.begin(), Ts.end()) << "] K | mixtures: " << mixtures.size() << "\n"
              << "Tolerances: T " << tolerance.T << " K (point already present), E/N " << tolerance.EN
              << " Td, interpolation gap " << tolerance.maxGap << " K\n";
    for (size_t j = 1; j < Ts.size(); j++) {
        if (std::abs(Ts[j] - Ts[j - 1]) > tolerance.maxGap) {
            std::cout << "[Note] T spacing larger than eedfTableMaxGap: temperatures between the points "
                         "will be computed during the simulation.\n";
            break;
        }
    }

    CppBOLOS::BoltzmannSolver bs;
    EEDF::setupSolver(bs, settings);
    buildTable(table, bs, settings, ENs, Ts, mixtures, processes);
    return 0;
}
