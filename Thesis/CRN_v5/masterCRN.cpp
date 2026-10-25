/*--------------------------------*- C++ -*----------------------------------*\
|     ____ _                    ____  _           _  ___                      |
|    / ___| |__   ___ _ __ ___ |  _ \| | __ _ ___| |/ (_)_ __                 |
|   | |   | '_ \ / _ \ '_ ` _ \| |_) | |/ _` / __| ' /| | '_ \                |
|   | |___| | | |  __/ | | | | |  __/| | (_| \__ \ . \| | | | |               |
|    \____|_| |_|\___|_| |_| |_|_|   |_|\__,_|___/_|\_\_|_| |_|               |
|                                                                             |
|   A Freeware for Unified Gas-Plasma Kinetics Simulation                     |
|   Version:      1.2 (February 2026)                                         |
|   License:      GNU LESSER GENERAL PUBLIC LICENSE, Version 2.1              |
|   Author:       Xiao Shao                                                   |
|   Organization: King Abdullah University of Science and Technology (KAUST)  |
|   Contact:      xiao.shao@kaust.edu.sa                                      |
|-----------------------------------------------------------------------------|
|   Reactor Network module made by: Cristian Casalanguida                     |
|   Contact:      cristian.casalanguida@studenti.polito.it                    |
|   Last edit:    2026-10-24 (CRN_v5: EEDF table, see eedfTable.h)           |
|                 CRN_v4: per-reactor Boltzmann state, fewer CVODE restarts,  |
|                 flowUpdateInterval, writeInterval                           |
\*---------------------------------------------------------------------------*/

// Python.h must come before any system headers that may conflict
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "cantera/core.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/ext/bolos/Logger.h"
#include "cantera/kinetics/Boltzmann.h"
#include <iostream>
#include "cantera/zerodim.h"

#include "cantera/numerics/Integrator.h"
#include <cmath>
#include <optional>
#include <filesystem>

using namespace Cantera;

static bool g_verbose = false;

#include  "plasmaReactor.h" // for original one "../../src/plasmaReactor.h"
#include "../../src/utilities.h"
#include <unordered_map>
#include "crnDict.h"
#include "eedfTable.h"


// Struct for mass flow controller entries
struct MassFlowControllerDef {
    std::string start;
    std::string end;
    double value = 0.0;
    double T = 0.0;
    std::vector<double> comp;
    double h = 0.0; // specific enthalpy of source [J/kg]
    // CRN_v4: source->end species index map and time integrals of the source state,
    // used to pass the time-averaged outflow between two flow updates
    std::vector<size_t> srcIdx;
    std::vector<double> sumYdt;
    double sumHdt = 0.0;
    double sumTdt = 0.0;
    double sumDt  = 0.0;
};

struct ReactorState {
    std::string name;
    std::shared_ptr<std::map<std::string, double>> boltzmannSpecies;
    std::shared_ptr<Solution> sol;
    std::shared_ptr<ThermoPhase> gas;
    std::unique_ptr<ChemPlasReactor> odes;
    std::unique_ptr<Integrator> integrator;
    std::ofstream outputFile;
    std::vector<int> indexList;
    std::vector<double> qf;
    std::vector<double> qr;
    std::vector<double> q;

    // CRN_v4: Boltzmann state owned by this reactor. The global BoltzmannRate::bsolver
    // is only used to (re)compute the EEDF; between updates the reactor's rate
    // constants are installed via BoltzmannRate::cachedSolutions.
    Eigen::VectorXd F0;          // last EEDF computed for this reactor, initial guess of the next one
    std::vector<double> bCache;  // [k_1 ... k_NumProcess, Te, E/N] as in BoltzmannRate::cacheBoltzmannSolutions
    double EN_b = 0.0;           // E/N of the current EEDF [Td]
    double kT_b = 0.0;           // gas temperature of the current EEDF [K]
    double TeSolver = 0.0;       // electron temperature from CppBOLOS [K]
    double elecPower = 0.0;      // power/N absorbed by electrons [eV m^3/s]
    bool needsReinit = true;     // CVODE must be restarted before the next integration
    double mdot_in = 0.0;        // [kg/s], updated with the flows
    double mdot_out = 0.0;       // [kg/s], updated with the flows
    std::ofstream pulseFile;     // one row per pulse (discharge reactors only)
    double neMaxPulse = 0.0;     // max electron density in the current pulse [#/cm^3]
    // End of the discharge at Ep: internal time of CVODE when the reactor is advanced with
    // step() (-1 = unknown, after integrate()), and pulse end already found by the root search
    double tInternal = -1.0;
    bool pulseEndFound = false;
    double EpDischarge = 0.0;    // energy deposited at the end of the discharge [mJ/cm^3]
};

struct ReservoirState {
    std::string name;
    std::shared_ptr<std::map<std::string, double>> boltzmannSpecies;
    std::shared_ptr<Solution> sol;
    std::shared_ptr<ThermoPhase> gas;
    std::ofstream outputFile;
    std::vector<int> indexList;
};

// Parse massFlowControllers block from OpenFOAM-style dict
std::map<std::string, MassFlowControllerDef> readMassFlowControllers(const std::string& fileName) {
    std::ifstream inFile(fileName);
    if (!inFile.is_open()) {
        throw std::runtime_error("Unable to open file: " + fileName);
    }

    auto normalizeToken = [](std::string token) {
        const auto commentPos = token.find("//");
        if (commentPos != std::string::npos) {
            token = token.substr(0, commentPos);
        }
        while (!token.empty() && (token.back() == ';' || token.back() == ',')) {
            token.pop_back();
        }
        return token;
    };

    std::map<std::string, MassFlowControllerDef> controllers;
    std::string line;
    bool inMFCBlock = false;

    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key)) continue;

        // Find massFlowControllers block
        if (!inMFCBlock && key == "massFlowControllers") {
            inMFCBlock = true;
            continue;
        }

        if (!inMFCBlock) continue;

        // End of massFlowControllers block
        if (line.find('}') != std::string::npos && inMFCBlock) {
            int braceCount = 0;
            for (char c : line) {
                if (c == '{') braceCount++;
                if (c == '}') braceCount--;
            }
            if (braceCount < 0) break;
        }

        // Read controller name
        if (key != "{" && key != "}") {
            std::string ctrlName = key;
            MassFlowControllerDef spec;

            // Find opening brace
            if (line.find('{') == std::string::npos) {
                while (std::getline(inFile, line) && line.find('{') == std::string::npos) {}
            }

            // Read controller definition
            while (std::getline(inFile, line)) {
                if (line.find('}') != std::string::npos) break;

                std::istringstream ls(line);
                std::string field;
                if (!(ls >> field)) continue;
                field = normalizeToken(field);

                if (field == "start") {
                    std::string val;
                    if (ls >> val) spec.start = normalizeToken(val);
                } else if (field == "end") {
                    std::string val;
                    if (ls >> val) spec.end = normalizeToken(val);
                } else if (field == "value") {
                    double val;
                    if (ls >> val) spec.value = val;
                }
            }
            controllers[ctrlName] = std::move(spec);
        }
    }

    return controllers;
}

// Tools for creating and advancing reactor states, used in main loop for better readability

ReactorState createReactorState(
    const std::string& name,
    const ReactorZoneDef& config,
    const std::string& plasmaHeatModel,
    const std::string& parameterPath,
    const std::string& controlDictPath,
    double EN,
    double runTime)
{
    ReactorState state;
    state.name = name;
    state.boltzmannSpecies = std::make_shared<std::map<std::string, double>>(config.boltzmannSpecies.value());
    const auto& boltzmannSpecies = *state.boltzmannSpecies;
    auto composition = *config.composition;
    const double initialTemperature = *config.initialTemperature;
    const double pressure = *config.pressure;
    const std::string outputPath = config.outputFile.value_or("../output/" + state.name + ".csv");

    BoltzmannRate::bsolver.set_kT (initialTemperature);
    BoltzmannRate::bsolver.set_EN(EN);
    BoltzmannRate::bsolver.set_density(boltzmannSpecies);
    BoltzmannRate::bsolver.init();

    // Solve EBE, serve as cache and preliminary convergence check
    // (BoltzmannRate objects created by newSolution need a valid global F0)
    BoltzmannRate::F0 = BoltzmannRate::bsolver.maxwell(2.0); // initial guess from Maxwell EEDF
    BoltzmannRate::F0 = BoltzmannRate::bsolver.converge(BoltzmannRate::F0, 200, 1e-5);
    //

    state.sol = newSolution(config.mechanismFile.value(), "", "mixture-averaged");

    state.gas = state.sol->thermo();

    state.odes = std::make_unique<ChemPlasReactor>(state.sol, *state.boltzmannSpecies, plasmaHeatModel);

    // User can assign electron mole fraction in 'Mixture' to override 'electronDens'
    const double gasNumberDensity = 1e-6 * pressure / initialTemperature / CppBOLOS::KB;
    if (composition.count(state.gas->speciesName(state.odes->electronIndex)) == 0) {
        composition[state.gas->speciesName(state.odes->electronIndex)] = config.initialElectronDensity.value() / gasNumberDensity;
    }
    state.gas->setState_TPX(initialTemperature, pressure, composition);
    state.odes->setConstPD(state.gas->pressure(), state.gas->density()); // Record pressure and density

    state.odes->nonThermal = readParameter<bool>(parameterPath, "nonThermal");
    state.odes->heatLoss = readParameter<bool>(parameterPath, "heatLoss");
    if (state.odes->heatLoss) {
        auto C0 = readParameter<double>(parameterPath, "C0");
        state.odes->T0 = initialTemperature;
        state.odes->C0 = C0;
        if (g_verbose) std::cout << "Start with constant-volume PAC reactor since heat loss model is on."<< std::endl;
        state.odes->constPressure = false;
    } else
    {
        if (readParameter<string>(parameterPath, "CP_or_CV") == "CP"){
            if (g_verbose) std::cout << "Start with constant-pressure PAC reactor"<< std::endl;
            state.odes->constPressure = true;
        } else {
            if (g_verbose) std::cout << "Start with constant-volume PAC reactor"<< std::endl;
            state.odes->constPressure = false;
        }
    }
    state.odes->inertSpIndex = state.odes->findSpeciesIndex(readParameter<string>(parameterPath, "inertSpecies"));

    // Create and initialize the ODE integrator

    state.integrator = std::unique_ptr<Integrator>(newIntegrator("CVODE"));

    // Tolerances must be set before initialize(): Cantera passes them to CVODE only there
    // (reinitialize() does not). Before 2026-10-25 they were set after it and never used, so
    // CVODE ran with Cantera's defaults (reltol 1e-9, abstol 1e-15), still used if 'odes' lacks them.
    auto Tolerances = readParameter<std::map<std::string, double>>(controlDictPath, "odes",
                                                                   std::map<std::string, double>{});
    const double reltol = Tolerances.count("reltol") ? Tolerances["reltol"] : 1e-9;
    const double abstol = Tolerances.count("abstol") ? Tolerances["abstol"] : 1e-15;
    state.integrator->setTolerances(reltol, abstol);
    state.integrator->initialize(runTime, *state.odes);

    auto kin = state.sol->kinetics();
    const int irxns = kin->nReactions();
    state.qf.resize(irxns);
    state.qr.resize(irxns);
    state.q.resize(irxns);

    state.outputFile.open(outputPath);
    if (!state.outputFile.is_open()) {
        std::filesystem::create_directories(std::filesystem::path(outputPath).parent_path());
        state.outputFile.open(outputPath);
    }
    if (!state.outputFile.is_open()) {
        throw std::runtime_error("Failed to open output file for reactor: " + name);
    }

    const std::vector<std::string> speciesNames = state.gas->speciesNames();
    state.outputFile << "Time(s), T_gas(K), p(Pa), N_gas(#/cm^3), MW(kg/kmol)";
    if (g_verbose) std::cout << "Writing information of " << name << ":\n" << "Time, T_gas(K), p(Pa), N_gas(#/cm^3), MW(kg/kmol)";
    for (const auto& sp : speciesNames) {
        if (g_verbose) std::cout << ", " << sp;
        const size_t index = state.gas->speciesIndex(sp);
        if (index < std::numeric_limits<size_t>::max()) {
            state.indexList.push_back(static_cast<int>(index));
        } else {
            throw std::runtime_error("No valid species name found: " + sp);
        }
        state.outputFile << ", " << sp;
    }
    state.outputFile << ", Power_input(W/cm^3), Energy_deposited(mJ/cm^3)";
    state.outputFile << std::endl;
    if (g_verbose) std::cout << std::endl;

    // Write initial values
    state.outputFile << std::setprecision(15) << runTime << std::setprecision(6) << ", " << state.gas->temperature() << ", " << state.gas->pressure() << ", " << 1e-6 * Avogadro * state.gas->molarDensity() << ", " << state.gas->meanMolecularWeight();
    if (g_verbose) std::cout << runTime << ", " << state.gas->temperature() << ", " << state.gas->pressure() << ", " << 1e-6 * Avogadro * state.gas->molarDensity() << ", " << state.gas->meanMolecularWeight();
    for (const auto index : state.indexList) {
        double dens = getNumberDens(state.gas, static_cast<size_t>(index));
        state.outputFile << ", " << dens;
        if (g_verbose) std::cout << ", " << dens;
    }
    state.outputFile << ", " << 0.0 << ", " << 0.0;
    state.outputFile << std::endl;
    if (g_verbose) std::cout << std::endl;

    if (config.discharge) {
        std::string pulsePath = outputPath;
        const auto dot = pulsePath.rfind(".csv");
        pulsePath = (dot != std::string::npos ? pulsePath.substr(0, dot) : pulsePath) + "_pulses.csv";
        state.pulseFile.open(pulsePath);
        if (!state.pulseFile.is_open()) {
            throw std::runtime_error("Failed to open pulse file for reactor: " + name);
        }
        state.pulseFile << "iPulse, Time(s), T_gas(K), p(Pa), n_e_end(#/cm^3), n_e_max(#/cm^3), "
                        << "Energy_deposited(mJ/cm^3), tau_NRP(s), Energy_discharge(mJ/cm^3)" << std::endl;
    }

    return state;
}

// CRN_v5: Boltzmann outputs of reactor s (Te, power/N, rate constants) for E/N, T and the
// current CppBOLOS mixture of s. They are interpolated from the EEDF table when it holds suitable
// points; otherwise the EEDF is computed as in CRN_v4 (initial guess, EN_ref, delta0) and the point
// is added to the table.
static size_t g_eedfFromTable = 0, g_eedfComputed = 0;

void computeEEDF(ReactorState& s, double EN_new, double T_new, double EN_ref, double delta0,
                 EEDF::Table& table, const EEDF::Settings& eedfSettings,
                 const std::vector<std::string>& processes)
{
    std::optional<EEDF::Point> p = table.interpolate(EN_new, T_new, *s.boltzmannSpecies, processes);
    if (p) {
        g_eedfFromTable++;
    } else {
        p = table.add(EEDF::solvePoint(BoltzmannRate::bsolver, eedfSettings, EN_new, T_new,
                                       *s.boltzmannSpecies, processes, s.F0, EN_ref, delta0));
        BoltzmannRate::F0 = s.F0;
        g_eedfComputed++;
    }

    s.EN_b = p->EN;
    s.kT_b = p->T;
    s.TeSolver = p->Te;
    s.elecPower = p->power;
    s.odes->setBoltzmannOutputs(s.elecPower, s.TeSolver);

    // Rate constants of all Boltzmann processes (indexed by the global processIndex)
    const size_t N = BoltzmannRate::NumProcess;
    s.bCache.assign(N + 2, 0.0);
    for (size_t i = 0; i < BoltzmannRate::ProcessesList.size(); i++) {
        s.bCache[i] = p->k.at(BoltzmannRate::ProcessesList[i]);
    }
    s.bCache[N] = std::max(p->Te, p->T); // enforce T_e >= T_gas, as BoltzmannRate::updateBoltzmannSolver
    s.bCache[N + 1] = p->EN;
}

// CRN_v4: install the Boltzmann rate constants of reactor s in the global cache
// read by BoltzmannRate::evalFromStruct (useCachedSolutions = true)
void activateReactor(const ReactorState& s)
{
    const size_t N = BoltzmannRate::NumProcess;
    BoltzmannRate::cachedSolutions = s.bCache;
    BoltzmannData::electronTemp = s.bCache[N];
    BoltzmannData::reducedEfield = s.bCache[N + 1];
    BoltzmannRate::updateOnce = true; // MultiRate caches re-read the rate constants
}

// Advance a discharging reactor from t0 to t1 with CVODE internal steps (step()), watching
// g = EpTarget - Ep (Ep = deposited plasma energy, y[2]). If g changes sign or becomes zero, the
// crossing time t* is located on the CVODE interpolant of the internal step (derivative(t, 0),
// no extra RHS evaluations), the reactor state is set to y(t*) and t* is returned. Otherwise the
// state is set to y(t1) and -1 is returned. CVODE must not be behind t0 (s.tInternal >= t0).
double integrateDischarge(ReactorState& s, double t0, double t1, double EpTarget)
{
    auto& integ = *s.integrator;
    double a = t0;
    double Ea = s.odes->depositedPlasmaEnergy(); // energy at t0, below the target
    double tn = s.tInternal;
    while (true) {
        // part of [t0, t1] covered by the last internal step: [a, b]
        const double b = std::min(tn, t1);
        if (b > a) {
            const double Eb = integ.derivative(b, 0)[2];
            if (Eb >= EpTarget) {
                // bracketed root: secant with bisection safeguard, down to 1e-6 of the target
                // or to adjacent representable times
                double lo = a, hi = b, Elo = Ea, Ehi = Eb;
                while (std::abs(Ehi - EpTarget) > 1e-6 * EpTarget && std::nextafter(lo, hi) < hi) {
                    double tm = lo + (hi - lo) * (EpTarget - Elo) / (Ehi - Elo);
                    tm = std::min(std::max(tm, lo + 0.01 * (hi - lo)), hi - 0.01 * (hi - lo));
                    if (!(tm > lo && tm < hi)) tm = lo + 0.5 * (hi - lo);
                    if (!(tm > lo && tm < hi)) break;
                    const double Em = integ.derivative(tm, 0)[2];
                    if (Em < EpTarget) {
                        lo = tm;
                        Elo = Em;
                    } else {
                        hi = tm;
                        Ehi = Em;
                    }
                }
                s.odes->updateState(integ.derivative(hi, 0));
                s.tInternal = -1.0;
                return hi;
            }
            a = b;
            Ea = Eb;
        }
        if (tn >= t1) break;
        tn = integ.step(t1);
    }
    s.odes->updateState(integ.derivative(t1, 0));
    s.tInternal = tn; // CVODE is at tn >= t1
    return -1.0;
}

// Write one row of the reactor output file
void writeReactorRow(ReactorState& state, double runTime, double power, double energy)
{
    // time with 15 significant digits (6 would give a 10 us resolution at t ~ 3 s, dt can be 1e-14 s)
    state.outputFile << std::setprecision(15) << runTime << std::setprecision(6) << ", " << state.gas->temperature() << ", " << state.gas->pressure() << ", "
                     << 1e-6 * Avogadro * state.gas->molarDensity() << ", " << state.gas->meanMolecularWeight();
    for (const auto index : state.indexList) {
        state.outputFile << ", " << getNumberDens(state.gas, static_cast<size_t>(index));
    }
    state.outputFile << ", " << power << ", " << energy << "\n";
}

ReservoirState createReservoirState(
    const std::string& name,
    const ReactorZoneDef& config,
    double runTime)
{
    ReservoirState state;
    state.name = name;
    state.boltzmannSpecies = std::make_shared<std::map<std::string, double>>(config.boltzmannSpecies.value());

    state.sol = newSolution(config.mechanismFile.value(), "", "mixture-averaged");
    state.gas = state.sol->thermo();
    state.gas->setState_TPX(config.initialTemperature.value(), config.pressure.value(), *config.composition);

    if (g_verbose) std::cout << "Reservoir '" << name << "': T=" << config.initialTemperature.value()
              << " K, p=" << config.pressure.value() << " Pa, M=" << state.gas->meanMolecularWeight() << " kg/kmol\n";

    const std::string outputPath = config.outputFile.value_or("../output/" + name + ".csv");
    state.outputFile.open(outputPath);
    if (!state.outputFile.is_open()) {
        std::filesystem::create_directories(std::filesystem::path(outputPath).parent_path());
        state.outputFile.open(outputPath);
    }
    if (!state.outputFile.is_open())
        throw std::runtime_error("Failed to open output file for reservoir: " + name);

    const std::vector<std::string>& speciesNames = state.gas->speciesNames();
    state.outputFile << "Time(s), T_gas(K), p(Pa), N_gas(#/cm^3), MW(kg/kmol)";
    for (const auto& sp : speciesNames) {
        size_t idx = state.gas->speciesIndex(sp);
        state.indexList.push_back(static_cast<int>(idx));
        state.outputFile << ", " << sp;
    }
    state.outputFile << std::endl;
    // Write initial reservoir state
    state.outputFile << std::setprecision(15) << runTime << std::setprecision(6) << ", " << state.gas->temperature() << ", "
                     << state.gas->pressure() << ", "
                     << 1e-6 * Avogadro * state.gas->molarDensity() << ", "
                     << state.gas->meanMolecularWeight();
    for (const int idx : state.indexList)
        state.outputFile << ", " << getNumberDens(state.gas, static_cast<size_t>(idx));
    state.outputFile << std::endl;

    return state;
}

void findTimeStep(
    ReactorState& state, ReactorZoneDef& config,
    double runTime,
    double startTime,
    double SMALL,
    double pulsePeriod, // TO VECTORIZE IF THE REACTORS ARE INDEPENDENT
    double pulseWidth, // TO VECTORIZE IF THE REACTORS ARE INDEPENDENT
    double pulseExten, // TO VECTORIZE IF THE REACTORS ARE INDEPENDENT
    uint8_t& dischargeOn, uint8_t& fastExpansion, 
    bool stayEN_DC, // TO VECTORIZE IF THE REACTORS ARE INDEPENDENT
    const std::string& CP_or_CV, // TO VECTORIZE IF THE REACTORS ARE INDEPENDENT
    int iPulse, int nPulses, // TO VECTORIZE IF THE REACTORS ARE INDEPENDENT
    double& dt,
    double dt1,
    double dt_max,
    double& dT,
    double dT_max,
    double& dTdt,
    double& nE_1, double& nE_2, double& nE_3,
    double& tau_NRP, double& disEnergy, 
    int& iStepInPulse, int NstepInPulse,
    int& iStepAfPulse, int NstepAfPulse,
    double EN, double EN_DC, // TO VECTORIZE IF THE REACTORS ARE INDEPENDENT
    double& EN_temp, const double EN_TOLERANCE, const double TGAS_TOLERANCE)
{
    // Function for dynamic time stepping

    // Ignition detection and time step adjustment based on dT/dt
    double dt_ig = dt_max; // ignition delay time, to be calculated based on dT/dt
    if (dT > dT_max && dTdt > 1000) {
        dt_ig = std::min(dT_max/dTdt, dt_max); // reduce time step to capture ignition
        std::cout << "\n!!!!!!!!!!IGNITION DETECTED!!!!!!!!!!" << std::endl;
    }

    EN_temp = std::max(EN_DC, 0.1);  // default E/N value

    // Time step adjustment based on pulse
    if (iPulse < nPulses && config.discharge) { // Still pulsing
        double dt3 = startTime + (iPulse + 1) * pulsePeriod - runTime; // Avoid overstepping
        dt_ig = std::min(dt3, dt_ig);

        if ((runTime + SMALL) < (startTime + iPulse * pulsePeriod + pulseExten)){ // During pulse or pulse extension
            if (g_verbose) std::cout << " Reactor " << state.name << " discharging: " << (config.discharge ? "ON" : "OFF") << "; runTime = " << runTime << "; iPulse = " << iPulse << std::endl;
            if ((runTime + SMALL) < (startTime + iPulse * pulsePeriod + pulseWidth) && dischargeOn) { // During discharge
                if (iPulse <= 1 || nE_1 < 1e10) {
                    if (g_verbose) std::cout << "Using fixed time step for first two pulses." << std::endl;
                    dt = dt1;
                } else {
                    if (g_verbose) std::cout << "Using log-dynamic time step." << std::endl;
                    double K = std::max(std::round(log10(nE_2 / nE_1) * 100) / 100, 0.2);
                    if (g_verbose) std::cout << "Growth order K = " << K << "; tau_NRP = " << tau_NRP << std::endl;
                    dt = logDynamicTimestep(K, tau_NRP, NstepInPulse, iStepInPulse);
                }

                EN_temp = EN; // Should be an output to update Boltzmann solver

                // FOR THE VECTORIZATION OF THE PULSES
                iStepInPulse += 1;
            } else { // After discharge, during pulse extension
                if (iStepAfPulse == 0) { // First step after the discharge is off; record it
                    // Already done at the crossing time when the discharge ended at Ep (pulseEndFound)
                    if (!state.pulseEndFound) {
                        // sim.reinitialize();
                        state.integrator->reinitialize(runTime, *state.odes);
                        state.tInternal = runTime;

                        // FOR THE VECTORIZATION OF THE PULSES
                        tau_NRP = (runTime - startTime - iPulse * pulsePeriod); // Duration of the nanosecond discharge
                        if (g_verbose) std::cout << "NSD terminated at " << disEnergy << " [mJ/cm^3]" << ", tau_NRP = " << tau_NRP << "\n";

                        nE_2 = getNumberDens(state.gas, state.odes->electronIndex); // Store n_e at the end of a nanosecond pulse
                        state.EpDischarge = 1e-3 * state.odes->depositedPlasmaEnergy();
                    }

                        // TODO: may need to update Boltzmann grid if E/N changes significantly

                        fastExpansion = true;
                    }
                double K = (iPulse <= 1 || nE_2 < 1e10) ? -1 : std::min(round(log10(nE_3 / nE_2) * 100) / 100, -0.1); // fixed K for first two pulses
                if (g_verbose) std::cout << "Decay order K = " << K << std::endl;
                dt = (iStepAfPulse >= NstepAfPulse) ? // Make dt unchanged for iStepAfPulse >= NstepAfPulse to cut off dt growth
                    dt : logDynamicTimestep(K, pulseExten, NstepAfPulse, iStepAfPulse);
                
                // Isentropic heat loss model. P, T, rho jump at 100 ns after the pulse
                if (runTime > (startTime + iPulse * pulsePeriod + 1e-7) && state.odes->heatLoss && fastExpansion)
                {
                    double T_new = state.gas->temperature() * pow( config.pressure.value()/state.gas->pressure(), 1 - state.gas->cv_mass()/state.gas->cp_mass() );
                    if (g_verbose) std::cout << "Temperature jump: " << state.gas->temperature() << " -> " << T_new << "\n";
                    state.gas->setState_TPY(T_new, config.pressure.value(), state.gas->massFractions());
                    state.odes->setConstPD(state.gas->pressure(), state.gas->density());
                    state.needsReinit = true; // state changed outside CVODE
                    if (CP_or_CV == "CP"){
                        state.odes->constPressure = true;
                    } else {
                        state.odes->constPressure = false; // Not necessary since heatLoss already changes the state, but just to be safe :)
                    }
                    fastExpansion = false; // only do once for a pulse
                }
                
                // FOR THE VECTORIZATION OF THE PULSES
                iStepAfPulse += 1;
                if (iStepAfPulse == NstepAfPulse) { // At the end of pulseExten, record n_e
                    nE_3 = getNumberDens(state.gas, state.odes->electronIndex);
                }
            }
        } else { // Between pulses, after pulse extension
            dt = std::min(dt3, dt_max);
        }
    } else { // No more pulsing, just use regular adaptive time stepping based on dT/dt
        if (not stayEN_DC) EN_temp = 0.1;
        dt = dt_max;
    }

    dt = std::min(dt, dt_ig); // Always ensure dt is small enough to capture ignition
    dt = std::max(dt, 1e-14);  // Make sure dt > 1e-14, avoid dt cascade
    
    if (dt < 1e-16) {
            std::cerr << "dt is too small !" << "\n";
            abort();
        }
    
    
    if (g_verbose) std::cout << "dt = " << dt << "; iStepInPulse = " << iStepInPulse << "; iStepAfPulse = " << iStepAfPulse
                << std::endl;
    
}

// ---------------------------------------------------------------------------
// RAII session for massflowsolver.mainSolver calls via the Python C API.
// Py_Initialize / module import happen once in the constructor;
// Py_Finalize happens once in the destructor.
// Use MassFlowSolverSession::call() repeatedly inside the time loop.
// ---------------------------------------------------------------------------
class MassFlowSolverSession {
public:
    explicit MassFlowSolverSession(const std::string& scriptDir) {
        Py_Initialize();

        PyObject* sysPath     = PySys_GetObject("path"); // borrowed ref
        PyObject* scriptDirObj = PyUnicode_FromString(scriptDir.c_str());
        PyList_Insert(sysPath, 0, scriptDirObj);
        Py_DECREF(scriptDirObj);

        module_ = PyImport_ImportModule("massflowsolver");
        if (!module_) {
            PyErr_Print();
            Py_Finalize();
            throw std::runtime_error("Python: failed to import massflowsolver");
        }

        func_ = PyObject_GetAttrString(module_, "mainSolver");
        if (!func_ || !PyCallable_Check(func_)) {
            PyErr_Print();
            Py_XDECREF(func_);
            Py_DECREF(module_);
            Py_Finalize();
            throw std::runtime_error("Python: mainSolver not found or not callable");
        }
    }

    ~MassFlowSolverSession() {
        Py_XDECREF(func_);
        Py_XDECREF(module_);
        Py_Finalize();
    }

    MassFlowSolverSession(const MassFlowSolverSession&) = delete;
    MassFlowSolverSession& operator=(const MassFlowSolverSession&) = delete;

    std::map<std::string, double> call(
        const std::map<std::string, double>& ox_Y,
        const std::map<std::string, double>& fuel_Y) const
    {
        auto buildPyDict = [](const std::map<std::string, double>& m) {
            PyObject* d = PyDict_New();
            for (const auto& [k, v] : m) {
                PyObject* val = PyFloat_FromDouble(v);
                PyDict_SetItemString(d, k.c_str(), val);
                Py_DECREF(val);
            }
            return d;
        };

        PyObject* pyOxY   = buildPyDict(ox_Y);
        PyObject* pyFuelY = buildPyDict(fuel_Y);
        PyObject* args    = PyTuple_Pack(2, pyOxY, pyFuelY);
        PyObject* result  = PyObject_CallObject(func_, args);
        Py_DECREF(args);
        Py_DECREF(pyOxY);
        Py_DECREF(pyFuelY);

        if (!result) {
            PyErr_Print();
            throw std::runtime_error("Python: mainSolver execution failed");
        }
        if (!PyDict_Check(result)) {
            Py_DECREF(result);
            throw std::runtime_error("Python: mainSolver did not return a dict");
        }

        std::map<std::string, double> massFlows;
        PyObject *key, *value;
        Py_ssize_t pos = 0;
        while (PyDict_Next(result, &pos, &key, &value))
            massFlows[PyUnicode_AsUTF8(key)] = PyFloat_AsDouble(value);
        Py_DECREF(result);

        return massFlows;
    }

private:
    PyObject* module_ = nullptr;
    PyObject* func_   = nullptr;
};

// Main function
int main(int argc, char *argv[]) {
    printHeader();
    std::cout << "|   Reactor Network module made by: Cristian Casalanguida                     |\n"
         << "|   Contact:      cristian.casalanguida@studenti.polito.it                    |\n"
         << "\\*---------------------------------------------------------------------------*/" << std::endl;
    
    CppBOLOS::currentLogLevel = CppBOLOS::LOG_WARNING;   // Default: only warnings

    // Parse command-line arguments
    std::string casePath = "..";
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-case" && i + 1 < argc) {
            casePath = argv[++i];
        } else if (arg == "-log") {
            g_verbose = true;
            CppBOLOS::currentLogLevel = CppBOLOS::LOG_INFO;
        }
    }

    // All paths derived from casePath
    const std::string controlDictPath = casePath + "/controlDict";
    const std::string parameterPath   = casePath + "/chemPlasProperties";
    const std::string crnDictPath     = casePath + "/crnDict";
    const std::string outputDir       = casePath + "/output";
    const std::string massFlowPath    = outputDir + "/massflow.csv";
    std::filesystem::create_directories(outputDir);
    
    std::ofstream massFlowFile(massFlowPath);
    /* -------------------------------- SET UP PARAMETERS -------------------------------- */

    // Quick sanity checks: ensure required case files exist before attempting to read them
    if (!std::filesystem::exists(controlDictPath)) {
        std::cerr << "Error: controlDict not found at '" << controlDictPath << "'.\n"
                  << "Use the -case <path> option or place a controlDict file in the current directory (" 
                  << std::filesystem::current_path() << ").\n";
        return 1;
    }
    if (!std::filesystem::exists(parameterPath)) {
        std::cerr << "Error: chemPlasProperties not found at '" << parameterPath << "'.\n"
                  << "Use the -case <path> option or place a chemPlasProperties file in the current directory ("
                  << std::filesystem::current_path() << ").\n";
        return 1;
    }
    if (!std::filesystem::exists(crnDictPath)) {
        std::cerr << "Error: crnDict not found at '" << crnDictPath << "'.\n"
                  << "Use the -case <path> option or place a crnDict file in the current directory ("
                  << std::filesystem::current_path() << ").\n";
        return 1;
    }
    // Read time control
    double runTime = readParameter<double>(controlDictPath, "startTime");
    double startTime = runTime;
    double t_end = readParameter<double>(controlDictPath, "endTime");
    const int nPulses = readParameter<int>(controlDictPath, "nPulses");

    double dt1 = readParameter<double>(controlDictPath, "dt1");
    double dt_max = readParameter<double>(controlDictPath, "dt_max");

    // Set T_gas, E/N, and other NRP parameters
    double Tgas = readParameter<double>(parameterPath, "Temperature");
    double P0 = readParameter<double>(parameterPath, "Pressure");
    double EN = readParameter<double>(parameterPath, "E/N(NRP)");
    double EN_DC = readParameter<double>(parameterPath, "E/N(DC)");
    bool stayEN_DC = readParameter<bool>(parameterPath, "stayDC");
    const double N_e_ini = readParameter<double>(parameterPath, "electronDens");
    double Ep = readParameter<double>(parameterPath, "pulseEnergy");
    const double f_NRP = readParameter<double>(parameterPath, "f_NRP");
    const double pulseWidth = readParameter<double>(parameterPath, "pulseWidth");
    double pulseExten = readParameter<double>(controlDictPath, "pulseExten");
    pulseExten = std::min(std::max(pulseExten, pulseWidth), 1.0/f_NRP); // Ensure pulse extension is not shorter than pulse width and not longer than pulse period
    const int NstepInPulse = readParameter<int>(controlDictPath, "NstepInPulse");
    const int NstepAfPulse = readParameter<int>(controlDictPath, "NstepAfPulse");
    double dT_max = readParameter<double>(controlDictPath, "dT_max");
    dT_max = std::max(dT_max, 1.0); // Ensure dT_max is at least 1 K to avoid dt cascade at the beginning
    const double TGAS_TOLERANCE = readParameter<double>(parameterPath, "TGAS_TOLERANCE");
    const double EN_TOLERANCE = readParameter<double>(parameterPath, "EN_TOLERANCE");
    const double dTdt_min = readParameter<double>(controlDictPath, "dTdt_min");
    const double kOutPressure = readParameter<double>(controlDictPath, "kOutPressure");
    // CRN_v4 (optional): interval between mass flow updates and between output rows [s].
    // 0 (default) = every time step, as in CRN_v3.
    const double flowUpdateInterval = readParameter<double>(controlDictPath, "flowUpdateInterval", 0.0);
    const double writeInterval = readParameter<double>(controlDictPath, "writeInterval", 0.0);
    // CRN_v4 (optional): restart CVODE at every time step as in CRN_v3 (1) or only when needed (0, default)
    const bool reinitEveryStep = readParameter<double>(controlDictPath, "reinitEveryStep", 0.0) > 0.5;

    // CRN input parameters: read from crnDict using OpenFOAM-style parser
    auto reactorZones = readReactors(crnDictPath);
    auto massFlowControllers = readMassFlowControllers(crnDictPath);
    massFlowFile << "Time(s)";
    for (const auto& [name, mfc] : massFlowControllers) {
        massFlowFile << ", " << name << "(kg/s)";
    }
    massFlowFile << std::endl;

    // Default initialization
    double dt = dt_max;                     // initial time step
    double tau_NRP = 0.0, nE_3 = 0.0, nE_1 = 0.0, nE_2 = 0.0;  // temporary variables used in the loop
    bool plasmaOn = true;                   // plasma switch flag
    bool dischargeOn = true;                // single nanosecond discharge switch flag
    int iPulse = 0;                         // pulse number index
    double disEnergy = 0.0;                 // initial discharge energy deposited
    const double pulsePeriod = 1.0/f_NRP;
    int iStepInPulse = 0, iStepAfPulse = 0; // step index in/after a pulse discharge
    double T_old = Tgas;
    double dTdt = 0.0, dTdt_max = 1.0;      // initial dT/dt and its maximum
    double time_max_dTdt = 0.0;             // ignition delay time (at maximum dT/dt)
    const double SMALL = 1.0e-16;
    bool fastExpansion = false;             // isentropic gas expansion flag
    t_end = std::max(runTime + pulsePeriod * nPulses, t_end); // overload t_end
    const auto endTimeForce = readParameter<double>(controlDictPath, "endTimeForce", t_end);
    t_end = std::min(t_end, endTimeForce);  // Terminate simulations is endTimeForce is specified

    /* ------------------------------- SET UP BOLTZMANN SOLVER ------------------------------- */
    if (g_verbose) std::cout << "\n========  SETTING BOLTZMANN SOLVER ... ========\n" << std::endl;

    // Species configured by CppBOLOS
    auto bSpecies = readParameter<std::vector<std::string>>(parameterPath, "BoltzmannSpecies");
    auto mixture = readParameter<std::map<std::string, double>>(parameterPath, "Mixture");
    const auto FUELNAME = readParameter<string>(parameterPath, "fuelName");
    const auto OXNAME = readParameter<string>(parameterPath, "oxName");

    std::map<std::string, double> BoltzmannSpecies;

    for (const auto& species : bSpecies) {
        // Check if species is in the mixture, if not assign 0.0
        double fraction = mixture.count(species) > 0 ? mixture[species] : 0.0;
        BoltzmannSpecies[species] = fraction;
        // Be careful of possible difference of species names in LXCat and YAML input,
        // especially for e/E, He/HE, Ar/AR. Always try to use unified names.
    }

    int outputCheckTemp = 0;
    // Set initialized default values for reactors not specified in the input but needed for Boltzmann solver
    for (auto& [name, spec] : reactorZones) {
        if (!spec.mechanismFile.has_value()) {
            spec.mechanismFile = readParameter<string>(parameterPath, "mechFile");
        }
        if (!spec.initialTemperature.has_value()) {
            spec.initialTemperature = Tgas;
        }
        if (!spec.pressure.has_value()) {
            spec.pressure = P0;
        }
        if (!spec.initialElectronDensity.has_value()) {
            spec.initialElectronDensity = N_e_ini;
        }
        if (!spec.composition.has_value()) {
            spec.composition = readParameter<std::map<std::string, double>>(parameterPath, "Mixture");
            spec.boltzmannSpecies = BoltzmannSpecies; // Use global Boltzmann species as default
        } else {
            // If composition is specified, calculate Boltzmann species based on it
            std::map<std::string, double> localBoltzmannSpecies;
            for (const auto& species : bSpecies) {
                // Check if species is in the mixture, if not assign 0.0
                double fraction = spec.composition->count(species) > 0 ? spec.composition->at(species) : 0.0;
                localBoltzmannSpecies[species] = fraction;
                // Be careful of possible difference of species names in LXCat and YAML input,
                // especially for e/E, He/HE, Ar/AR. Always try to use unified names.
            }
            spec.boltzmannSpecies = std::move(localBoltzmannSpecies);
        }
        if (!spec.outputFile.has_value()) {
            spec.outputFile = outputDir + "/" + name + ".csv";
        }
        if (spec.type != "reservoir" && spec.volume == 0.0) {
            std::cerr << "Error: volume must be specified and greater than 0 for non-reservoir reactor '" << name << "'.\n";
            outputCheckTemp = 1;
        }
    }
    if (outputCheckTemp) {
        std::cerr << "Please fix the above issues in the input files and try again.\n";
        return 1;
    }

    /* ----------------------------- SET UP GAS PHASE AND REACTOR  -------------------------------- */
    if (g_verbose) std::cout << "\n========  SETTING GAS PHASE & REACTORS, CREATE  & INIT ODE INTEGRATOR ========\n" << std::endl;
    /* ------------------------------ CREATE & INIT ODE INTEGRATOR --------------------------------- */
    //  - the default settings for CVodesIntegrator are used:
    //     solution method: BDF_Method
    //     problem type: DENSE + NOJAC
    //     relative tolerance: 1.0e-9
    //     absolute tolerance: 1.0e-15
    //     max step size: +inf

    // Read plasma heating model
    std::string plasmaHeatModel = " ";
    try {
        plasmaHeatModel = readParameter<string>(parameterPath, "plasmaHeatModel");
    } catch (const std::runtime_error& e) {
        std::cerr << "Warning: " << e.what() << std::endl;
    }

    std::string CP_or_CV = readParameter<string>(parameterPath, "CP_or_CV");

    std::vector<ReactorState> reactors;
    std::vector<ReservoirState> reservoirs;
    // Reserve only for reactors that will be created (exclude "reservoir" types)
    size_t n_non_reservoir = 0;
    for (const auto& [name, cfg] : reactorZones) {
        if (cfg.type != "reservoir") ++n_non_reservoir;
    }
    if (n_non_reservoir == 0) {
    std::cerr << "Error: no active reactors found (all entries are of type 'reservoir')." << std::endl;
    return 1;
    }
    reactors.reserve(n_non_reservoir);

    std::vector<double> nE_1_vec(n_non_reservoir, nE_1);
    std::vector<double> nE_2_vec(n_non_reservoir, nE_2);
    std::vector<double> nE_3_vec(n_non_reservoir, nE_3);
    std::vector<double> EN_vec(n_non_reservoir, EN);
    std::vector<uint8_t> dischargeOn_vec(n_non_reservoir, dischargeOn);
    std::vector<uint8_t> plasmaOn_vec(n_non_reservoir, plasmaOn);
    std::vector<uint8_t> fastExpansion_vec(n_non_reservoir, fastExpansion);
    std::vector<double> disEnergy_vec(n_non_reservoir, disEnergy);
    std::vector<int> iPulse_vec(n_non_reservoir, iPulse);
    std::vector<int> iStepInPulse_vec(n_non_reservoir, iStepInPulse);
    std::vector<int> iStepAfPulse_vec(n_non_reservoir, iStepAfPulse);
    std::vector<double> tau_NRP_vec(n_non_reservoir, tau_NRP);
    std::vector<double> T_old_vec(n_non_reservoir, T_old);
    std::vector<double> dT_vec(n_non_reservoir, 0.0);
    std::vector<double> dTdt_vec(n_non_reservoir, dTdt);
    std::vector<double> dTdt_max_vec(n_non_reservoir, dTdt_max);
    std::vector<double> time_max_dTdt_vec(n_non_reservoir, time_max_dTdt);

    // Initialize Boltzmann solver grid and cross-sections once for all reactors
    EEDF::Settings eedfSettings;
    {
        eedfSettings.csDataFile = readParameter<string>(parameterPath, "csDataFile");
        auto gridSize = readParameter<std::map<std::string, double>>(parameterPath, "gridSize");
        eedfSettings.gridType = readParameter<string>(parameterPath, "gridType");
        eedfSettings.gridStart = gridSize["start"];
        eedfSettings.gridEnd = gridSize["end"];
        eedfSettings.gridPoints = (int)gridSize["points"];
        // EEDF convergence of converge(): default as CRN_v4; tighter values make the EEDF
        // independent of the initial guess (history of the reactor)
        eedfSettings.rtol = readParameter<double>(parameterPath, "eedfRtol", 1e-5);
        eedfSettings.maxItr = (int)readParameter<double>(parameterPath, "eedfMaxIterations", 200.0);
        EEDF::setupSolver(BoltzmannRate::bsolver, eedfSettings);
        LOG_INFO("\nA total of " + std::to_string(BoltzmannRate::bsolver.number_of_targets()) +
                 " targets have been loaded:\n" + BoltzmannRate::bsolver.targetNames());
    }

    int idxState_temp = 0;
    int idxReservoir_temp = 0;
    for (auto& [name, config] : reactorZones) {
        if (config.type != "reservoir") {
            reactors.push_back(createReactorState(name, config, plasmaHeatModel, parameterPath, controlDictPath, EN, runTime));
            config.indexState = idxState_temp++;
        } else {
            reservoirs.push_back(createReservoirState(name, config, runTime));
            config.indexState = idxReservoir_temp++;
        }
    }

    // CRN_v4: per-reactor Boltzmann state. All Solutions (hence all BoltzmannRate
    // objects) exist now, so NumProcess is final and the caches can be sized.
    // The cached path of BoltzmannRate::evalFromStruct ignores branch_ratio:
    // refuse mechanisms that use it.
    for (const auto& [name, config] : reactorZones) {
        std::ifstream mech(config.mechanismFile.value());
        std::stringstream buffer;
        buffer << mech.rdbuf();
        if (buffer.str().find("branch_ratio") != std::string::npos) {
            std::cerr << "Error: mechanism '" << config.mechanismFile.value() << "' uses branch_ratio, "
                      << "not supported by the cached Boltzmann rates of CRN_v4/v5. Use CRN_v3." << std::endl;
            return 1;
        }
    }

    // CRN_v5: EEDF table. Optional file of the case (default <case>/eedfTable.yaml, "none" = in memory
    // only); it is read at start and every EEDF computed during the run is appended to it.
    std::vector<std::string> boltzmannProcesses; // distinct process strings of the mechanisms
    for (const auto& process : BoltzmannRate::ProcessesList) {
        if (std::find(boltzmannProcesses.begin(), boltzmannProcesses.end(), process) == boltzmannProcesses.end()) {
            boltzmannProcesses.push_back(process);
        }
    }
    std::string eedfTablePath = readParameter<string>(parameterPath, "eedfTable", std::string("<case>/eedfTable.yaml"));
    if (eedfTablePath.rfind("<case>", 0) == 0) eedfTablePath.replace(0, 6, casePath);
    if (eedfTablePath == "none") eedfTablePath.clear();
    EEDF::Tolerance eedfTolerance;
    eedfTolerance.EN = EN_TOLERANCE;
    eedfTolerance.T = TGAS_TOLERANCE;
    eedfTolerance.maxGap = readParameter<double>(parameterPath, "eedfTableMaxGap", eedfTolerance.maxGap);
    eedfTolerance.Xrel = readParameter<double>(parameterPath, "eedfTableXrel", eedfTolerance.Xrel);
    eedfTolerance.Xabs = readParameter<double>(parameterPath, "eedfTableXabs", eedfTolerance.Xabs);
    EEDF::Table eedfTable(eedfTablePath, eedfSettings, eedfTolerance);
    std::cout << "EEDF table: " << (eedfTablePath.empty() ? "in memory only" : eedfTablePath)
              << " (" << eedfTable.nLoaded() << " points loaded)" << std::endl;

    for (auto& state : reactors) {
        const ReactorZoneDef& config = reactorZones.at(state.name);
        computeEEDF(state, EN, config.initialTemperature.value(), EN, 1e14, // same as the initial converge() of CRN_v3
                    eedfTable, eedfSettings, boltzmannProcesses);
    }
    BoltzmannRate::useCachedSolutions = true;

    // Map end species -> source species once for every mass flow controller
    for (auto& [name, mfc] : massFlowControllers) {
        const ReactorZoneDef& configStart = reactorZones.at(mfc.start);
        const ReactorZoneDef& configEnd = reactorZones.at(mfc.end);
        shared_ptr<ThermoPhase> srcGas = (configStart.type != "reservoir")
            ? reactors[configStart.indexState].gas : reservoirs[configStart.indexState].gas;
        shared_ptr<ThermoPhase> endGas = (configEnd.type != "reservoir")
            ? reactors[configEnd.indexState].gas : reservoirs[configEnd.indexState].gas;
        mfc.srcIdx.resize(endGas->nSpecies());
        for (size_t k = 0; k < endGas->nSpecies(); k++) {
            mfc.srcIdx[k] = srcGas->speciesIndex(endGas->speciesName(k));
        }
        mfc.sumYdt.assign(endGas->nSpecies(), 0.0);
    }

    MassFlowSolverSession pySession(casePath);

    /* ------------------------------------- RUN THE SIMULATION ---------------------------------- */
    std::cout << "\n======== RUNNING SIMULATION ========\n" << std::endl;

    clock_t t0 = clock(); // save start time
    double dTdtcheck = std::abs(dTdt_min) * 10; // initial value to enter the loop
    const double TIME_TOL = 1e-12;    // tolerance on time comparisons for flow updates and output
    double nextFlowUpdate = runTime;  // CRN_v4: first flow update at the first step
    double nextWrite = runTime + writeInterval;
    bool lastStepWritten = true;      // initial state already written in createReactorState
    size_t step = 0;
    int lastPrintedPulse = -1;
    clock_t lastPrint = clock();
    // Main time loop
    // const auto dTbelowTeq = readParameter<double>(controlDictPath, "dTbelowTeq");
    // while (runTime < t_end && gas->temperature() < (T_eq - dTbelowTeq)) {
    while (runTime < t_end && dTdtcheck > dTdt_min) { // dTdt_min < 0 disables the dT/dt stop

        for (size_t idx = 0; idx < reactors.size(); idx++) {
            T_old_vec[idx] = reactors[idx].gas->temperature();
        }

        // HERE WE CALCULATE THE MASS FLOW CONTROLLERS VALUES.
        // CRN_v4: done every flowUpdateInterval (every step if 0). Between two updates the
        // source state passed downstream is the time average of the source reactor over
        // the previous interval, so a downstream reactor sees the pulse-averaged outflow.
        const bool flowsUpdated = (flowUpdateInterval <= 0.0) || (runTime + TIME_TOL >= nextFlowUpdate);
        if (flowsUpdated) {
            std::map<std::string, double> yFuel, yOx;
            for (const auto& state : reactors) {
                yFuel[state.name] = state.gas->elementalMassFraction(state.gas->elementIndex(FUELNAME));
                yOx[state.name] = state.gas->elementalMassFraction(state.gas->elementIndex(OXNAME));
            }
            for (const auto& state : reservoirs) {
                yFuel[state.name] = state.gas->elementalMassFraction(state.gas->elementIndex(FUELNAME));
                yOx[state.name] = state.gas->elementalMassFraction(state.gas->elementIndex(OXNAME));
            }

            std::map<std::string, double> massFlowRates = pySession.call(yOx, yFuel);

            massFlowFile << std::setprecision(15) << runTime << std::setprecision(6);
            for (auto& [name, mfc] : massFlowControllers) {
                auto it = massFlowRates.find(name);
                if (it == massFlowRates.end())
                    throw std::runtime_error("massFlowSolver: no key '" + name + "' in Python output; check controller names in crnDict match massflowsolver.py");
                mfc.value = it->second;
                massFlowFile << ", " << mfc.value;

                const ReactorZoneDef& configStart = reactorZones.at(mfc.start);
                shared_ptr<ThermoPhase> srcGas = (configStart.type != "reservoir")
                    ? reactors[configStart.indexState].gas
                    : reservoirs[configStart.indexState].gas;

                const size_t nSpec = mfc.srcIdx.size();
                mfc.comp.assign(nSpec, 0.0);
                if (mfc.sumDt > 0.0) {
                    for (size_t k = 0; k < nSpec; k++)
                        mfc.comp[k] = mfc.sumYdt[k] / mfc.sumDt;
                    mfc.T = mfc.sumTdt / mfc.sumDt;
                    mfc.h = mfc.sumHdt / mfc.sumDt;
                } else {
                    for (size_t k = 0; k < nSpec; k++)
                        mfc.comp[k] = (mfc.srcIdx[k] != Cantera::npos) ? srcGas->massFraction(mfc.srcIdx[k]) : 0.0;
                    mfc.T = srcGas->temperature();
                    mfc.h = srcGas->enthalpy_mass();
                }
                std::fill(mfc.sumYdt.begin(), mfc.sumYdt.end(), 0.0);
                mfc.sumHdt = mfc.sumTdt = mfc.sumDt = 0.0;
            }
            massFlowFile << "\n";

            // For each reservoir that receives inflows, compute the mixed state from
            // all contributing mass flows (ignoring the existing reservoir content).
            {
                struct ResAccum {
                    std::vector<double> sumYm; // Σ(Y_src[k] * mdot) per species k
                    double sumM  = 0.0;        // Σ(mdot) [kg/s]
                    double sumHM = 0.0;        // Σ(h_src * mdot) [J/s]
                };
                std::map<std::string, ResAccum> resAccum;

                for (const auto& [name, mfc] : massFlowControllers) {
                    ReactorZoneDef& configEnd = reactorZones.at(mfc.end);
                    if (configEnd.type != "reservoir" || mfc.value <= 0.0) continue;

                    ReservoirState& resEnd = reservoirs[configEnd.indexState];
                    const size_t nSpec = resEnd.gas->nSpecies();
                    ResAccum& accum = resAccum[mfc.end];
                    if (accum.sumYm.empty()) accum.sumYm.assign(nSpec, 0.0);

                    for (size_t k = 0; k < nSpec; k++)
                        accum.sumYm[k] += mfc.comp[k] * mfc.value;
                    accum.sumM  += mfc.value;
                    accum.sumHM += mfc.h * mfc.value;
                }

                for (auto& [resName, accum] : resAccum) {
                    if (accum.sumM <= 0.0) continue;
                    ReactorZoneDef& configRes = reactorZones.at(resName);
                    ReservoirState& resState  = reservoirs[configRes.indexState];
                    const size_t nSpec = resState.gas->nSpecies();

                    double h_mix = accum.sumHM / accum.sumM;
                    std::vector<double> Y_mix(nSpec);
                    for (size_t k = 0; k < nSpec; k++)
                        Y_mix[k] = accum.sumYm[k] / accum.sumM;

                    resState.gas->setMassFractions(Y_mix.data());
                    resState.gas->setState_HP(h_mix, configRes.pressure.value());
                }
            }

            // Flow terms of every reactor; CVODE is restarted only when they change
            for (auto& state : reactors) {
                const ReactorZoneDef& config = reactorZones.at(state.name);
                std::vector<ChemPlasReactor::FlowTerm> inflows;
                state.mdot_in = 0.0;
                state.mdot_out = 0.0;
                for (auto& [name, mfc] : massFlowControllers) {
                    if (mfc.end == state.name) {
                        inflows.push_back({mfc.value / config.volume, mfc.comp, mfc.T, mfc.h});
                        state.mdot_in += mfc.value;
                    }
                    if (mfc.start == state.name) {
                        state.mdot_out += mfc.value;
                    }
                }
                state.odes->setFlows(inflows, state.mdot_out / config.volume);
                state.needsReinit = true;
            }

            if (flowUpdateInterval > 0.0) {
                while (nextFlowUpdate <= runTime + TIME_TOL) nextFlowUpdate += flowUpdateInterval;
            }
        }

        double dt_old = dt;
        dt = dt_max; // reset dt to max at the beginning of each loop, will be updated in findTimeStep
        for (size_t idx = 0; idx < reactors.size(); idx++) {
            ReactorState& state = reactors[idx];
            ReactorZoneDef& config = reactorZones.at(state.name);
            double dt_element = dt_old;
            double& EN_temp = EN_vec[idx];
            double& nE_1 = nE_1_vec[idx];
            double& nE_2 = nE_2_vec[idx];
            double& nE_3 = nE_3_vec[idx];
            uint8_t& dischargeOn = dischargeOn_vec[idx];
            double& disEnergy = disEnergy_vec[idx];
            int& iPulse = iPulse_vec[idx];
            int& iStepInPulse = iStepInPulse_vec[idx];
            int& iStepAfPulse = iStepAfPulse_vec[idx];
            double& tau_NRP = tau_NRP_vec[idx];
            uint8_t& fastExpansion = fastExpansion_vec[idx];
            double& dT = dT_vec[idx];
            double& dTdt = dTdt_vec[idx];

            findTimeStep(state, config, runTime, startTime, SMALL, pulsePeriod, pulseWidth, pulseExten, dischargeOn, fastExpansion, stayEN_DC,
                        CP_or_CV, iPulse, nPulses, dt_element, dt1, dt_max, dT, dT_max, dTdt,
                        nE_1, nE_2, nE_3, tau_NRP, disEnergy,
                        iStepInPulse, NstepInPulse,
                        iStepAfPulse, NstepAfPulse,
                        EN, EN_DC, EN_temp, EN_TOLERANCE, TGAS_TOLERANCE);

            dt = std::min(dt, dt_element); // Take the minimum dt among reactors to ensure stability
        }

        // CRN_v4: progress line once per pulse or every 2 s of CPU time (not every step)
        step++;
        {
            const int pulseNow = *std::max_element(iPulse_vec.begin(), iPulse_vec.end());
            const clock_t now = clock();
            if (pulseNow != lastPrintedPulse || static_cast<double>(now - lastPrint) / CLOCKS_PER_SEC > 2.0) {
                const auto oldFlags = std::cout.flags();
                const auto oldPrecision = std::cout.precision();
                std::cout << "\rrunTime [s]: " << std::scientific << std::setprecision(9) << runTime
                          << " | iPulse: " << pulseNow << " | steps: " << step << "   " << std::flush;
                std::cout.flags(oldFlags);
                std::cout.precision(oldPrecision);
                lastPrintedPulse = pulseNow;
                lastPrint = now;
            }
        }

        /*----------------------  Advance time step --------------------*/
        const double runTimeOld = runTime;
        runTime += dt;

        // CRN_v4: write output every writeInterval (every step if 0) and at the end time
        const bool doWrite = (writeInterval <= 0.0) || (runTime + TIME_TOL >= nextWrite) || (runTime + TIME_TOL >= t_end);
        if (doWrite && writeInterval > 0.0) {
            while (nextWrite <= runTime + TIME_TOL) nextWrite += writeInterval;
        }
        lastStepWritten = doWrite;

        for (size_t idx = 0; idx < reactors.size(); idx++) {
            ReactorState& state = reactors[idx];
            ReactorZoneDef& config = reactorZones.at(state.name);
            double& disEnergy = disEnergy_vec[idx];
            uint8_t& dischargeOn = dischargeOn_vec[idx];
            uint8_t& plasmaOn = plasmaOn_vec[idx];
            int& iPulse = iPulse_vec[idx];
            int& iStepInPulse = iStepInPulse_vec[idx];
            int& iStepAfPulse = iStepAfPulse_vec[idx];
            double& T_old = T_old_vec[idx];
            double& dT = dT_vec[idx];
            double& dTdt = dTdt_vec[idx];
            double& dTdt_max = dTdt_max_vec[idx];
            double& time_max_dTdt = time_max_dTdt_vec[idx];
            double& nE_1 = nE_1_vec[idx];
            double& nE_2 = nE_2_vec[idx];
            double& EN_temp = EN_vec[idx];

            // CRN_v4: install this reactor's Boltzmann rate constants; the EEDF is
            // recomputed only when E/N, T or the CppBOLOS mixture change beyond tolerance
            activateReactor(state);
            bool updateEEDF = state.odes->updateBoltzmannMixture();
            double EN_new = state.EN_b;
            double T_new = state.kT_b;
            if (std::abs(EN_temp - state.EN_b) > EN_TOLERANCE) {
                EN_new = EN_temp;
                updateEEDF = true;
            }
            if (std::abs(state.gas->temperature() - state.kT_b) > TGAS_TOLERANCE) {
                T_new = state.gas->temperature();
                updateEEDF = true;
            }
            if (updateEEDF) {
                if (g_verbose) std::cout << "Updating EEDF ..." << "\n";
                computeEEDF(state, EN_new, T_new, EN_temp, 1E20/(EN_temp*EN_temp), eedfTable, eedfSettings, boltzmannProcesses);
                activateReactor(state);
                state.needsReinit = true;
            }

            if (iPulse == 0 || reinitEveryStep) {
                state.needsReinit = true;
            }

            if (g_verbose) std::cout << "EN = " << state.EN_b
                        << "; Tgas = " << state.gas->temperature()
                        << "; Te = " << state.TeSolver
                        << "; n_e = " << getNumberDens(state.gas, state.odes->electronIndex)
                        << "\n";

            // CRN_v4: restart CVODE only when something outside the ODE changed
            // (new EEDF, new flows, energy reset at a new pulse, heat-loss jump)
            // Discharge on in this step: same condition as findTimeStep
            const bool discharging = config.discharge && iPulse < nPulses && dischargeOn
                                     && (runTimeOld + SMALL < startTime + iPulse * pulsePeriod + pulseWidth);
            // integrate() needs CVODE at runTimeOld; after step() it can be ahead
            if (!discharging && state.tInternal > runTimeOld) state.needsReinit = true;
            if (discharging && state.tInternal < runTimeOld) state.needsReinit = true;

            if (state.needsReinit) {
                state.integrator->reinitialize(runTimeOld, *state.odes);
                state.needsReinit = false;
                state.tInternal = runTimeOld;
            }

            if (discharging) {
                // Discharge ended exactly at Ep: internal steps only for this reactor and only now
                const double tStar = integrateDischarge(state, runTimeOld, runTime, 1e3 * Ep); // Ep [mJ/cm^3] -> J/m^3
                if (tStar >= 0.0) {
                    // Discharge off at t*: E/N back to the base value, new EEDF, CVODE restarted at t*,
                    // then a step to runTime with the discharge off
                    tau_NRP_vec[idx] = tStar - (startTime + iPulse * pulsePeriod);
                    nE_2 = getNumberDens(state.gas, state.odes->electronIndex);
                    state.EpDischarge = 1e-3 * state.odes->depositedPlasmaEnergy();
                    dischargeOn = false;
                    state.pulseEndFound = true;
                    if (g_verbose) std::cout << "NSD terminated at Ep, tau_NRP = " << tau_NRP_vec[idx] << "\n";
                    EN_temp = std::max(EN_DC, 0.1);
                    computeEEDF(state, EN_temp, state.gas->temperature(), EN_temp, 1E20/(EN_temp*EN_temp),
                            eedfTable, eedfSettings, boltzmannProcesses);
                    activateReactor(state);
                    state.integrator->reinitialize(tStar, *state.odes);
                    // CVODE cannot integrate over a few ulps of t: below 100 ulps the state at runTime is
                    // the one at t* (that interval is integrated with the next step anyway)
                    if (runTime - tStar > 100 * std::numeric_limits<double>::epsilon() * std::abs(runTime)) {
                        state.integrator->integrate(runTime);
                        state.odes->updateState(state.integrator->solution());
                    }
                    state.tInternal = -1.0;
                }
            } else {
                state.integrator->integrate(runTime);
                state.odes->updateState(state.integrator->solution());
                state.tInternal = -1.0;
            }

            // For CV reactors: update m_density to account for net mass flux.
            // For CP reactors: density follows the EOS automatically — no update needed.
            if (!state.odes->constPressure) {
                double net_mdot_V = (state.mdot_in - state.mdot_out) / config.volume;
                double rho_new = state.gas->density() + net_mdot_V * dt;
                state.gas->setState_TD(state.gas->temperature(), std::max(rho_new, 1e-6));
                state.odes->setConstPD(state.gas->pressure(), std::max(rho_new, 1e-6));
            }

            dT = std::abs(state.gas->temperature() - T_old);

            // Check plasma power and energy deposited
            const double ne = getNumberDens(state.gas, state.odes->electronIndex);
            state.neMaxPulse = std::max(state.neMaxPulse, ne);
            double power_temp = state.elecPower
                                * Avogadro * state.gas->molarDensity()* ElectronCharge * ne;
            // unit: eV m^3/s * #/kmol * kmol/m^3 * coulomb(J/eV) * #/cm^3 = J/s/cm^3
            disEnergy = 1e-3 * state.odes->depositedPlasmaEnergy(); // mJ/cm^3

            if (doWrite) writeReactorRow(state, runTime, power_temp, disEnergy);
            if (g_verbose) std::cout << "Power_input [W/cm^3]: " << power_temp << " | Energy deposited [mJ/cm^3]: " << disEnergy << "\n";

            // Check if terminate nanosecond discharge, cut off by designated power deposition Ep
            dischargeOn = (disEnergy < Ep);

            // Update dT/dt
            dTdt =  dT / dt;
            // Update dTdt_max and its time
            if (dTdt > dTdt_max && iStepAfPulse > 1) {
                dTdt_max = dTdt;
                time_max_dTdt = runTime;
            }

            // Trigger plasmaOn only ONCE for a new pulse
            if (runTime + SMALL >= ((iPulse + 1) * pulsePeriod + startTime)) { // Entering next pulse period
                if (state.pulseFile.is_open()) {
                    state.pulseFile << iPulse << ", " << std::setprecision(15) << runTime << std::setprecision(6) << ", " << state.gas->temperature() << ", "
                                    << state.gas->pressure() << ", " << ne << ", " << state.neMaxPulse << ", "
                                    << disEnergy << ", " << tau_NRP_vec[idx] << ", " << state.EpDischarge << "\n";
                }
                state.neMaxPulse = ne;
                iPulse += 1;
                if (iPulse == nPulses) {
                    std::cout << "\n >>>>>>>> NRP Discharges Finished for reactor " << state.name << " <<<<<<<< \n" << std::endl;
                    state.integrator->reinitialize(runTime, *state.odes);
                    plasmaOn = false;
                }
                if (plasmaOn) {
                    if (g_verbose) std::cout << "\n**** Start Pulse Loop No. " << iPulse << " at runTime = " << runTime << " for reactor " << state.name << " ****\n";
                    if (state.odes->heatLoss) {
                        state.odes->setConstPD(state.gas->pressure(), state.gas->density());
                        state.odes->constPressure = false; // switch to constant volume
                    }

                    disEnergy = 0.0, dischargeOn = true;
                    state.odes->resetDepositedPlasmaEnergy();
                    state.needsReinit = true; // Ep is part of the CVODE state vector
                    iStepInPulse = 0, iStepAfPulse = 0;
                    state.pulseEndFound = false;
                    nE_1 = ne; // n_e at the beginning of a new pulse
                }
            }
        }

        // CRN_v4: time integrals of the reactor outflows, averaged at the next flow update
        if (flowUpdateInterval > 0.0) {
            for (auto& [name, mfc] : massFlowControllers) {
                const ReactorZoneDef& configStart = reactorZones.at(mfc.start);
                if (configStart.type == "reservoir") continue;
                const auto& srcGas = reactors[configStart.indexState].gas;
                for (size_t k = 0; k < mfc.srcIdx.size(); k++) {
                    if (mfc.srcIdx[k] != Cantera::npos)
                        mfc.sumYdt[k] += srcGas->massFraction(mfc.srcIdx[k]) * dt;
                }
                mfc.sumHdt += srcGas->enthalpy_mass() * dt;
                mfc.sumTdt += srcGas->temperature() * dt;
                mfc.sumDt  += dt;
            }
        }

        // Write updated reservoir states at current time step
        if (doWrite) {
            for (auto& res : reservoirs) {
                res.outputFile << std::setprecision(15) << runTime << std::setprecision(6) << ", " << res.gas->temperature() << ", "
                               << res.gas->pressure() << ", "
                               << 1e-6 * Avogadro * res.gas->molarDensity() << ", "
                               << res.gas->meanMolecularWeight();
                for (const int idx : res.indexList)
                    res.outputFile << ", " << getNumberDens(res.gas, static_cast<size_t>(idx));
                res.outputFile << "\n";
            }
        }

        dTdtcheck = *std::max_element(dTdt_vec.begin(), dTdt_vec.end());

        if (g_verbose) std::cout << "-------------" << "\n";
    }
    // CRN_v4: make sure the final state is in the output files
    if (!lastStepWritten) {
        for (auto& state : reactors) {
            const double ne = getNumberDens(state.gas, state.odes->electronIndex);
            const double power = state.elecPower * Avogadro * state.gas->molarDensity() * ElectronCharge * ne;
            writeReactorRow(state, runTime, power, 1e-3 * state.odes->depositedPlasmaEnergy());
        }
        for (auto& res : reservoirs) {
            res.outputFile << std::setprecision(15) << runTime << std::setprecision(6) << ", " << res.gas->temperature() << ", "
                           << res.gas->pressure() << ", "
                           << 1e-6 * Avogadro * res.gas->molarDensity() << ", "
                           << res.gas->meanMolecularWeight();
            for (const int idx : res.indexList)
                res.outputFile << ", " << getNumberDens(res.gas, static_cast<size_t>(idx));
            res.outputFile << "\n";
        }
    }
    // End of the main time loop
    std::cout << "\n";

    for (size_t idx = 0; idx < reactors.size(); idx++) {
            ReactorState& state = reactors[idx];
            double& time_max_dTdt = time_max_dTdt_vec[idx];
             std::cout << "Reactor " << state.name << ": ";
             if (time_max_dTdt > 0) {
                 std::cout << "Ignition Delay Time [s]: " << time_max_dTdt << std::endl;
             } else {
                 std::cout << "No ignition detected within the simulation time." << std::endl;
             }
    }
    
    //std::cout << "\nIgnition Delay Time [s]: " << time_max_dTdt << std::endl;

    clock_t t1 = clock(); // save end time
    double elapsed_time = static_cast<double>(t1 - t0) / CLOCKS_PER_SEC;

    if (runTime >= t_end) {
        std::cout << "\nSimulation terminated at end time " << t_end << " s." << std::endl;
    } else if (dTdtcheck <= dTdt_min) {
        std::cout << "\nSimulation terminated at runTime = " << runTime << " s due to max(dT/dt) <= " << dTdt_min << " K/s." << std::endl;
    } else {
        std::cout << "\nSimulation terminated at runTime = " << runTime << " s due to unknown reason." << std::endl;
    }
    std::cout << "EEDF: " << g_eedfFromTable << " from table, " << g_eedfComputed << " computed"
              << (eedfTablePath.empty() ? "" : " and added to " + eedfTablePath) << std::endl;
    std::cout << "\n======== ChemPlasKin finished. Execution time: " << elapsed_time << " (sec) ========\n" <<std::endl;

    return 0;
}