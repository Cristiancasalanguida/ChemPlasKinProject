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
|   Last edit:    2026-05-20                                                  |
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
#include  "plasmaReactor.h" // for original one "../../src/plasmaReactor.h"
#include "../../src/utilities.h"
#include "../../include/kinetics/BoltzmannGuard.h"


// Struct used for reading reactor entries
struct ReactorZoneDef {
    std::string type;
    bool discharge = false;
    double volume = 0.0;
    int indexState = -1;
    std::optional<std::string> mechanismFile; // Filled from global default if not provided
    std::optional<std::string> outputFile; // Filled from global default if not provided
    std::optional<double> initialTemperature; // Filled from global default if not provided
    std::optional<double> pressure; // Filled from global default if not provided
    std::optional<double> initialElectronDensity; // Filled from global default if not provided
    std::optional<std::map<std::string, double>> composition; // Filled from global default if not provided
    std::optional<std::map<std::string, double>> boltzmannSpecies; // Calculated from composition when present
};

// Struct for mass flow controller entries
struct MassFlowControllerDef {
    std::string start;
    std::string end;
    double value = 0.0;
    double T = 0.0;
    std::vector<double> comp;
    double h = 0.0; // specific enthalpy of source [J/kg]
};

struct ReactorState {
    std::string name;
    std::shared_ptr<std::map<std::string, double>> boltzmannSpecies;
    std::shared_ptr<Solution> sol;
    std::shared_ptr<ThermoPhase> gas;
    BoltzmannSnapshot boltzmannState;
    std::unique_ptr<ChemPlasReactor> odes;
    std::unique_ptr<Integrator> integrator;
    std::ofstream outputFile;
    std::vector<int> indexList;
    std::vector<double> qf;
    std::vector<double> qr;
    std::vector<double> q;
};

struct ReservoirState {
    std::string name;
    std::shared_ptr<std::map<std::string, double>> boltzmannSpecies;
    std::shared_ptr<Solution> sol;
    std::shared_ptr<ThermoPhase> gas;
    std::ofstream outputFile;
    std::vector<int> indexList;
};

std::string makeUniqueReactorName(const std::string& baseName, const std::map<std::string, ReactorZoneDef>& reactors) {
    if (reactors.find(baseName) == reactors.end()) {
        return baseName;
    }

    int suffix = 1;
    while (reactors.find(baseName + std::to_string(suffix)) != reactors.end()) {
        ++suffix;
    }

    return baseName + std::to_string(suffix);
}

// Parse reactors block from OpenFOAM-style dict
std::map<std::string, ReactorZoneDef> readReactors(const std::string& fileName) {
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

    std::map<std::string, ReactorZoneDef> reactors;
    std::string line;
    bool inReactorsBlock = false;

    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key)) continue;

        // Find reactors block
        if (!inReactorsBlock && key == "reactors") {
            inReactorsBlock = true;
            continue;
        }

        if (!inReactorsBlock) continue;

        // End of reactors block
        if (line.find('}') != std::string::npos && inReactorsBlock) {
            // Check if this closes the reactors block
            int braceCount = 0;
            for (char c : line) {
                if (c == '{') braceCount++;
                if (c == '}') braceCount--;
            }
            if (braceCount < 0) break; // Closing brace for reactors block
        }

        // Read reactor name and open brace
        if (key != "{" && key != "}") {
            std::string reactorName = key;
            ReactorZoneDef spec;

            // Find opening brace
            if (line.find('{') == std::string::npos) {
                while (std::getline(inFile, line) && line.find('{') == std::string::npos) {}
            }

            // Read reactor definition until closing brace
            while (std::getline(inFile, line)) {
                if (line.find('}') != std::string::npos) break;

                std::istringstream ls(line);
                std::string field;
                if (!(ls >> field)) continue;
                field = normalizeToken(field);

                if (field == "type") {
                    std::string val;
                    if (ls >> val) {
                        spec.type = normalizeToken(val);
                    }
                } else if (field == "discharge") {
                    std::string val;
                    if (ls >> val) {
                        val = normalizeToken(val);
                        spec.discharge = (val == "yes" || val == "true" || val == "1");
                    }
                } else if (field == "mechanismFile") {
                    std::string val;
                    if (ls >> val) spec.mechanismFile = normalizeToken(val);
                } else if (field == "outputFile") {
                    std::string val;
                    if (ls >> val) spec.outputFile = normalizeToken(val);
                } else if (field == "initialTemperature") {
                    double val;
                    if (ls >> val) spec.initialTemperature = val;
                } else if (field == "pressure") {
                    double val;
                    if (ls >> val) spec.pressure = val;
                } else if (field == "initialElectronDensity") {
                    double val;
                    if (ls >> val) spec.initialElectronDensity = val;
                } else if (field == "volume") {
                    double val;
                    if (ls >> val) spec.volume = val;
                } else if (field == "composition") {
                    // Find opening brace for composition
                    if (line.find('{') == std::string::npos) {
                        while (std::getline(inFile, line) && line.find('{') == std::string::npos) {}
                    }
                    // Read composition entries
                    std::map<std::string, double> comp;
                    while (std::getline(inFile, line)) {
                        if (line.find('}') != std::string::npos) break;
                        std::istringstream cs(line);
                        std::string species;
                        double frac;
                        if (cs >> species >> frac) {
                            comp[species] = frac;
                        }
                    }
                    spec.composition = std::move(comp);
                }
            }
            reactorName = makeUniqueReactorName(reactorName, reactors);
            reactors[reactorName] = std::move(spec);
        }
    }

    return reactors;
}

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
    BoltzmannRate::F0 = BoltzmannRate::bsolver.maxwell(2.0); // initial guess from Maxwell EEDF
    BoltzmannRate::F0 = BoltzmannRate::bsolver.converge(BoltzmannRate::F0, 200, 1e-5);
    //

    state.boltzmannState = BoltzmannSnapshot::capture();
    state.boltzmannState.density = boltzmannSpecies;

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
        std::cout << "Start with constant-volume PAC reactor since heat loss model is on."<< std::endl;
        state.odes->constPressure = false;
    } else
    {
        if (readParameter<string>(parameterPath, "CP_or_CV") == "CP"){
            std::cout << "Start with constant-pressure PAC reactor"<< std::endl;
            state.odes->constPressure = true;
        } else {
            std::cout << "Start with constant-volume PAC reactor"<< std::endl;
            state.odes->constPressure = false;
        }
    }
    state.odes->inertSpIndex = state.odes->findSpeciesIndex(readParameter<string>(parameterPath, "inertSpecies"));

    // Create and initialize the ODE integrator

    state.integrator = std::unique_ptr<Integrator>(newIntegrator("CVODE"));
    state.integrator->initialize(runTime, *state.odes);

    auto Tolerances = readParameter<std::map<std::string, double>>(controlDictPath, "odes");
    state.integrator->setTolerances(Tolerances["reltol"], Tolerances["abstol"]);

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
    std::cout << "Writing information of " << name << ":\n" << "Time, T_gas(K), p(Pa), N_gas(#/cm^3), MW(kg/kmol),";
    for (const auto& sp : speciesNames) {
        std::cout << ", " << sp;
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
    std::cout << std::endl;

    // Write initial values
    state.outputFile << runTime << ", " << state.gas->temperature() << ", " << state.gas->pressure() << ", " << 1e-6 * Avogadro * state.gas->molarDensity() << ", " << state.gas->meanMolecularWeight();
    std::cout << runTime << ", " << state.gas->temperature() << ", " << state.gas->pressure() << ", " << 1e-6 * Avogadro * state.gas->molarDensity() << ", " << state.gas->meanMolecularWeight();
    for (const auto index : state.indexList) {
        double dens = getNumberDens(state.gas, static_cast<size_t>(index));
        state.outputFile << ", " << dens;
        std::cout << ", " << dens;
    }
    state.outputFile << ", " << 0.0 << ", " << 0.0;
    state.outputFile << std::endl;
    std::cout << std::endl;

    return state;
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

    std::cout << "Reservoir '" << name << "': T=" << config.initialTemperature.value()
              << " K, p=" << config.pressure.value() << " Pa, M=" << state.gas->meanMolecularWeight() << " kg/kmol\n";

    const std::string outputPath = config.outputFile.value_or("");
    if (!outputPath.empty()) {
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
        // Reservoir state is constant: write initial values once
        state.outputFile << runTime << ", " << state.gas->temperature() << ", "
                         << state.gas->pressure() << ", "
                         << 1e-6 * Avogadro * state.gas->molarDensity() << ", " 
                         << state.gas->meanMolecularWeight();
        for (const int idx : state.indexList)
            state.outputFile << ", " << getNumberDens(state.gas, static_cast<size_t>(idx));
        state.outputFile << std::endl;
    }

    return state;
}

void findTimeStep(
    ReactorState& state, ReactorZoneDef& config,
    double runTime,
    double SMALL,
    double pulsePeriod, // VETTORIALIZZARE SE I REATTORI SONO INDIPENDENTI
    double pulseWidth, // VETTORIALIZZARE SE I REATTORI SONO INDIPENDENTI
    double pulseExten, // VETTORIALIZZARE SE I REATTORI SONO INDIPENDENTI
    uint8_t& dischargeOn, uint8_t& fastExpansion, 
    bool stayEN_DC, // VETTORIALIZZARE SE I REATTORI SONO INDIPENDENTI
    const std::string& CP_or_CV, // VETTORIALIZZARE SE I REATTORI SONO INDIPENDENTI
    int iPulse, int nPulses, // VETTORIALIZZARE SE I REATTORI SONO INDIPENDENTI
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
    double EN, double EN_DC, // VETTORIALIZZARE SE I REATTORI SONO INDIPENDENTI
    double& EN_temp, const double EN_TOLERANCE, const double TGAS_TOLERANCE)
{
    // Function for dynamic time stepping

    // Ignition detection and time step adjustment based on dT/dt
    double dt_ig = dt_max; // ignition delay time, to be calculated based on dT/dt
    if (dT > dT_max && dTdt > 1000) {
        dt_ig = std::min(dT_max/dTdt, dt_max); // reduce time step to capture ignition
        std::cout << "!!!!!!!!!!IGNITION DETECTED!!!!!!!!!!" << std::endl;
    }

    EN_temp = std::max(EN_DC, 0.1);  // default E/N value

    // Time step adjustment based on pulse
    if (iPulse < nPulses) { // Still pulsing
        double dt3 = (iPulse + 1) * pulsePeriod - runTime; // Avoid overstepping
        dt_ig = std::min(dt3, dt_ig);

        if ((runTime + SMALL) < (iPulse * pulsePeriod + pulseExten)){ // During pulse or pulse extension
            std::cout << " Reactor " << state.name << " discharging: " << (config.discharge ? "ON" : "OFF") << "; runTime = " << runTime << "; iPulse = " << iPulse << std::endl;
            if ((runTime + SMALL) < (iPulse * pulsePeriod + pulseWidth) && dischargeOn && config.discharge) { // During discharge
                if (iPulse <= 1 || nE_1 < 1e10) {
                    std::cout << "Using fixed time step for first two pulses." << std::endl;
                    dt = dt1;
                } else {
                    std::cout << "Using log-dynamic time step." << std::endl;
                    double K = std::max(std::round(log10(nE_2 / nE_1) * 100) / 100, 0.2);
                    std::cout << "Growth order K = " << K << "; tau_NRP = " << tau_NRP << std::endl;
                    dt = logDynamicTimestep(K, tau_NRP, NstepInPulse, iStepInPulse);
                }

                EN_temp = EN; // Should be an output to update Boltzmann solver

                // PER VETTORIALIZZAZIONE DEGLI IMPULSI
                iStepInPulse += 1;
            } else { // After discharge, during pulse extension
                if (iStepAfPulse == 0) { // First step after the discharge is off; record it
                        // sim.reinitialize();
                        state.integrator->reinitialize(runTime, *state.odes);

                        // PER VETTORIALIZZAZIONE DEGLI IMPULSI
                        tau_NRP = (runTime - iPulse * pulsePeriod); // Duration of the nanosecond discharge
                        std::cout << "NSD terminated at " << disEnergy << " [mJ/cm^3]" << ", tau_NRP = " << tau_NRP << "\n";

                        nE_2 = getNumberDens(state.gas, state.odes->electronIndex); // Store n_e at the end of a nanosecond pulse

                        // TODO: may need to update Boltzmann grid if E/N changes significantly

                        fastExpansion = true;
                    }
                double K = (iPulse <= 1 || nE_2 < 1e10) ? -1 : std::min(round(log10(nE_3 / nE_2) * 100) / 100, -0.1); // fixed K for first two pulses
                std::cout << "Decay order K = " << K << std::endl;
                dt = (iStepAfPulse >= NstepAfPulse) ? // Make dt unchanged for iStepAfPulse >= NstepAfPulse to cut off dt growth
                    dt : logDynamicTimestep(K, pulseExten, NstepAfPulse, iStepAfPulse);
                
                // Isentropic heat loss model. P, T, rho jump at 100 ns after the pulse
                if (runTime > (iPulse * pulsePeriod + 1e-7) && state.odes->heatLoss && fastExpansion)
                {
                    double T_new = state.gas->temperature() * pow( config.pressure.value()/state.gas->pressure(), 1 - state.gas->cv_mass()/state.gas->cp_mass() );
                    std::cout << "Temperature jump: " << state.gas->temperature() << " -> " << T_new << "\n";
                    state.gas->setState_TPY(T_new, config.pressure.value(), state.gas->massFractions());
                    state.odes->setConstPD(state.gas->pressure(), state.gas->density());
                    if (CP_or_CV == "CP"){
                        state.odes->constPressure = true;
                    } else {
                        state.odes->constPressure = false; // Not necessary since heatLoss already changes the state, but just to be safe :)
                    }
                    fastExpansion = false; // only do once for a pulse
                }
                
                // PER VETTORIALIZZAZIONE DEGLI IMPULSI
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

    dt = std::max(dt, 1e-14);  // Make sure dt > 1e-14, avoid dt cascade
    dt = std::min(dt, dt_ig); // Always ensure dt is small enough to capture ignition
    
    if (dt < 1e-16) {
            std::cerr << "dt is too small !" << "\n";
            abort();
        }
    
    
    std::cout << "dt = " << dt << "; iStepInPulse = " << iStepInPulse << "; iStepAfPulse = " << iStepAfPulse
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
    MassFlowSolverSession() {
        Py_Initialize();

        PyObject* sysPath  = PySys_GetObject("path"); // borrowed ref
        PyObject* scriptDir = PyUnicode_FromString(MASSFLOW_SCRIPT_DIR);
        PyList_Insert(sysPath, 0, scriptDir);
        Py_DECREF(scriptDir);

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
    
    CppBOLOS::currentLogLevel = CppBOLOS::LOG_INFO;      // Default log level
    std::string controlDictPath = "../controlDict";       // Default path
    std::string parameterPath = "../chemPlasProperties";  // Default path
    std::string crnDictPath = "../crnDict";        // Default path
    const std::string massFlowPath = "../output/massflow.csv";
    std::ofstream massFlowFile(massFlowPath);

    // Map string to LogLevel
    std::map<std::string, CppBOLOS::LogLevel> logLevels = {
            {"NONE", CppBOLOS::LOG_NONE},
            {"WARNING", CppBOLOS::LOG_WARNING},
            {"INFO", CppBOLOS::LOG_INFO},
            {"DEBUG", CppBOLOS::LOG_DEBUG}
    };

    // Parse command-line arguments
    std::string casePath = "..";
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-case" && i + 1 < argc) {
            casePath = argv[++i];
            controlDictPath = casePath + "/controlDict";
            parameterPath = casePath + "/chemPlasProperties";
            crnDictPath = casePath + "/crnDict";
        } else if (arg == "-log" && i + 1 < argc) {
            std::string logLevelStr = argv[++i];
            auto it = logLevels.find(logLevelStr);
            if (it != logLevels.end()) {
                CppBOLOS::currentLogLevel = it->second;
            } else {
                std::cerr << "Unknown log level: " << logLevelStr << std::endl;
                return 1;
            }
        }
    }

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
    t_end = std::max(pulsePeriod * nPulses, t_end); // overload t_end
    const auto endTimeForce = readParameter<double>(controlDictPath, "endTimeForce", t_end);
    t_end = std::min(t_end, endTimeForce);  // Terminate simulations is endTimeForce is specified

    /* ------------------------------- SET UP BOLTZMANN SOLVER ------------------------------- */
    std::cout << "\n========  SETTING BOLTZMANN SOLVER ... ========\n" << std::endl;

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
    std::cout << "\n========  SETTING GAS PHASE & REACTORS, CREATE  & INIT ODE INTEGRATOR ========\n" << std::endl;
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
    {
        auto csDataFile = readParameter<string>(parameterPath, "csDataFile");
        auto ss = CppBOLOS::clean_file(csDataFile);
        std::vector<CppBOLOS::Collision> collisions = CppBOLOS::parse(ss);
        auto gridSize = readParameter<std::map<std::string, double>>(parameterPath, "gridSize");
        std::string gridType = readParameter<string>(parameterPath, "gridType");
        BoltzmannRate::bsolver.set_grid(gridType, gridSize["start"], gridSize["end"], (int)gridSize["points"]);
        BoltzmannRate::bsolver.load_collisions(collisions);
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

    MassFlowSolverSession pySession;

    /* ------------------------------------- RUN THE SIMULATION ---------------------------------- */
    std::cout << "\n======== RUNNING SIMULATION ========\n" << std::endl;

    clock_t t0 = clock(); // save start time
    double dTdtcheck = dTdt_min * 10; // initial value to enter the loop
    // Main time loop
    // const auto dTbelowTeq = readParameter<double>(controlDictPath, "dTbelowTeq");
    // while (runTime < t_end && gas->temperature() < (T_eq - dTbelowTeq)) {
    while (runTime < t_end && (runTime < 2.0 || dTdtcheck > dTdt_min)) { // Do at least two seconds
        std::ostringstream oss; // Format output digit
        oss << "\nrunTime [s]: " << std::scientific << std::setprecision(9) << runTime;
        // TIMESTEP FUNCTION TO ITERATE FROM HERE

        for (size_t idx = 0; idx < reactors.size(); idx++) {
            ReactorState& state = reactors[idx];
            double& T_old = T_old_vec[idx];
            T_old = state.gas->temperature();
        }

        // HERE WE CALCULATE THE MASS FLOW CONTROLLERS VALUES.
        std::map<std::string, double> yFuel, yOx;
        for (const auto& state : reactors) {
            
            // yFuel[state.name] = state.gas->massFraction(state.gas->speciesIndex(FUELNAME));
            // yOx[state.name] = state.gas->massFraction(state.gas->speciesIndex(OXNAME));
            
            yFuel[state.name] = state.gas->elementalMassFraction(state.gas->elementIndex(FUELNAME));
            yOx[state.name] = state.gas->elementalMassFraction(state.gas->elementIndex(OXNAME));
        }

        for (const auto& state : reservoirs) {
            
            // yFuel[state.name] = state.gas->massFraction(state.gas->speciesIndex(FUELNAME));
            // yOx[state.name] = state.gas->massFraction(state.gas->speciesIndex(OXNAME));
            
            yFuel[state.name] = state.gas->elementalMassFraction(state.gas->elementIndex(FUELNAME));
            yOx[state.name] = state.gas->elementalMassFraction(state.gas->elementIndex(OXNAME));
        }

        std::map<std::string, double> massFlowRates = pySession.call(yOx, yFuel);

        massFlowFile << runTime;
        for (auto& [name, mfc] : massFlowControllers) {
            auto it = massFlowRates.find(name);
            if (it == massFlowRates.end())
                throw std::runtime_error("massFlowSolver: no key '" + name + "' in Python output; check controller names in crnDict match massflowsolver.py");
            mfc.value = it->second;
            massFlowFile << ", " << mfc.value;
            
            ReactorZoneDef& configStart = reactorZones.at(mfc.start);
            shared_ptr<ThermoPhase> srcGas = (configStart.type != "reservoir")
                ? reactors[configStart.indexState].gas
                : reservoirs[configStart.indexState].gas;

            ReactorZoneDef& configEnd = reactorZones.at(mfc.end);
            shared_ptr<ThermoPhase> endGas = (configEnd.type != "reservoir")
                ? reactors[configEnd.indexState].gas
                : reservoirs[configEnd.indexState].gas;

            const size_t nSpec = endGas->nSpecies();
            std::vector<double> Y_src(nSpec, 0.0);
            for (size_t k = 0; k < nSpec; k++) {
                size_t k_src = srcGas->speciesIndex(endGas->speciesName(k));
                Y_src[k] = (k_src != Cantera::npos) ? srcGas->massFraction(k_src) : 0.0;
            }

            mfc.T = srcGas->temperature();
            mfc.comp = Y_src;
            mfc.h = srcGas->enthalpy_mass();
            
        }
        massFlowFile << std::endl;

        /* //REMOVED FOR INTEGRATOR
        // Update density and fractions with mass flow controllers.
        // Two-pass approach: pass 1 reads source states, updates source reactors immediately
        // (density-only removal), and accumulates contributions per destination reactor.
        // Pass 2 applies all accumulated contributions simultaneously, so multiple inflows
        // to the same destination reactor are mixed correctly in a single operation.
        struct EndAccum {
            std::vector<double> sumYmass;             // Σ(Y_src[k] * mass_transferred)
            double sumMass    = 0.0;                  // Σ(mass_transferred)            [kg]
            std::map<std::string, double> sumBoltzXn; // Σ(X_boltz_src * n_transferred)
            double sumN       = 0.0;                  // Σ(n_transferred)               [kmol]
            double sumCpMassT = 0.0;                  // Σ(cp_src * mass_transferred * T_src) [J]
            double sumCpMass  = 0.0;                  // Σ(cp_src * mass_transferred)   [J/K]
        };
        std::map<std::string, EndAccum> endAccum;

        // --- Pass 1: read source states, update source reactors, accumulate for destinations ---
        for (auto& [name, mfc] : massFlowControllers) {
            double mass_transferred = mfc.value * dt; // [kg]
            if (mass_transferred < 0.0)
                std::cerr << "[Warning] Negative mass transfer for controller '" << name
                          << "' (" << mass_transferred << " kg). Check massflowsolver sign convention.\n";

            ReactorZoneDef& configStart = reactorZones.at(mfc.start);
            ReactorZoneDef& configEnd   = reactorZones.at(mfc.end);

            double T_src = 0.0, cp_src = 0.0, n_transferred = 0.0;

            if (configStart.type != "reservoir") {
                ReactorState& stateStart = reactors[configStart.indexState];
                double nStart     = stateStart.gas->molarDensity() * configStart.volume; // [kmol]
                double mass_start = nStart * stateStart.gas->meanMolecularWeight();                // [kg]
                T_src  = stateStart.gas->temperature();
                cp_src = stateStart.gas->cp_mass();

                if (configEnd.type == "reservoir") {
                    mass_transferred += kOutPressure * (stateStart.gas->pressure() - configEnd.pressure.value()) * dt;
                }
                n_transferred    = mass_transferred / stateStart.gas->meanMolecularWeight();
                n_transferred    = std::min(n_transferred, nStart);
                if (n_transferred == nStart) {
                    std::cerr << "[Warning] Mass transfer from '" << mfc.start << "' fully depletes the source reactor. Consider reducing the mass flow rate or increasing the time resolution.\n";
                }
                mass_transferred = n_transferred * stateStart.gas->meanMolecularWeight();

                // Accumulate into destination only if it is a reactor (not a reservoir)
                if (configEnd.type != "reservoir") {
                    const ReactorState& stateEnd = reactors[configEnd.indexState];
                    const size_t nSpec = stateEnd.gas->nSpecies();
                    EndAccum& accum = endAccum[mfc.end];
                    if (accum.sumYmass.empty()) accum.sumYmass.assign(nSpec, 0.0);

                    for (size_t k = 0; k < nSpec; k++) {
                        size_t k_src = stateStart.gas->speciesIndex(stateEnd.gas->speciesName(k));
                        accum.sumYmass[k] += ((k_src != Cantera::npos) ? stateStart.gas->massFraction(k_src) : 0.0)
                                            * mass_transferred;
                    }
                    for (const auto& [sp, X_src] : *stateStart.boltzmannSpecies)
                        accum.sumBoltzXn[sp] += X_src * n_transferred;
                    accum.sumMass    += mass_transferred;
                    accum.sumN       += n_transferred;
                    accum.sumCpMassT += cp_src * mass_transferred * (T_src - stateEnd.gas->temperature());
                    accum.sumCpMass  += cp_src * mass_transferred;
                }
                //
                //const size_t nSpec = stateStart.gas->nSpecies();
                //EndAccum& accum = endAccum[mfc.start];
                //if (accum.sumYmass.empty()) accum.sumYmass.assign(nSpec, 0.0);
                //
                //for (size_t k = 0; k < nSpec; k++) {
                //    size_t k_src = stateStart.gas->speciesIndex(stateStart.gas->speciesName(k));
                //    accum.sumYmass[k] -= ((k_src != Cantera::npos) ? stateStart.gas->massFraction(k_src) : 0.0)
                //                        * mass_transferred;
                //}
                //for (const auto& [sp, X_src] : *stateStart.boltzmannSpecies)
                //    accum.sumBoltzXn[sp] -= X_src * n_transferred;
                //accum.sumMass    -= mass_transferred;
                //accum.sumN       -= n_transferred;
                //accum.sumCpMassT -= cp_src * mass_transferred * T_src;
                //accum.sumCpMass  -= cp_src * mass_transferred;
                // Always update source density (regardless of destination type)
                double new_density_start = (mass_start - mass_transferred) / configStart.volume;
                stateStart.gas->setState_TD(T_src, new_density_start);
                stateStart.odes->setConstPD(stateStart.gas->pressure(), new_density_start);
                stateStart.integrator->reinitialize(runTime, *stateStart.odes);
            } else {

                if (configEnd.type != "reservoir") {
                    // Reservoir source: infinite supply at fixed state, no density update
                    const ReservoirState& stateStart = reservoirs[configStart.indexState];
                    T_src  = stateStart.gas->temperature();
                    cp_src = stateStart.gas->cp_mass();
                    n_transferred = mass_transferred / stateStart.gas->meanMolecularWeight(); // no clamping

                    const ReactorState& stateEnd = reactors[configEnd.indexState];
                    const size_t nSpec = stateEnd.gas->nSpecies();
                    EndAccum& accum = endAccum[mfc.end];
                    if (accum.sumYmass.empty()) accum.sumYmass.assign(nSpec, 0.0);

                    for (size_t k = 0; k < nSpec; k++) {
                        size_t k_src = stateStart.gas->speciesIndex(stateEnd.gas->speciesName(k));
                        accum.sumYmass[k] += ((k_src != Cantera::npos) ? stateStart.gas->massFraction(k_src) : 0.0)
                                            * mass_transferred;
                    }
                    for (const auto& [sp, X_src] : *stateStart.boltzmannSpecies)
                        accum.sumBoltzXn[sp] += X_src * n_transferred;
                    accum.sumMass    += mass_transferred;
                    accum.sumN       += n_transferred;
                    accum.sumCpMassT += cp_src * mass_transferred * T_src;
                    accum.sumCpMass  += cp_src * mass_transferred;
                }
            }
        }

        // --- Pass 2: apply all accumulated contributions to each reactor ---
        for (auto& [endName, accum] : endAccum) {
            if (accum.sumMass <= 0.0) continue;

            ReactorZoneDef& configEnd = reactorZones.at(endName);
            ReactorState& stateEnd    = reactors[configEnd.indexState];
            const size_t nSpec        = stateEnd.gas->nSpecies();

            double nEnd     = stateEnd.gas->molarDensity() * configEnd.volume; // [kmol]
            double mass_end = nEnd * stateEnd.gas->meanMolecularWeight();                // [kg]
            double T_end    = stateEnd.gas->temperature();
            double cp_end   = stateEnd.gas->cp_mass();

            // Mix Boltzmann species mole fractions (mole-weighted)
            for (auto& [sp, X_end] : *stateEnd.boltzmannSpecies) {
                double sumX = accum.sumBoltzXn.count(sp) ? accum.sumBoltzXn.at(sp) : 0.0;
                X_end = (X_end * nEnd + sumX) / (nEnd + accum.sumN);
            }

            // Mix temperature via enthalpy balance
            T_end = (cp_end * mass_end * T_end + accum.sumCpMassT)
                  / (cp_end * mass_end + accum.sumCpMass);

            // Mix all K species by mass fraction
            double mass_total = mass_end + accum.sumMass;
            std::vector<double> Y_new(nSpec);
            for (size_t k = 0; k < nSpec; k++)
                if (mass_total > 0.0) {
                    Y_new[k] = (stateEnd.gas->massFraction(k) * mass_end + accum.sumYmass[k]) / mass_total;
                } else {
                    mass_total = SMALL; // prevent negative or zero mass_total due to numerical issues
                    Y_new[k] = stateEnd.gas->massFraction(k); // fallback to old composition if mass_total is zero or negative (should not happen)
                }
            double new_density_end = mass_total / configEnd.volume;
            stateEnd.gas->setMassFractions_NoNorm(Y_new.data());
            stateEnd.gas->setState_TD(T_end, new_density_end);
            stateEnd.odes->setConstPD(stateEnd.gas->pressure(), new_density_end);
            stateEnd.boltzmannState.density = *stateEnd.boltzmannSpecies;
            stateEnd.integrator->reinitialize(runTime, *stateEnd.odes);
        }
        // REMOVED FOR INTEGRATOR */
    
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

            std::cout << oss.str() << " | reactor " << state.name << ", iPulse: " << iPulse << "\n";
            /*
            {
                BoltzmannSolverGuard guard(state.boltzmannState);

                findTimeStep(state, config, runTime, SMALL, pulsePeriod, pulseWidth, pulseExten, dischargeOn, fastExpansion, stayEN_DC,
                            CP_or_CV, iPulse, nPulses, dt_element, dt1, dt_max, dT, dT_max, dTdt,
                            nE_1, nE_2, nE_3, tau_NRP, disEnergy,
                            iStepInPulse, NstepInPulse,
                            iStepAfPulse, NstepAfPulse,
                            EN, EN_DC, EN_temp, EN_TOLERANCE, TGAS_TOLERANCE);
            }

            state.boltzmannState.density = *state.boltzmannSpecies;
            */

            findTimeStep(state, config, runTime, SMALL, pulsePeriod, pulseWidth, pulseExten, dischargeOn, fastExpansion, stayEN_DC,
                        CP_or_CV, iPulse, nPulses, dt_element, dt1, dt_max, dT, dT_max, dTdt,
                        nE_1, nE_2, nE_3, tau_NRP, disEnergy,
                        iStepInPulse, NstepInPulse,
                        iStepAfPulse, NstepAfPulse,
                        EN, EN_DC, EN_temp, EN_TOLERANCE, TGAS_TOLERANCE);

            dt = std::min(dt, dt_element); // Take the minimum dt among reactors to ensure stability
        }
        

        /*----------------------  Advance time step --------------------*/
        const double runTimeOld = runTime;
        runTime += dt;

        
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
            double& EN_temp = EN_vec[idx];

            // T_old = state.gas->temperature();

            {
                BoltzmannSolverGuard guard(state.boltzmannState);
                
                //
                // Update Boltzmann solver as needed and advance time step.
                // This could be costly and should be optimized
                bool updateBoltzmannSolver = false;
                if (state.odes->updateBoltzmannMixture()) {
                    BoltzmannRate::bsolver.set_density(*state.boltzmannSpecies);
                    updateBoltzmannSolver = true;
                }
                if (std::abs(EN_temp - BoltzmannRate::bsolver.get_EN()) > EN_TOLERANCE) {
                    BoltzmannRate::bsolver.set_EN(EN_temp);
                    updateBoltzmannSolver = true;
                }
                if(std::abs(state.gas->temperature() - BoltzmannRate::bsolver.get_kT()) > TGAS_TOLERANCE) {
                    BoltzmannRate::bsolver.set_kT(state.gas->temperature());
                    updateBoltzmannSolver = true;
                }
                if (updateBoltzmannSolver) {

                    {
                        BoltzmannRate::bsolver.init();

                        std::cout << "Updating EEDF ..." << "\n";
                        if (EN_temp < 1) {
                            BoltzmannRate::F0 = BoltzmannRate::bsolver.maxwell(state.gas->temperature()*8.6173E-5);
                        }

                        BoltzmannRate::updateBoltzmannSolver(200, 1e-5, 1E20/(EN_temp*EN_temp));

                        state.integrator->reinitialize(runTimeOld, *state.odes);
                    }
                }

                if (iPulse == 0) {
                    state.integrator->reinitialize(runTimeOld, *state.odes);
                }

                std::cout << "EN = " << BoltzmannRate::bsolver.get_EN()
                            << "; Tgas = " << state.gas->temperature()
                            << "; Te = " << BoltzmannRate::bsolver.get_Te()
                            << "; n_e = " << getNumberDens(state.gas, state.odes->electronIndex)
                            << "\n";
                //

                // NEW PART TO INTEGRATE MASS FLOW CONTROLLERS
                std::vector<ChemPlasReactor::FlowTerm> inflows;
                double mdot_in = 0.0;
                double mdot_out = 0.0;
                for (auto& [name, mfc]: massFlowControllers) {
                    if (mfc.end == state.name) {
                        inflows.push_back({mfc.value / config.volume, mfc.comp, mfc.T, mfc.h});
                        mdot_in += mfc.value;
                    }
                    if (mfc.start == state.name) {
                        mdot_out += mfc.value;
                    }
                }
                state.odes->setFlows(inflows, mdot_out / config.volume);
                state.integrator->reinitialize(runTimeOld, *state.odes);
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                state.integrator->integrate(runTime);
                state.odes->updateState(state.integrator->solution());

                
                // For CV reactors: update m_density to account for net mass flux.
                // For CP reactors: density follows the EOS automatically — no update needed.
                if (!state.odes->constPressure) {
                    double net_mdot_V = (mdot_in - mdot_out) / config.volume;
                    double rho_new = state.gas->density() + net_mdot_V * dt;
                    state.gas->setState_TD(state.gas->temperature(), std::max(rho_new, 1e-6));
                    state.odes->setConstPD(state.gas->pressure(), std::max(rho_new, 1e-6));
                }
                

                dT = std::abs(state.gas->temperature() - T_old);

                // Write values
                state.outputFile << runTime << ", " << state.gas->temperature() << ", " << state.gas->pressure() << ", " << 1e-6 * Avogadro * state.gas->molarDensity() << ", " << state.gas->meanMolecularWeight();
                // std::cout << runTime << ", " << state.gas->temperature() << ", " << state.gas->pressure() << ", " << 1e-6 * Avogadro * state.gas->molarDensity() << ", " << state.gas->meanMolecularWeight();
                for (const auto index : state.indexList) {
                    double dens = getNumberDens(state.gas, static_cast<size_t>(index));
                    state.outputFile << ", " << dens;
                    // std::cout << ", " << dens;
                }
                
                std::cout << std::endl;

                // Check plasma power and energy deposited
                double power_temp = BoltzmannRate::bsolver.elec_power()
                                    * Avogadro * state.gas->molarDensity()* ElectronCharge * getNumberDens(state.gas, state.odes->electronIndex);
                // unit: eV m^3/s * #/kmol * kmol/m^3 * coulomb(J/eV) * #/cm^3 = J/s/cm^3
                disEnergy = 1e-3 * state.odes->depositedPlasmaEnergy(); // mJ/cm^3
                
                state.outputFile << ", " << power_temp << ", " << disEnergy;
                state.outputFile << std::endl;
                std::cout << "Power_input [W/cm^3]: " << power_temp << " | Energy deposited [mJ/cm^3]: " << disEnergy << "\n";

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
                if (runTime + SMALL >= (iPulse + 1) * pulsePeriod ) { // Entering next pulse period
                    iPulse += 1;
                    if (iPulse == nPulses) {
                        std::cout << "\n >>>>>>>> NRP Discharges Finished for reactor " << state.name << " <<<<<<<< \n" << std::endl;
                        state.integrator->reinitialize(runTime, *state.odes);
                        plasmaOn = false;
                    }
                    if (plasmaOn) {
                        std::cout << "\n**** Start Pulse Loop No. " << iPulse + 1 << " at runTime = " << runTime << " for reactor " << state.name << " ****\n";
                        if (state.odes->heatLoss) {
                            state.odes->setConstPD(state.gas->pressure(), state.gas->density());
                            state.odes->constPressure = false; // switch to constant volume
                        }

                        disEnergy = 0.0, dischargeOn = true;
                        state.odes->resetDepositedPlasmaEnergy();
                        iStepInPulse = 0, iStepAfPulse = 0;
                        nE_1 = getNumberDens(state.gas, state.odes->electronIndex); // n_e at the beginning of a new pulse
                    }
                }
            }
            state.boltzmannState.density = *state.boltzmannSpecies;
        }

        dTdtcheck = *std::max_element(dTdt_vec.begin(), dTdt_vec.end());
        
        std::cout << "-------------" << "\n";
    }
    // End of the main time loop

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
    std::cout << "\n======== ChemPlasKin finished. Execution time: " << elapsed_time << " (sec) ========\n" <<std::endl;

    return 0;
}