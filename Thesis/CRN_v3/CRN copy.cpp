/*--------------------------------*- C++ -*----------------------------------*\
|     ____ _                    ____  _           _  ___                      |
|    / ___| |__   ___ _ __ ___ |  _ \| | __ _ ___| |/ (_)_ __                 |
|   | |   | '_ \ / _ \ '_ ` _ \| |_) | |/ _` / __| ' /| | '_ \                |
|   | |___| | | |  __/ | | | | |  __/| | (_| \__ \ . \| | | | |               |
|    \____|_| |_|\___|_| |_| |_|_|   |_|\__,_|___/_|\_\_|_| |_|               |
|                                                                             |
|   A Freeware for Unified Gas-Plasma Kinetics Simulation                     |
|   Version:      1.2 (Feb 2026)                                           |
|   License:      GNU LESSER GENERAL PUBLIC LICENSE, Version 2.1              |
|   Author:       Xiao Shao                                                   |
|   Organization: King Abdullah University of Science and Technology (KAUST)  |
|   Contact:      xiao.shao@kaust.edu.sa                                      |
\*---------------------------------------------------------------------------*/

// Test case of air plasma kinetics using external profiles of electron density and electric field
// Reference: https://doi.org/10.1088/0022-3727/46/46/464010; https://doi.org/10.1016/j.combustflame.2022.111990

#include "cantera/core.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/ext/bolos/Logger.h"
#include "cantera/kinetics/Boltzmann.h"
#include "cantera/zerodim.h"
#include "cantera/numerics/Integrator.h"

using namespace Cantera;
#include "plasmaReactor.h" // "../../src/plasmaReactor.h" for original one
#include "../../src/utilities.h"

struct ReactorConfig {
    std::string name;
    std::string mechanismFile;
    std::string outputFile;
    double initialTemperature;
    double pressure;
    double initialElectronDensity;
    Composition boltzmannSpecies;
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
};

ReactorState createReactorState(
    const ReactorConfig& config,
    const std::string& plasmaHeatModel,
    bool thermalEffect,
    const std::string& inertSpecies)
{
    ReactorState state;
    state.name = config.name;
    state.boltzmannSpecies = std::make_shared<std::map<std::string, double>>(config.boltzmannSpecies);
    state.sol = newSolution(config.mechanismFile);
    state.gas = state.sol->thermo();

    const double gasNumberDensity = 1e-6 * config.pressure / config.initialTemperature / CppBOLOS::KB;
    Composition compMap = config.boltzmannSpecies;
    compMap["Electron"] = config.initialElectronDensity / gasNumberDensity;
    state.gas->setMoleFractionsByName(compMap);
    state.gas->setState_TP(config.initialTemperature, config.pressure);

    state.odes = std::make_unique<ChemPlasReactor>(state.sol, *state.boltzmannSpecies, plasmaHeatModel);
    state.odes->inertSpIndex = state.odes->findSpeciesIndex(inertSpecies);
    state.odes->nonThermal = !thermalEffect;
    state.odes->constPressure = false;

    state.integrator = std::unique_ptr<Integrator>(newIntegrator("CVODE"));
    state.integrator->initialize(0.0, *state.odes);

    auto kin = state.sol->kinetics();
    const int irxns = kin->nReactions();
    state.qf.resize(irxns);
    state.qr.resize(irxns);
    state.q.resize(irxns);

    state.outputFile.open(config.outputFile);
    if (!state.outputFile.is_open()) {
        throw std::runtime_error("Failed to open output file for reactor: " + config.name);
    }

    const std::vector<std::string> speciesNames = state.gas->speciesNames();
    state.outputFile << "Time(s), T_gas(K), N_gas(#/cm^3)";
    std::cout << "Writing information of " << config.name << ":\n" << "Time, T_gas(K), N_gas(#/cm^3)";
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
    state.outputFile << std::endl;
    std::cout << std::endl;

    return state;
}

void advanceReactorState(
    ReactorState& state,
    double runTime,
    double dt,
    const std::vector<std::pair<double, double>>& Vp_t_data,
    const std::vector<std::pair<double, double>>& Ne_t_data,
    bool printReactionRates)
{
    auto gas = state.gas;
    auto kin = state.sol->kinetics();

    double EN = std::max(interpolate(Vp_t_data, runTime * 1e9) * 1000 / (4e-3 * Avogadro * gas->molarDensity()) * 1e21, 0.1);
    double nE = interpolate(Ne_t_data, runTime * 1e9);

    BoltzmannRate::bsolver.set_kT(gas->temperature());
    BoltzmannRate::bsolver.set_EN(EN);
    if (state.odes->updateBoltzmannMixture()) {
        BoltzmannRate::bsolver.set_density(*state.boltzmannSpecies);
    }

    BoltzmannRate::bsolver.init();
    BoltzmannRate::updateBoltzmannSolver(200, 1e-5, 1E20 / (EN * EN));

    state.odes->imposeNe(nE);
    state.integrator->reinitialize(runTime, *state.odes);

    const double nextTime = runTime + dt;
    state.integrator->integrate(nextTime);
    state.odes->updateState(state.integrator->solution());
    state.odes->nSteps = 0;

    state.outputFile << nextTime << ", " << gas->temperature() << ", " << 1e-6 * Avogadro * gas->molarDensity();
    for (const auto index : state.indexList) {
        state.outputFile << ", " << getNumberDens(gas, static_cast<size_t>(index));
    }
    state.outputFile << std::endl;

    const double power_temp = BoltzmannRate::bsolver.elec_power()
        * Avogadro * gas->molarDensity() * ElectronCharge * getNumberDens(gas, state.odes->electronIndex);
    const double disEnergy = 1e-3 * state.odes->depositedPlasmaEnergy();

    std::cout << "[" << state.name << "] EN = " << EN
              << ", Tgas = " << gas->temperature()
              << "[K]; Te = " << BoltzmannRate::bsolver.get_Te()
              << "[K]; power = " << power_temp
              << " W/cm^3; deposited = " << disEnergy << " mJ/cm^3\n";

    if (printReactionRates) {
        const int irxns = kin->nReactions();
        kin->getFwdRatesOfProgress(&state.qf[0]);
        kin->getRevRatesOfProgress(&state.qr[0]);
        kin->getNetRatesOfProgress(&state.q[0]);

        writelog("{:30s} {:>14s} {:>14s} {:>14s}  {:s}\n",
                 "Reaction", "Forward", "Reverse", "Net", "Unit");
        for (int i = 0; i < irxns; i++) {
            const auto& rxn = kin->reaction(i);
            writelog("{:30s} {:14.5g} {:14.5g} {:14.5g}  kmol/m3/s\n",
                     rxn->equation(), state.qf[i], state.qr[i], state.q[i]);
        }
    }
}

int main()
{
    CppBOLOS::currentLogLevel = CppBOLOS::LOG_INFO;

    // Set up time control and pulse number
    double runTime = 0.0;
    double dt = 0.1E-9;
    double endTime = 100E-9;
    bool thermalEffect = true;
    bool printReactionRates = false;

    /* ------------------ READ EXPERIMENTAL DATA ------------------ */
    std::vector<std::pair<double, double>> Vp_t_data;
    std::vector<std::pair<double, double>> Ne_t_data;
    // Read data from CSV
    readCSV("../Vp-t.csv", Vp_t_data);
    readCSV("../ne-t.csv", Ne_t_data);

    // Example usage of interpolate function
    double queryPoint = 5.0;  // Example query point
    try {
        double interpolatedVp = interpolate(Vp_t_data, queryPoint); // linearinterpolate(Vp_t_data, queryPoint); for linear interpolation
        double interpolatedNe = interpolate(Ne_t_data, queryPoint); // linearinterpolate(Vp_t_data, queryPoint); for linear interpolation
        std::cout << "Interpolated Vp(kV) value at t = " << queryPoint << " is: " << interpolatedVp << std::endl;
        std::cout << "Interpolated Ne(#/cm^3) value at t = " << queryPoint << " is: " << interpolatedNe << std::endl;
    } catch (const std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
    }

    /* ------------------------------- SET UP BOLTZMANN SOLVER ------------------------------- */
    std::cout << "\n========  SETTING BOLTZMANN SOLVER ... ========\n" << std::endl;    // Species configured by CppBOLOS

    // Species configured by CppBOLOS
    std::map<std::string, double> BoltzmannSpecies = {
            // Be careful of possible difference of species names in LXCat and *yaml input,
            // especially for e/E, He/HE, Ar/AR. Always try to use unified names.
            {"N2", 0.774},
            {"O2", 0.186},
            {"O", 0.04},
            {"N2(A3)", 0}, {"N2(B3)", 0}, {"N2(a1)", 0}, {"N2(C3)", 0},
            {"O2(a1)", 0}, {"O(1D)", 0},
            {"N", 0}, {"N(2D)", 0},
            {"NO", 0},
            {"O3", 0},
            {"H2", 0}, {"H", 0},
            {"H2O", 0},
//            {"CH4", 0},
//            {"Ar", 0},
    };

    // Read cross-section data
    std::string CS_data_file = "../../../data/LXCat/bolsigdb_air_NH3_H2.dat";
    const auto ss = CppBOLOS::clean_file(CS_data_file);
    std::vector<CppBOLOS::Collision> collisions = CppBOLOS::parse(ss); // parse collision data

    // Set up grid. This can affect accuracy.
    BoltzmannRate::bsolver.set_grid("QuadraticGrid", 0, 60, 150);
    BoltzmannRate::bsolver.load_collisions(collisions);
    LOG_INFO("\nA total of " + std::to_string(BoltzmannRate::bsolver.number_of_targets()) +
             " targets have been loaded:\n" + BoltzmannRate::bsolver.targetNames());

    // Set T_gas, E/N, density and Initialization
    double Tgas = 1500; // gas teperature [K]
    double EN = 150; // reduced electirc field [Td]
    double nE = 1e13; // initial electron number density [#/cm^3]
    BoltzmannRate::bsolver.set_kT (Tgas);
    BoltzmannRate::bsolver.set_EN(EN);
    BoltzmannRate::bsolver.set_density(BoltzmannSpecies);
    BoltzmannRate::bsolver.init();

    // Solve EBE, serve as cache and preliminary convergence check
    BoltzmannRate::F0 = BoltzmannRate::bsolver.maxwell(4.0); // initial guess from Maxwell EEDF

    try{
        BoltzmannRate::F0 = BoltzmannRate::bsolver.converge(BoltzmannRate::F0, 200, 1e-5);
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    double mean_energy = BoltzmannRate::bsolver.mean_energy(BoltzmannRate::F0);
    std::cout << "mean energy: " << mean_energy << "eV (" << BoltzmannRate::bsolver.get_Te() << "K)" << std::endl;

    /* ----------------------------- SET UP GAS PHASE AND REACTORS  -------------------------------- */
    std::cout << "\n========  SETTING GAS PHASE & REACTORS ... ========\n" << std::endl;

    std::cout << "counting Boltzmann Processes: " << BoltzmannRate::NumProcess << std::endl;

    const double P0 = OneAtm;
    const std::string mechanismFile = "../../../data/PAC_kinetics/gri30_plasma.yaml";
    const double gas_density = 1e-6 * P0 / Tgas / CppBOLOS::KB; // [#/cm^3]

    std::vector<ReactorConfig> reactorConfigs = {
        {
            "zone_1",
            mechanismFile,
            "../outputAirCRN_zone1.csv",
            1500.0,
            P0,
            1e13,
            Composition{
                {"N2", 0.774},
                {"O2", 0.186},
                {"O", 0.04},
                {"N2(A3)", 0}, {"N2(B3)", 0}, {"N2(a1)", 0}, {"N2(C3)", 0},
                {"O2(a1)", 0}, {"O(1D)", 0},
                {"N", 0}, {"N(2D)", 0},
                {"NO", 0},
                {"O3", 0},
                {"H2", 0.0}, {"H", 0.0}, {"H2O", 0.0}
            }
        },
        {
            "zone_2",
            mechanismFile,
            "../outputAirCRN_zone2.csv",
            1650.0,
            P0,
            8e12,
            Composition{
                {"N2", 0.760},
                {"O2", 0.190},
                {"O", 0.030},
                {"N2(A3)", 0}, {"N2(B3)", 0}, {"N2(a1)", 0}, {"N2(C3)", 0},
                {"O2(a1)", 0}, {"O(1D)", 0},
                {"N", 0}, {"N(2D)", 0},
                {"NO", 0},
                {"O3", 0},
                {"H2", 0.010}, {"H", 0.0}, {"H2O", 0.0}
            }
        }
    };

    std::vector<ReactorState> reactors;
    reactors.reserve(reactorConfigs.size());

    for (const auto& config : reactorConfigs) {
        reactors.push_back(createReactorState(config, " ", thermalEffect, "N2"));
    }

    std::vector<std::string> species_names = reactors.front().gas->speciesNames();

    for (auto& reactor : reactors) {
        auto gas = reactor.gas;
        std::cout << "[" << reactor.name << "] initial density: "
                  << getNumberDens(gas, gas->speciesIndex("Electron"))
                  << ", N2+ mole fraction: " << gas->moleFraction("N2+")
                  << ", nSpecies: " << gas->nSpecies() << std::endl;
    }

    std::cout << "Time\t";
    for (const auto& name : species_names) {
        std::cout << name << "\t";
    }
    std::cout << std::endl;

    /* ------------------------------------- RUN THE SIMULATION ---------------------------------- */
    std::cout << "\n======== RUNNING SIMULATION ========\n" << std::endl;

    clock_t t0 = clock(); // save start time

    // Main time loop
    while (runTime < endTime){

        // Relax dt after 25ns
        if (runTime > 25e-9) {
            dt = 5e-10;
        }

        std::cout << "\nrunTime [s]: "  << runTime << ", dt = " << dt << std::endl;

        // Each reactor advances in the same time window; the Boltzmann solver is reused sequentially.
        try {
            for (auto& reactor : reactors) {
                advanceReactorState(reactor, runTime, dt, Vp_t_data, Ne_t_data, printReactionRates);
            }
        } catch (const std::runtime_error& e) {
            std::cerr << e.what() << std::endl;
            return -1;
        }

        runTime += dt;

        std::cout << "-------------" << "\n";
    }

    clock_t t1 = clock(); // save end time
    double elapsed_time = static_cast<double>(t1 - t0) / CLOCKS_PER_SEC;
    std::cout << "\n======== ChemPlasKin finished. Execution time: " << elapsed_time << " (sec) ========\n" <<std::endl;

    return 0;

}
