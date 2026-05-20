#pragma once
#include "cantera/kinetics/Boltzmann.h"
#include <map>
#include <string>

/**
 * Snapshot dello stato numerico del BoltzmannSolver.
 * Copiabile perché contiene solo scalari, mappe e VectorXd.
 */
struct BoltzmannSnapshot {
    double EN   = 0.0;
    double kT   = 0.0;
    std::map<std::string, double> density = {};
    Eigen::VectorXd F0;

    /// Cattura lo stato attuale dal solver globale
    static BoltzmannSnapshot capture() {
        BoltzmannSnapshot s;
        s.EN      = BoltzmannRate::bsolver.get_EN();
        s.kT      = BoltzmannRate::bsolver.get_kT();
        s.F0      = BoltzmannRate::F0;
        return s;
    }

    /// Installa questo snapshot nel solver globale e lo reinizializza
    void install() const {
        BoltzmannRate::bsolver.set_EN(EN);
        BoltzmannRate::bsolver.set_kT(kT);
        BoltzmannRate::bsolver.set_density(density);
        BoltzmannRate::bsolver.init();
        BoltzmannRate::F0 = F0;
        // Signal all MultiRate caches to recompute on next eval() call,
        // since the global Boltzmann state has changed.
        BoltzmannRate::updateOnce = true;
    }
};

/**
 * RAII guard: installa lo snapshot di un reattore nel solver globale
 * e ne salva lo stato aggiornato alla distruzione.
 *
 * Non ripristina uno stato precedente: lo stato globale tra due guard
 * non è mai usato, quindi non c'è nulla da ripristinare. Ogni reattore
 * installa il proprio stato tramite install() all'inizio del guard.
 *
 * Nota: BoltzmannSnapshot::capture() non cattura density, quindi un
 * backup catturato sarebbe sempre con density={} e chiamerebbe
 * set_density({}) + init() corrompendo le tabelle interne di CppBOLOS.
 *
 * Uso:
 *   {
 *     BoltzmannSolverGuard guard(state.snapshot);
 *     state.integrator->integrate(runTime);
 *   } // stato del reattore aggiornato automaticamente
 */
class BoltzmannSolverGuard {
public:
    explicit BoltzmannSolverGuard(BoltzmannSnapshot& reactorSnapshot)
        : m_reactorSnapshot(reactorSnapshot)
    {
        reactorSnapshot.install();
    }

    ~BoltzmannSolverGuard() {
        m_reactorSnapshot.EN = BoltzmannRate::bsolver.get_EN();
        m_reactorSnapshot.kT = BoltzmannRate::bsolver.get_kT();
        m_reactorSnapshot.F0 = BoltzmannRate::F0;
    }

    BoltzmannSolverGuard(const BoltzmannSolverGuard&) = delete;
    BoltzmannSolverGuard& operator=(const BoltzmannSolverGuard&) = delete;

private:
    BoltzmannSnapshot& m_reactorSnapshot;
};