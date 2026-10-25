// CRN_v8: CVODES integrator with a user-supplied Jacobian.
//
// Same set-up as Cantera's CVodesIntegrator with the "DENSE" linear solver (BDF, Newton, dense LAPACK
// solver, scalar tolerances, at most 20000 steps per integrate(), linear solver rebuilt at every restart),
// so that without a Jacobian (mode Fd) it takes exactly the same steps. In mode Analytic CVODE calls
// ChemPlasReactor::jacobian() instead of building the Jacobian by finite differences (N right-hand-side
// evaluations per Jacobian). With checkEvery > 0 every checkEvery-th analytic Jacobian is compared with
// a finite-difference one (diagnostics only, slow).
#ifndef CRN_CVODE_JAC_H
#define CRN_CVODE_JAC_H

#include "cantera/numerics/Integrator.h"
#include "cantera/numerics/sundials_headers.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/AnyMap.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "plasmaReactor.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace CRN
{

class CVodeJac : public Cantera::Integrator
{
public:
    enum class Mode { Fd, Analytic };

    explicit CVodeJac(Mode mode, int checkEvery = 0) : m_mode(mode), m_checkEvery(checkEvery) {}

    ~CVodeJac() override {
        if (m_mem) CVodeFree(&m_mem);
        if (m_ls) SUNLinSolFree(m_ls);
        if (m_A) SUNMatDestroy(m_A);
        if (m_y) N_VDestroy_Serial(m_y);
        if (m_dky) N_VDestroy_Serial(m_dky);
    }

    void setTolerances(double reltol, double abstol) override {
        m_reltol = reltol;
        m_abstol = abstol;
    }

    void initialize(double t0, Cantera::FuncEval& func) override {
        m_func = &func;
        m_reactor = dynamic_cast<Cantera::ChemPlasReactor*>(&func);
        m_neq = func.neq();
        func.clearErrors();
        if (m_y) N_VDestroy_Serial(m_y);
        if (m_dky) N_VDestroy_Serial(m_dky);
        m_y = N_VNew_Serial(static_cast<sunindextype>(m_neq));
        m_dky = N_VNew_Serial(static_cast<sunindextype>(m_neq));
        N_VConst(0.0, m_y);
        N_VConst(0.0, m_dky);
        func.getState(NV_DATA_S(m_y));
        m_time = m_tInteg = t0;
        if (m_mem) CVodeFree(&m_mem);
        m_mem = CVodeCreate(CV_BDF);
        if (!m_mem) throw Cantera::CanteraError("CVodeJac::initialize", "CVodeCreate failed.");
        check(CVodeInit(m_mem, rhs, t0, m_y), "CVodeInit");
        CVodeSetErrHandlerFn(m_mem, errHandler, this);
        check(CVodeSStolerances(m_mem, m_reltol, m_abstol), "CVodeSStolerances");
        check(CVodeSetUserData(m_mem, this), "CVodeSetUserData");
        applyOptions();
    }

    void reinitialize(double t0, Cantera::FuncEval& func) override {
        accumulateStats(); // CVodeReInit resets the counters
        m_func = &func;
        m_reactor = dynamic_cast<Cantera::ChemPlasReactor*>(&func);
        m_time = m_tInteg = t0;
        func.getState(NV_DATA_S(m_y));
        func.clearErrors();
        check(CVodeReInit(m_mem, t0, m_y), "CVodeReInit");
        applyOptions();
    }

    void integrate(double tout) override {
        if (tout == m_time) return;
        if (tout < m_time) {
            throw Cantera::CanteraError("CVodeJac::integrate", "Cannot integrate backwards in time: {} < {}", tout, m_time);
        }
        int nsteps = 0;
        while (m_tInteg < tout) {
            if (nsteps >= m_maxsteps) {
                throw Cantera::CanteraError("CVodeJac::integrate", "Maximum number of timesteps ({}) taken "
                                            "without reaching output time ({}).", nsteps, tout);
            }
            cvodeStep(tout, "integrate");
            nsteps++;
        }
        check(CVodeGetDky(m_mem, tout, 0, m_y), "CVodeGetDky");
        m_time = tout;
    }

    double step(double tout) override {
        cvodeStep(tout, "step");
        m_time = m_tInteg;
        return m_time;
    }

    double* solution() override { return NV_DATA_S(m_y); }
    double& solution(size_t k) override { return NV_DATA_S(m_y)[k]; }

    double* derivative(double tout, int n) override {
        check(CVodeGetDky(m_mem, tout, n, m_dky), "CVodeGetDky");
        return NV_DATA_S(m_dky);
    }

    int nEquations() const override { return static_cast<int>(m_neq); }

    // Counters since the last restart, plus "total_*": since the start of the run
    Cantera::AnyMap solverStats() const override {
        long steps = 0, rhs = 0, jac = 0, linRhs = 0, nli = 0, ncf = 0, etf = 0;
        CVodeGetNumSteps(m_mem, &steps);
        CVodeGetNumRhsEvals(m_mem, &rhs);
        CVodeGetNumJacEvals(m_mem, &jac);
        CVodeGetNumLinRhsEvals(m_mem, &linRhs);
        CVodeGetNonlinSolvStats(m_mem, &nli, &ncf);
        CVodeGetNumErrTestFails(m_mem, &etf);
        Cantera::AnyMap s;
        s["steps"] = steps;
        s["rhs_evals"] = rhs;
        s["jac_evals"] = jac;
        s["lin_rhs_evals"] = linRhs;
        s["nonlinear_iters"] = nli;
        s["nonlinear_conv_fails"] = ncf;
        s["err_test_fails"] = etf;
        s["total_steps"] = m_tot[0] + steps;
        s["total_rhs_evals"] = m_tot[1] + rhs;
        s["total_jac_evals"] = m_tot[2] + jac;
        s["total_lin_rhs_evals"] = m_tot[3] + linRhs;
        s["total_nonlinear_iters"] = m_tot[4] + nli;
        s["total_nonlinear_conv_fails"] = m_tot[5] + ncf;
        s["total_err_test_fails"] = m_tot[6] + etf;
        s["restarts"] = m_restarts;
        s["jac_checks"] = m_nChecks;
        s["jac_check_max_error"] = m_checkMax;
        s["jac_time"] = m_jacTime; // [s] wall time spent in ChemPlasReactor::jacobian
        return s;
    }

    // Analytic Jacobian in use (mode Analytic and model supported by ChemPlasReactor::jacobian)
    bool analytic() const {
        return m_mode == Mode::Analytic && m_reactor && m_reactor->analyticJacobianSupported();
    }

private:
    static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* data) {
        auto* self = static_cast<CVodeJac*>(data);
        return self->m_func->evalNoThrow(t, NV_DATA_S(y), NV_DATA_S(ydot));
    }

    static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* data, N_Vector, N_Vector, N_Vector) {
        auto* self = static_cast<CVodeJac*>(data);
        try {
            self->m_nJac++;
            const auto t0 = std::chrono::steady_clock::now();
            const int flag = self->m_reactor->jacobian(t, NV_DATA_S(y), NV_DATA_S(fy), SM_DATA_D(J),
                                                       static_cast<size_t>(SM_LDATA_D(J) / SM_COLUMNS_D(J)));
            self->m_jacTime += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (self->m_checkEvery > 0 && self->m_nJac % self->m_checkEvery == 0) self->checkJacobian(t, y, fy, J);
            return flag;
        } catch (std::exception& err) {
            Cantera::writelog("CVodeJac: Jacobian failed: {}\n", err.what());
            return 1; // recoverable: CVODE retries with a smaller step
        }
    }

    static void errHandler(int, const char*, const char*, char* msg, void* data) {
        static_cast<CVodeJac*>(data)->m_error = msg;
    }

    void applyOptions() {
        // as CVodesIntegrator::applyOptions for "DENSE" (SUNDIALS < 6 with LAPACK)
        if (m_ls) SUNLinSolFree(m_ls);
        if (m_A) SUNMatDestroy(m_A);
        const sunindextype N = static_cast<sunindextype>(m_neq);
        m_A = SUNDenseMatrix(N, N);
        m_ls = SUNLapackDense(m_y, m_A);
        if (!m_A || !m_ls) throw Cantera::CanteraError("CVodeJac::applyOptions", "dense linear solver not created");
        check(CVDlsSetLinearSolver(m_mem, m_ls, m_A), "CVDlsSetLinearSolver");
        if (analytic()) check(CVodeSetJacFn(m_mem, jac), "CVodeSetJacFn"); // after the solver: it resets the Jacobian
        CVodeSetMaxNumSteps(m_mem, m_maxsteps);
    }

    void cvodeStep(double tout, const char* where) {
        const int flag = CVode(m_mem, tout, m_y, &m_tInteg, CV_ONE_STEP);
        if (flag != CV_SUCCESS) {
            std::string errs = m_func->getErrors();
            throw Cantera::CanteraError(std::string("CVodeJac::") + where, "CVodes error encountered. Error code: {}\n{}\n{}",
                                        flag, m_error, errs.empty() ? "" : "Exceptions caught during RHS evaluation:\n" + errs);
        }
    }

    void check(int flag, const char* fn) const {
        if (flag != CV_SUCCESS) throw Cantera::CanteraError("CVodeJac", "{} failed with flag {}: {}", fn, flag, m_error);
    }

    void accumulateStats() {
        long v[7] = {0, 0, 0, 0, 0, 0, 0};
        CVodeGetNumSteps(m_mem, &v[0]);
        CVodeGetNumRhsEvals(m_mem, &v[1]);
        CVodeGetNumJacEvals(m_mem, &v[2]);
        CVodeGetNumLinRhsEvals(m_mem, &v[3]);
        CVodeGetNonlinSolvStats(m_mem, &v[4], &v[5]);
        CVodeGetNumErrTestFails(m_mem, &v[6]);
        for (int i = 0; i < 7; i++) m_tot[i] += v[i];
        m_restarts++;
    }

    // Diagnostics: finite-difference Jacobian at (t, y), as CVODE builds it (increment
    // max(sqrt(eps)|y_j|, sqrt(eps)/w_j) with the error weights w), and the largest difference
    // |J_analytic - J_fd|_ij * (1/w_i) / (1/w_j) relative to the largest scaled entry of J_fd.
    void checkJacobian(double t, N_Vector y, N_Vector fy, SUNMatrix J) {
        const size_t N = m_neq;
        N_Vector w = N_VClone(y);
        CVodeGetErrWeights(m_mem, w);
        const double* yv = NV_DATA_S(y);
        const double* f0 = NV_DATA_S(fy);
        const double* wv = NV_DATA_S(w);
        const double* Ja = SM_DATA_D(J);
        const size_t ld = static_cast<size_t>(SM_LDATA_D(J) / SM_COLUMNS_D(J));
        std::vector<double> yp(yv, yv + N), fp(N), Jfd(N * N);
        const double srur = std::sqrt(std::numeric_limits<double>::epsilon());
        for (size_t j = 0; j < N; j++) {
            const double inc = std::max(srur * std::abs(yv[j]), srur / wv[j]);
            yp[j] = yv[j] + inc;
            m_func->eval(t, yp.data(), fp.data(), nullptr);
            yp[j] = yv[j];
            for (size_t i = 0; i < N; i++) Jfd[i + j * N] = (fp[i] - f0[i]) / inc;
        }
        double ref = 0.0, err = 0.0;
        size_t ie = 0, je = 0;
        for (size_t j = 0; j < N; j++) {
            for (size_t i = 0; i < N; i++) {
                const double s = wv[j] / wv[i]; // scaled entry: J_ij * (1/w_i) / (1/w_j)
                ref = std::max(ref, std::abs(Jfd[i + j * N]) * s);
                const double d = std::abs(Ja[i + j * ld] - Jfd[i + j * N]) * s;
                if (d > err) {
                    err = d;
                    ie = i;
                    je = j;
                }
            }
        }
        N_VDestroy(w);
        const double rel = err / std::max(ref, 1e-300);
        m_nChecks++;
        m_checkMax = std::max(m_checkMax, rel);
        std::printf("\n[Jacobian check %ld] t = %.12g: max scaled |J - J_fd| / max scaled |J_fd| = %.3e at (row %zu, col %zu):"
                    " y_col = %.3e, J = %.4e, J_fd = %.4e\n",
                    m_nChecks, t, rel, ie, je, yv[je], Ja[ie + je * ld], Jfd[ie + je * N]);
    }

    Mode m_mode;
    int m_checkEvery;
    void* m_mem = nullptr;
    SUNMatrix m_A = nullptr;
    SUNLinearSolver m_ls = nullptr;
    N_Vector m_y = nullptr;
    N_Vector m_dky = nullptr;
    Cantera::FuncEval* m_func = nullptr;
    Cantera::ChemPlasReactor* m_reactor = nullptr;
    size_t m_neq = 0;
    double m_reltol = 1e-9, m_abstol = 1e-15;
    double m_time = 0.0, m_tInteg = 0.0;
    int m_maxsteps = 20000;
    std::string m_error;
    long m_tot[7] = {0, 0, 0, 0, 0, 0, 0};
    long m_restarts = 0, m_nJac = 0, m_nChecks = 0;
    double m_checkMax = 0.0;
    double m_jacTime = 0.0;
};

} // namespace CRN

#endif
