// CRN_v7: steady state of the reactor network.
//
// The state of every node (reactor or reservoir) is sampled at the checks: flow updates on a pulse
// boundary, one every w pulse periods (clustering, steadyStateCluster). For reactors T and Y_k are
// averaged over the block of w periods since the last check (the state passed downstream; the
// instantaneous short-lived species at the end of a period are at the end of their decay and jump
// between pulses without any drift); p, mdot_in and mdot_out are taken at the check. The derivative
// between two blocks, w periods apart, damps the pulse-to-pulse noise by ~w. Reservoirs are at steady
// state by definition and are not checked. For a quantity x of the node (T, p, mdot_in, Y_k) the change between
// two checks gives the relative rate projected over the time scale tau of the node (steadyStateTau, or
// the residence time m/mdot_out; without outflow the filling time m/mdot_in; for a closed reactor the
// window of W pulse periods):
//
//     d = |x1 - x0| / max(|x0|, |x1|) * tau / (t1 - t0)
//
// The mass balance enters the same way: with a net flow the reactor fills or empties, and the relative
// change of its mass m over tau is d_m = |mdot_in - mdot_out| * tau / m (= 1 for a reactor that only
// fills or only empties: never at steady state).
//
// d is the relative change of x over tau if the trend of the last interval goes on, i.e. an estimate of
// the relative distance from the steady state for an approach with time constant ~tau. One criterion
// covers both cases: small quantities are measured relative to themselves (mass fractions with
// max(|x0|, |x1|) <= steadyStateYmin are skipped, round-off level), and the derivative is always taken
// over a whole check interval (>= one pulse period), never over a small time step.
//
// Nodes are grouped in strongly connected components of the flow graph (a group = one reactor, or the
// reactors of a recirculation loop). At every check the counter n of a group is increased when d < tol
// for every quantity of every node of the group and every node upstream of the group has converged;
// otherwise it is reset (a single node of the group above tol resets it). A group converges when
// n * w >= W, i.e. after W pulse periods below tol: its reactors are no longer integrated (the outflow
// keeps the last averaged state) and its downstream groups go on. A converged group is released (and
// integrated again) by release() when the inflows of one of its reactors change.

#ifndef CRN_STEADY_STATE_H
#define CRN_STEADY_STATE_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

namespace SteadyState
{

struct Settings {
    int W = 0;            // pulse periods below tol to converge (0 = steady state monitor off)
    int w = 1;            // pulse periods per check (clustering)
    double tol = 1e-3;    // threshold on d
    double Ymin = 1e-12;  // mass fractions below this in both samples are not checked
    double tau = 0.0;     // time scale of d [s]; 0 = residence (or filling) time of each reactor
    bool stop = true;     // stop the simulation when every group has converged
};

struct Node {
    std::string name;
    bool reservoir = false;
    size_t index = 0;                // index in reactors or reservoirs
    std::vector<size_t> upstream;    // nodes with a mass flow into this node
    int group = -1;
    // last sample: [T, p, mdot_in, mdot_out, Y_0 .. Y_K-1]
    std::vector<double> x0;
    double t0 = -1.0;
    double d = 0.0;                  // largest d of the last check
    std::string worst;               // quantity with the largest d
    bool converged = false;
};

class Monitor {
public:
    Settings set;
    std::vector<Node> nodes;
    std::vector<std::vector<size_t>> groups; // node indices of every group, upstream groups first
    std::vector<int> counter;                // per group

    // Groups from the flow graph (Tarjan): call after the nodes and their upstream lists are set
    void buildGroups() {
        const size_t n = nodes.size();
        std::vector<std::vector<size_t>> down(n);
        for (size_t v = 0; v < n; v++)
            for (size_t u : nodes[v].upstream) down[u].push_back(v);
        std::vector<int> order(n, -1), low(n, 0);
        std::vector<bool> onStack(n, false);
        std::vector<size_t> stack;
        int next = 0;
        std::vector<std::vector<size_t>> sccs; // Tarjan gives the downstream components first
        std::function<void(size_t)> visit = [&](size_t v) {
            order[v] = low[v] = next++;
            stack.push_back(v);
            onStack[v] = true;
            for (size_t w : down[v]) {
                if (order[w] < 0) {
                    visit(w);
                    low[v] = std::min(low[v], low[w]);
                } else if (onStack[w]) {
                    low[v] = std::min(low[v], order[w]);
                }
            }
            if (low[v] == order[v]) {
                std::vector<size_t> comp;
                size_t w;
                do {
                    w = stack.back();
                    stack.pop_back();
                    onStack[w] = false;
                    comp.push_back(w);
                } while (w != v);
                sccs.push_back(comp);
            }
        };
        for (size_t v = 0; v < n; v++)
            if (order[v] < 0) visit(v);
        groups.assign(sccs.rbegin(), sccs.rend());
        counter.assign(groups.size(), 0);
        for (size_t g = 0; g < groups.size(); g++)
            for (size_t v : groups[g]) nodes[v].group = static_cast<int>(g);
        for (auto& nd : nodes)
            if (nd.reservoir) nd.converged = true; // at steady state by definition
    }

    // d of node v from its new sample x1 at time t1 (tau [s] > 0); the sample is stored for the next check.
    // netRate = (mdot_in - mdot_out) / m [1/s] (mass balance). names: names of the components of x (for
    // the worst quantity). Returns false at the first sample.
    bool evaluate(size_t v, const std::vector<double>& x1, double t1, double tau, double netRate,
                  const std::vector<std::string>& names) {
        Node& nd = nodes[v];
        const bool first = nd.t0 < 0.0;
        nd.d = 0.0;
        nd.worst.clear();
        if (!first) {
            const double scale = tau / (t1 - nd.t0);
            const double dMass = std::abs(netRate) * tau;
            if (dMass > 0.0) {
                nd.d = dMass;
                nd.worst = "mass";
            }
            for (size_t i = 0; i < x1.size(); i++) {
                const double ref = std::max(std::abs(nd.x0[i]), std::abs(x1[i]));
                const double xmin = (i >= 4) ? set.Ymin : 0.0; // only mass fractions have a floor
                if (ref <= xmin || ref == 0.0) continue;
                const double di = std::abs(x1[i] - nd.x0[i]) / ref * scale;
                if (di > nd.d) {
                    nd.d = di;
                    nd.worst = names[i];
                }
            }
        }
        nd.x0 = x1;
        nd.t0 = t1;
        return !first;
    }

    // Counters and convergence of the groups after evaluate() of every node not converged yet.
    // Returns the groups that converged at this check.
    std::vector<size_t> update(const std::vector<bool>& evaluated) {
        std::vector<size_t> newly;
        for (size_t g = 0; g < groups.size(); g++) {
            if (groupConverged(g)) continue;
            bool ok = true;
            for (size_t v : groups[g]) {
                if (nodes[v].reservoir) continue;
                if (!evaluated[v] || nodes[v].d >= set.tol) ok = false; // a single reactor resets the group
                for (size_t u : nodes[v].upstream)
                    if (nodes[u].group != static_cast<int>(g) && !nodes[u].converged) ok = false;
            }
            counter[g] = ok ? counter[g] + 1 : 0;
            if (counter[g] * set.w >= set.W) {
                for (size_t v : groups[g]) nodes[v].converged = true;
                newly.push_back(g);
            }
        }
        return newly;
    }

    // Reactors of group g converged (reservoirs are at steady state by definition)
    bool groupConverged(size_t g) const {
        for (size_t v : groups[g])
            if (!nodes[v].converged) return false;
        return true;
    }

    // The group of node v is no longer converged (new inflows): counter and samples restart.
    // Returns the reactor nodes of the group.
    std::vector<size_t> release(size_t v) {
        const size_t g = static_cast<size_t>(nodes[v].group);
        std::vector<size_t> members;
        counter[g] = 0;
        for (size_t u : groups[g]) {
            if (nodes[u].reservoir) continue;
            nodes[u].converged = false;
            nodes[u].t0 = -1.0;
            members.push_back(u);
        }
        return members;
    }

    bool allConverged() const {
        for (const auto& nd : nodes)
            if (!nd.converged) return false;
        return true;
    }
};

} // namespace SteadyState

#endif
