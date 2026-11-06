// kks_parabolic.hpp — KKS with CALPHAD-fitted parabolic free energies, for a
// C++ phase-field solver. Header-only, no dependencies.
//
// Coefficients come from Python (KKSModel.export_parabolic_csv) or from the
// server (POST /kks/parameters, see kks_example.cpp). Per temperature T:
//   f_a(c) = G_a + mu (c - ca) + A_a/2 (c - ca)^2      phase a: phi = 1
//   f_b(c) = G_b + mu (c - cb) + A_b/2 (c - cb)^2      phase b: phi = 0
// (ca, cb, mu, G_a, G_b): the tie line at T; A = d2G/dc2 there. Coefficients
// are interpolated linearly in T, so T can be a field (e.g. a moving
// temperature gradient in additive manufacturing).
//
//   auto kks = calphad::KKSParabolic::from_csv("nial_gamma_liquid.csv", 7.0e-6);
//   double ca, cb, mu;
//   kks.split(c, phi, T, ca, cb, mu);              // KKS constraint
//   double dF = kks.driving_force(ca, cb, mu, T);  // J/m^3
//   // dphi/dt = M [eps2*lap(phi) - W*dg(phi) - dh(phi)*dF]
//   // dc/dt   = div[ D(phi) (h grad(ca) + (1-h) grad(cb)) ]
#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace calphad {

struct KKSParabolic {
    std::vector<double> T, ca, cb, mu, Ga, Gb, Aa, Ab;  // one entry per temperature
    double Vm = 0.0;                                    // m^3/mol

    struct Coef { double ca, cb, mu, Ga, Gb, Aa, Ab; };

    // CSV written by KKSModel.export_parabolic_csv: lines starting with '#' are
    // comments; header T_K,c_a_eq,c_b_eq,mu_eq,G_a_eq,G_b_eq,A_a,A_b
    static KKSParabolic from_csv(const std::string& path, double Vm) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("cannot open " + path);
        KKSParabolic k;
        k.Vm = Vm;
        std::string line;
        bool header = false;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            if (!header) { header = true; continue; }
            std::stringstream ss(line);
            std::string tok;
            std::vector<double> v;
            while (std::getline(ss, tok, ',')) v.push_back(std::stod(tok));
            if (v.size() != 8) throw std::runtime_error("bad CSV row: " + line);
            k.add(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
        }
        if (k.T.empty()) throw std::runtime_error("no rows in " + path);
        return k;
    }

    void add(double T_, double ca_, double cb_, double mu_, double Ga_, double Gb_,
             double Aa_, double Ab_) {
        T.push_back(T_); ca.push_back(ca_); cb.push_back(cb_); mu.push_back(mu_);
        Ga.push_back(Ga_); Gb.push_back(Gb_); Aa.push_back(Aa_); Ab.push_back(Ab_);
    }

    // coefficients at temperature Tq (linear interpolation, clamped to the table)
    Coef at(double Tq) const {
        if (T.size() == 1 || Tq <= T.front()) return row(0);
        if (Tq >= T.back()) return row(T.size() - 1);
        size_t i = std::upper_bound(T.begin(), T.end(), Tq) - T.begin() - 1;
        double w = (Tq - T[i]) / (T[i + 1] - T[i]);
        Coef a = row(i), b = row(i + 1);
        auto L = [w](double p, double q) { return p + w * (q - p); };
        return {L(a.ca, b.ca), L(a.cb, b.cb), L(a.mu, b.mu), L(a.Ga, b.Ga),
                L(a.Gb, b.Gb), L(a.Aa, b.Aa), L(a.Ab, b.Ab)};
    }

    // interpolation and double-well functions
    static double h(double p) { return p * p * p * (6 * p * p - 15 * p + 10); }
    static double dh(double p) { return 30 * p * p * (1 - p) * (1 - p); }
    static double g(double p) { return p * p * (1 - p) * (1 - p); }
    static double dg(double p) { return 2 * p * (1 - p) * (1 - 2 * p); }

    // free energies, J/mol
    static double f_a(double c, const Coef& k) {
        double d = c - k.ca; return k.Ga + k.mu * d + 0.5 * k.Aa * d * d;
    }
    static double f_b(double c, const Coef& k) {
        double d = c - k.cb; return k.Gb + k.mu * d + 0.5 * k.Ab * d * d;
    }

    // KKS split: c = h c_a + (1-h) c_b with equal diffusion potential (closed form)
    void split(double c, double phi, double Tq, double& c_a, double& c_b, double& mu_out) const {
        Coef k = at(Tq);
        double hh = h(phi), r = k.Aa / k.Ab;
        c_a = (c - (1 - hh) * (k.cb - r * k.ca)) / (hh + (1 - hh) * r);
        c_b = k.cb + r * (c_a - k.ca);
        mu_out = k.mu + k.Aa * (c_a - k.ca);  // J/mol
    }

    // f_a(c_a) - f_b(c_b) - mu (c_a - c_b), in J/m^3 (negative favours phase a)
    double driving_force(double c_a, double c_b, double mu_, double Tq) const {
        if (Vm <= 0) throw std::runtime_error("set Vm > 0");
        Coef k = at(Tq);
        return (f_a(c_a, k) - f_b(c_b, k) - mu_ * (c_a - c_b)) / Vm;
    }

    // d(mu)/dc at fixed phi, J/mol (for dc/dt = div(D/dmu_dc grad mu))
    double dmu_dc(double phi, double Tq) const {
        Coef k = at(Tq);
        double hh = h(phi);
        return 1.0 / (hh / k.Aa + (1 - hh) / k.Ab);
    }

    // sigma [J/m^2], interface width 2*lambda [m] -> eps^2 [J/m], W [J/m^3] (KKS 1999)
    static void interface_parameters(double sigma, double width, double& eps2, double& W) {
        eps2 = 3.0 * sigma * width / 2.2;
        W = 13.2 * sigma / width;
    }

private:
    Coef row(size_t i) const { return {ca[i], cb[i], mu[i], Ga[i], Gb[i], Aa[i], Ab[i]}; }
};

}  // namespace calphad
