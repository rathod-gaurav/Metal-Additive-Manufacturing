// kks_example.cpp — fetch KKS parameters from the server once, then run the
// KKS split and driving force locally (as a phase-field solver would).
//   ./kks_example [host] [port]
// Server:  python calphad_server.py --tdb alni_dupin_2001.tdb
#include <cstdio>
#include <iostream>

#include "calphad_client.hpp"
#include "kks_parabolic.hpp"

using calphad::json;

int main(int argc, char** argv) {
    std::string host = argc > 1 ? argv[1] : "127.0.0.1";
    int port = argc > 2 ? std::stoi(argv[2]) : 8000;
    try {
        calphad::CalphadClient c(host, port, "alni_dupin_2001");

        // 1. once, at start-up: gamma (FCC_A1, phi=1) / liquid, Ni-10Al, 1700-1720 K
        json body = {{"solute", "AL"}, {"phase_a", "FCC_A1"}, {"phase_b", "LIQUID"},
                     {"mode", "parabolic"}, {"c_alloy", 0.10}, {"Vm", 7.0e-6}};
        std::vector<double> temps;
        for (double T = 1700; T <= 1720; T += 2.5) temps.push_back(T);
        body["T"] = temps;
        json p = c.kks_parameters(body);

        calphad::KKSParabolic kks;
        kks.Vm = 7.0e-6;
        for (auto& r : p["tie_lines"]) {
            if (r.is_null()) continue;
            kks.add(r["T"], r["c_a"], r["c_b"], r["mu"], r["G_a"], r["G_b"], r["A_a"], r["A_b"]);
        }
        std::printf("loaded %zu tie lines, %.0f-%.0f K\n", kks.T.size(), kks.T.front(), kks.T.back());
        auto k = kks.at(1710.0);
        std::printf("T=1710 K: c_gamma=%.5f c_L=%.5f k=%.4f A_gamma=%.0f A_L=%.0f J/mol\n",
                    k.ca, k.cb, k.ca / k.cb, k.Aa, k.Ab);

        // 2. in the time loop: KKS split + driving force at a few points
        double cs[] = {0.10, 0.12, 0.14};
        double phis[] = {0.0, 0.5, 1.0};
        for (double cc : cs)
            for (double ph : phis) {
                double ca, cb, mu;
                kks.split(cc, ph, 1710.0, ca, cb, mu);
                double dF = kks.driving_force(ca, cb, mu, 1710.0);
                std::printf("  c=%.2f phi=%.1f -> c_a=%.5f c_b=%.5f mu=%.1f J/mol dF=%.4e J/m^3\n",
                            cc, ph, ca, cb, mu, dF);
            }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
