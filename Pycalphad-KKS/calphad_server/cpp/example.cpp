// example.cpp — call the CALPHAD server from C++.
//   ./calphad_example [host] [port]
// Start the server first:  python calphad_server.py --tdb cuni_cost507.tdb
#include <cstdio>
#include <iostream>

#include "calphad_client.hpp"

using calphad::json;

int main(int argc, char** argv) {
    std::string host = argc > 1 ? argv[1] : "127.0.0.1";
    int port = argc > 2 ? std::stoi(argv[2]) : 8000;

    try {
        calphad::CalphadClient c(host, port, "cuni_cost507");
        std::cout << "server: " << c.health().dump() << "\n\n";

        // API 1 — Cu-33Ni (at.%) at 1500 K
        json eq = c.equilibrium({{"NI", 0.33}}, 1500.0, "CU");
        std::cout << "API 1  equilibrium at 1500 K\n";
        for (auto& ph : eq["results"]["phases"])
            std::printf("  %-8s fraction %.4f   x_Ni %.4f   G %.1f J/mol\n",
                        ph["phase"].get<std::string>().c_str(), ph["fraction"].get<double>(),
                        ph["composition"]["NI"].get<double>(), ph["GM"].get<double>());

        // Tie lines at several temperatures in one call (KKS inputs vs T)
        json steps = c.equilibrium({{"NI", 0.5}}, std::vector<double>{1550, 1575, 1600}, "CU");
        for (auto& r : steps["results"]) {
            std::printf("  T=%.0f:", r["T"].get<double>());
            for (auto& ph : r["phases"])
                std::printf("  %s x_Ni=%.4f", ph["phase"].get<std::string>().c_str(),
                            ph["composition"]["NI"].get<double>());
            std::printf("\n");
        }

        // API 2 — every phase's G at the alloy composition
        json pe = c.phase_energies({{"NI", 0.33}}, 1500.0, "CU");
        std::cout << "\nAPI 2  phase energies at x_Ni = 0.33\n";
        for (auto& ph : pe["phases"])
            std::printf("  %-8s G %.2f J/mol   mu_Ni - mu_Cu = %.2f J/mol\n",
                        ph["phase"].get<std::string>().c_str(), ph["GM"].get<double>(),
                        ph["chemical_potentials"]["NI"].get<double>() -
                            ph["chemical_potentials"]["CU"].get<double>());

        // API 3 — G(x) curves
        json s = c.energy_surface(1500.0, {{"step", 0.1}});
        std::cout << "\nAPI 3  energy surface, step 0.1\n";
        for (auto& ph : s["phases"]) {
            std::printf("  %-8s", ph["phase"].get<std::string>().c_str());
            for (size_t i = 0; i < ph["GM"].size(); ++i)
                std::printf(" (%.1f, %.0f)", ph["X"][i][1].get<double>(), ph["GM"][i].get<double>());
            std::printf("\n");
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
