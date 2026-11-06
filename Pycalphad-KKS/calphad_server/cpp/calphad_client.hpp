// calphad_client.hpp — C++ client for calphad_server.py
//
// Header-only. Needs two header-only libraries (CMakeLists.txt downloads them):
//   httplib.h  https://github.com/yhirose/cpp-httplib   (MIT)
//   json.hpp   https://github.com/nlohmann/json          (MIT)
//
//   CalphadClient c("127.0.0.1", 8000, "cuni_cost507");
//   auto r = c.equilibrium({{"NI", 0.33}}, 1500.0, "CU");
//   for (auto& ph : r["results"]["phases"])
//       std::cout << ph["phase"] << " " << ph["fraction"] << "\n";
#pragma once

#include <map>
#include <stdexcept>
#include <string>

#include "httplib.h"
#include "json.hpp"

namespace calphad {

using json = nlohmann::json;

class CalphadClient {
public:
    CalphadClient(const std::string& host = "127.0.0.1", int port = 8000,
                  const std::string& db = "")
        : cli_(host, port), db_(db) {
        cli_.set_read_timeout(600, 0);  // big energy surfaces can take a while
    }

    // ---- generic -----------------------------------------------------------
    json get(const std::string& path) {
        auto res = cli_.Get(path);
        return parse(res, path);
    }
    json post(const std::string& path, json body) {
        if (!db_.empty() && !body.contains("db")) body["db"] = db_;
        auto res = cli_.Post(path, body.dump(), "application/json");
        return parse(res, path);
    }

    // ---- databases -----------------------------------------------------------
    json health() { return get("/health"); }
    json databases() { return get("/databases"); }
    json info(const std::string& db = "") { return get("/databases/" + (db.empty() ? db_ : db)); }
    json load_database_path(const std::string& path, const std::string& name = "") {
        json b = {{"path", path}};
        if (!name.empty()) b["name"] = name;
        return post("/databases", b);
    }

    // ---- API 1: equilibrium ------------------------------------------------------
    // composition: solutes (with `balance`) or all elements; mole fractions by default.
    json equilibrium(const std::map<std::string, double>& composition, double T,
                     const std::string& balance = "", json extra = json::object()) {
        return post("/equilibrium", alloy_body(composition, T, balance, extra));
    }
    // several temperatures in one call (much faster than a loop)
    json equilibrium(const std::map<std::string, double>& composition,
                     const std::vector<double>& temps,
                     const std::string& balance = "", json extra = json::object()) {
        return post("/equilibrium", alloy_body(composition, temps, balance, extra));
    }

    // ---- API 2: Gibbs energy of each phase -----------------------------------------
    json phase_energies(const std::map<std::string, double>& composition, double T,
                        const std::string& balance = "", json extra = json::object()) {
        return post("/phase_energies", alloy_body(composition, T, balance, extra));
    }

    // ---- API 3: Gibbs energy surfaces ------------------------------------------------
    // extra may hold "elements", "phases", "mode", "step", "pdens", "max_points", "P"
    json energy_surface(double T, json extra = json::object()) {
        extra["T"] = T;
        return post("/energy_surface", extra);
    }

    // ---- KKS inputs (binary, two phases) ------------------------------------------
    // body: {"solute","phase_a","phase_b","T":[...], "mode","c_alloy","Vm",...}
    json kks_parameters(json body) { return post("/kks/parameters", body); }

private:
    httplib::Client cli_;
    std::string db_;

    template <typename TT>
    json alloy_body(const std::map<std::string, double>& composition, const TT& T,
                    const std::string& balance, json extra) {
        extra["composition"] = composition;
        extra["T"] = T;
        if (!balance.empty()) extra["balance"] = balance;
        return extra;
    }

    static json parse(const httplib::Result& res, const std::string& path) {
        if (!res)
            throw std::runtime_error("calphad server unreachable (" + path + "): " +
                                     httplib::to_string(res.error()));
        json j = json::parse(res->body, nullptr, /*allow_exceptions=*/false);
        if (res->status != 200) {
            std::string msg = j.is_object() && j.contains("detail") ? j["detail"].dump()
                                                                   : res->body;
            throw std::runtime_error("HTTP " + std::to_string(res->status) + " " + path +
                                     ": " + msg);
        }
        return j;
    }
};

}  // namespace calphad
