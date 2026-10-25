// crnDict reader shared by ChemPlasKin (masterCRN.cpp) and EEDFTable (eedfTableBuilder.cpp)
#ifndef CRN_DICT_H
#define CRN_DICT_H

#include <fstream>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

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

inline std::string makeUniqueReactorName(const std::string& baseName, const std::map<std::string, ReactorZoneDef>& reactors) {
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
inline std::map<std::string, ReactorZoneDef> readReactors(const std::string& fileName) {
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

#endif // CRN_DICT_H
