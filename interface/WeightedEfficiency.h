#pragma once

#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <cp3_llbb/Framework/interface/BinnedValuesJSONParser.h>

#include <memory>
#include <random>
#include <vector>

class WeightedEfficiency {
    public:
        WeightedEfficiency(const std::vector<edm::ParameterSet>& parts);
        /**
         * Randomly select one set of HLT efficiencies from the ones
         * available, according to the fraction of integrated luminosity used to
         * compute the efficiency. A set of efficiencies evaluated on more luminosity
         * will be sampled more than one evaluated on less
         */
        std::vector<float> get(const Parameters&);

    private:
        std::mt19937 random_generator;
        std::unique_ptr<std::discrete_distribution<int>> probability_distribution;
        std::vector<BinnedValues> efficiencies;
};
