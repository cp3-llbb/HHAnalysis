#include <cp3_llbb/HHAnalysis/interface/WeightedEfficiency.h>

WeightedEfficiency::WeightedEfficiency(const std::vector<edm::ParameterSet>& parts):
    random_generator(42) {

    std::vector<double> weights;

    for (const auto& p: parts) {
        double weight = p.getUntrackedParameter<double>("weight");
        weights.push_back(weight);

        BinnedValuesJSONParser parser(p.getUntrackedParameter<edm::FileInPath>("file").fullPath());
        efficiencies.push_back(std::move(parser.get_values()));
    }

    probability_distribution.reset(new std::discrete_distribution<>(weights.begin(), weights.end()));
}

std::vector<float> WeightedEfficiency::get(const Parameters& parameters) {
    return efficiencies[(*probability_distribution)(random_generator)].get(parameters);
}
