#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cp3_llbb/HHAnalysis/interface/Categories.h>

// ***** ***** *****
// Dilepton categories
// ***** ***** *****

const std::vector<HH::Lepton>& DileptonCategory::getLeptons(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hh_analyzer.leptons;
}

const std::vector<HH::Dilepton>& DileptonCategory::getDileptons(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hh_analyzer.ll;
}

// ***** ***** *****
// Dilepton Mu-Mu category
// ***** ***** *****
bool MuMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    return (muons.p4.size() >= 2);
};

bool MuMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HH::Dilepton>& ll = getDileptons(analyzers);
    bool isMuMu = false;
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isMuMu) isMuMu = true;
    }
    return (isMuMu);
};

void MuMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger_Mu17_Mu8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
    manager.new_cut("fire_trigger_Mu17_TkMu8", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
    manager.new_cut("fire_trigger_IsoMu27", "HLT_IsoMu27_v*");
};

void MuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") != std::string::npos) manager.pass_cut("fire_trigger_Mu17_Mu8");
        if (path.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") != std::string::npos) manager.pass_cut("fire_trigger_Mu17_TkMu8");
        if (path.find("HLT_IsoMu27_v") != std::string::npos) manager.pass_cut("fire_trigger_IsoMu27");
    }
}

// ***** ***** *****
// Dilepton El-El category
// ***** ***** *****
bool ElElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    return (electrons.p4.size() >= 2);
};

bool ElElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HH::Dilepton>& ll = getDileptons(analyzers);
    bool isElEl = false;
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isElEl) isElEl = true;
    }
    return (isElEl);
};

void ElElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger_Ele17_Ele12", "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ*");
    //manager.new_cut("fire_trigger_Ele23_WPLoose", "HLT_Ele23_WPLoose_Gsf_v*");
};

void ElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) manager.pass_cut("fire_trigger_Ele17_Ele12");
        //if (path.find("HLT_Ele23_WPLoose_Gsf_v") != std::string::npos) manager.pass_cut("fire_trigger_Ele23_WPLoose");
    }
}

// ***** ***** *****
// Dilepton El-Mu category
// ***** ***** *****
bool ElMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    return ((electrons.p4.size() + muons.p4.size()) >= 2);
};

bool ElMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HH::Dilepton>& ll = getDileptons(analyzers);
    bool isElMu = false;
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isElMu) isElMu = true;
    }
    return (isElMu);
};

void ElMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger_Mu8_Ele17", "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_*");
};

void ElMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) manager.pass_cut("fire_trigger_Mu8_Ele17");
    }
}

// ***** ***** *****
// Dilepton Mu-El category
// ***** ***** *****
bool MuElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    return ((electrons.p4.size() + muons.p4.size()) >= 2);
};

bool MuElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HH::Dilepton>& ll = getDileptons(analyzers);
    bool isMuEl = false;
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isMuEl) isMuEl = true;
    }
    return (isMuEl);
};

void MuElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger_Mu17_Ele12", "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_*");
};

void MuElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) manager.pass_cut("fire_trigger_Mu17_Ele12");
    }
}

