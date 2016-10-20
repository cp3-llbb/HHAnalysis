#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cp3_llbb/HHAnalysis/interface/Categories.h>

// ***** ***** *****
// Dilepton categories
// ***** ***** *****

const std::vector<HH::Lepton>& DileptonCategory::getLeptons(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>(m_analyzer_name);
    return hh_analyzer.leptons;
}

const std::vector<HH::Dilepton>& DileptonCategory::getDileptons(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>(m_analyzer_name);
    return hh_analyzer.ll;
}

const std::vector<HH::DileptonMetDijet>& DileptonCategory::getDileptonMetDijets(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>(m_analyzer_name);
    return hh_analyzer.llmetjj;
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
    const std::vector<HH::DileptonMetDijet>& llmetjj = getDileptonMetDijets(analyzers);
    bool isMuMu = false;
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isMuMu) isMuMu = true;
    }
    return (isMuMu && llmetjj.size() > 0);
};

void MuMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger_MuMu", "HLT_Mu*");
};

void MuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Mu") != std::string::npos) manager.pass_cut("fire_trigger_MuMu");
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
    const std::vector<HH::DileptonMetDijet>& llmetjj = getDileptonMetDijets(analyzers);
    bool isElEl = false;
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isElEl) isElEl = true;
    }
    return (isElEl && llmetjj.size() > 0);
};

void ElElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger_EleEle", "HLT_Ele*");
};

void ElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Ele") != std::string::npos) manager.pass_cut("fire_trigger_EleEle");
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
    const std::vector<HH::DileptonMetDijet>& llmetjj = getDileptonMetDijets(analyzers);
    bool isElMu = false;
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isElMu) isElMu = true;
    }
    return (isElMu && llmetjj.size() > 0);
};

void ElMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger_MuEle", "HLT_Mu*Ele*");
};

void ElMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        // admittedly this could be either Mu-Ele or Ele-Mu trigger paths
        if (path.find("HLT_Mu") != std::string::npos
            && path.find("Ele") != std::string::npos)
            manager.pass_cut("fire_trigger_MuEle");
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
    const std::vector<HH::DileptonMetDijet>& llmetjj = getDileptonMetDijets(analyzers);
    bool isMuEl = false;
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isMuEl) isMuEl = true;
    }
    return (isMuEl && llmetjj.size() > 0);
};

void MuElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("fire_trigger_MuEle", "HLT_Mu*Ele*");
};

void MuElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Mu") != std::string::npos
            && path.find("Ele") != std::string::npos)
            manager.pass_cut("fire_trigger_MuEle");
    }
}

