#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cp3_llbb/HHAnalysis/interface/Categories.h>

// NB : All di-lepton categories require also two jets, as we do not want to store events without at least two leptons and two jets.
//      The criteria "having two b-jets" is passed as a cut in the category in question.

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

const unsigned int DileptonCategory::getNJets(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hh_analyzer.nJets;
}

const unsigned int DileptonCategory::getNBJets(const AnalyzersManager& analyzers) const {
    const HHAnalyzer& hh_analyzer = analyzers.get<HHAnalyzer>("hh_analyzer");
    return hh_analyzer.nBJetsM;
}

// ***** ***** *****
// Dilepton Mu-Mu category
// ***** ***** *****
bool MuMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const JetsProducer& alljets = producers.get<JetsProducer>("jets");
    return (muons.p4.size() >= 2 && alljets.p4.size() >= 2);
};

bool MuMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HH::Dilepton>& ll = getDileptons(analyzers);
    bool isMuMu = false;
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isMuMu) isMuMu = true;
    }
    const unsigned int nJets = getNJets(analyzers);
    return (isMuMu && nJets >= 2);
};

void MuMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
    manager.new_cut("ll_mass_lowerZcut", "mll > 85");
    manager.new_cut("has_two_bJets", "nBJetM >= 2");
    manager.new_cut("fire_trigger_Mu17_Mu8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
    manager.new_cut("fire_trigger_Mu17_TkMu8", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
    manager.new_cut("fire_trigger_IsoMu27", "HLT_IsoMu27_v*");
};

void MuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HH::Dilepton>& ll = getDileptons(analyzers);
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isMuMu) {
            if (ll[idilep].p4.M() > m_mll_cut) manager.pass_cut("ll_mass");
            if (ll[idilep].p4.M() > m_mll_lowerZcut) manager.pass_cut("ll_mass_lowerZcut");
        }
    }
    const unsigned int nBJets = getNBJets(analyzers);
    if (nBJets >=2) manager.pass_cut("has_two_bJets"); 
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") != std::string::npos) manager.pass_cut("fire_trigger_Mu17_Mu8");
        if (path.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") != std::string::npos) manager.pass_cut("fire_trigger_Mu17_TkMu8");
    }
}

// ***** ***** *****
// Dilepton El-El category
// ***** ***** *****
bool ElElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const JetsProducer& alljets = producers.get<JetsProducer>("jets");
    return (electrons.p4.size() >= 2  && alljets.p4.size() >= 2);
};

bool ElElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HH::Dilepton>& ll = getDileptons(analyzers);
    bool isElEl = false;
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isElEl) isElEl = true;
    }
    const unsigned int nJets = getNJets(analyzers);
    return (isElEl && nJets >= 2);
};

void ElElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
    manager.new_cut("ll_mass_lowerZcut", "mll > 85");
    manager.new_cut("has_two_bJets", "nBJet >= 2");
    manager.new_cut("fire_trigger_Ele17_Ele12", "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ*");
    manager.new_cut("fire_trigger_Ele32_eta2p1_WPTight_Gsf", "HLT_Ele32_eta2p1_WPTight_Gsf_v*");
};

void ElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HH::Dilepton>& ll = getDileptons(analyzers);
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isElEl) {
            if (ll[idilep].p4.M() > m_mll_cut) manager.pass_cut("ll_mass");
            if (ll[idilep].p4.M() > m_mll_lowerZcut) manager.pass_cut("ll_mass_lowerZcut");
        }
    }
    const unsigned int nBJets = getNBJets(analyzers);
    if (nBJets >=2) manager.pass_cut("has_two_bJets"); 
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos) manager.pass_cut("fire_trigger_Ele17_Ele12");
    }
}

// ***** ***** *****
// Dilepton El-Mu category
// ***** ***** *****
bool ElMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
    const JetsProducer& alljets = producers.get<JetsProducer>("jets");
    return ((electrons.p4.size() + muons.p4.size() >= 2)  && alljets.p4.size() >= 2);
};

bool ElMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HH::Dilepton>& ll = getDileptons(analyzers);
    bool isElMu = false;
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isElMu) isElMu = true;
    }
    const unsigned int nJets = getNJets(analyzers);
    return (isElMu && nJets >= 2);
};

void ElMuCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
    manager.new_cut("ll_mass_lowerZcut", "mll > 85");
    manager.new_cut("has_two_bJets", "nBJet >= 2");
    manager.new_cut("fire_trigger_Mu8_Ele17", "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_*");
};

void ElMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HH::Dilepton>& ll = getDileptons(analyzers);
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isElMu) {
            if (ll[idilep].p4.M() > m_mll_cut) manager.pass_cut("ll_mass");
            if (ll[idilep].p4.M() > m_mll_lowerZcut) manager.pass_cut("ll_mass_lowerZcut");
        }
    }
    const unsigned int nBJets = getNBJets(analyzers);
    if (nBJets >=2) manager.pass_cut("has_two_bJets"); 
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
    const JetsProducer& alljets = producers.get<JetsProducer>("jets");
    return ((electrons.p4.size() + muons.p4.size() >= 2)  && alljets.p4.size() >= 2);
};

bool MuElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HH::Dilepton>& ll = getDileptons(analyzers);
    bool isMuEl = false;
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isMuEl) isMuEl = true;
    }
    const unsigned int nJets = getNJets(analyzers);
    return (isMuEl && nJets >= 2);
};

void MuElCategory::register_cuts(CutManager& manager) {
    manager.new_cut("ll_mass", "mll > 20");
    manager.new_cut("ll_mass_lowerZcut", "mll > 85");
    manager.new_cut("has_two_bJets", "nBJet >= 2");
    manager.new_cut("fire_trigger_Mu17_Ele12", "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_*");
};

void MuElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
    const std::vector<HH::Dilepton>& ll = getDileptons(analyzers);
    for (unsigned int idilep = 0; idilep < ll.size(); idilep++) 
    {
        if (ll[idilep].isMuEl) {
            if (ll[idilep].p4.M() > m_mll_cut) manager.pass_cut("ll_mass");
            if (ll[idilep].p4.M() > m_mll_lowerZcut) manager.pass_cut("ll_mass_lowerZcut");
        }
    }
    const unsigned int nBJets = getNBJets(analyzers);
    if (nBJets >=2) manager.pass_cut("has_two_bJets"); 
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    for (const std::string& path: hlt.paths) 
    {
        if (path.find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos) manager.pass_cut("fire_trigger_Mu17_Ele12");
    }
}

