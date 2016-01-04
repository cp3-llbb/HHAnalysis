#ifndef HHANALYZER_H
#define HHANALYZER_H

#include <cp3_llbb/Framework/interface/Analyzer.h>
#include <cp3_llbb/Framework/interface/Category.h>
#include <cp3_llbb/HHAnalysis/interface/Types.h>

#include <Math/VectorUtil.h>

using namespace HHAnalysis;

class HHAnalyzer: public Framework::Analyzer {
    public:
        HHAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config)
        {
            m_muonLooseIsoCut = config.getUntrackedParameter<double>("muonLooseIsoCut", 0.25);
            m_muonTightIsoCut = config.getUntrackedParameter<double>("muonTightIsoCut", 0.15);
            m_muonEtaCut = config.getUntrackedParameter<double>("muonEtaCut", 2.4);
            m_leadingMuonPtCut = config.getUntrackedParameter<double>("leadingMuonPtCut", 20);
            m_subleadingMuonPtCut = config.getUntrackedParameter<double>("subleadingMuonPtCut", 10);

            m_electronIsoCut_EB_Loose = config.getUntrackedParameter<double>("electronIsoCut_EB_Loose", 0.0893);
            m_electronIsoCut_EE_Loose = config.getUntrackedParameter<double>("electronIsoCut_EE_Loose", 0.121);
            m_electronIsoCut_EB_Tight = config.getUntrackedParameter<double>("electronIsoCut_EB_Tight", 0.0354);
            m_electronIsoCut_EE_Tight = config.getUntrackedParameter<double>("electronIsoCut_EE_Tight", 0.0646);
            m_electronEtaCut = config.getUntrackedParameter<double>("electronEtaCut", 2.5);
            m_leadingElectronPtCut = config.getUntrackedParameter<double>("leadingElectronPtCut", 20);
            m_subleadingElectronPtCut = config.getUntrackedParameter<double>("subleadingElectronPtCut", 15);
            m_electron_loose_wp_name = config.getUntrackedParameter<std::string>("electrons_loose_wp_name", "cutBasedElectronID-Spring15-50ns-V1-standalone-loose");
            m_electron_tight_wp_name = config.getUntrackedParameter<std::string>("electrons_tight_wp_name", "cutBasedElectronID-Spring15-50ns-V1-standalone-tight");

            m_jetEtaCut = config.getUntrackedParameter<double>("jetEtaCut", 2.4);
            m_jetPtCut = config.getUntrackedParameter<double>("jetPtCut", 20);
            m_jet_bDiscrName = config.getUntrackedParameter<std::string>("discr_name", "pfCombinedInclusiveSecondaryVertexV2BJetTags");
            m_jet_bDiscrCut_loose = config.getUntrackedParameter<double>("discr_cut_loose", 0.605);
            m_jet_bDiscrCut_medium = config.getUntrackedParameter<double>("discr_cut_medium", 0.89);
            m_jet_bDiscrCut_tight = config.getUntrackedParameter<double>("discr_cut_tight", 0.97);
            m_minDR_l_j_Cut = config.getUntrackedParameter<double>("minDR_l_j_Cut", 0.3);

            m_hltDRCut = config.getUntrackedParameter<double>("hltDRCut", std::numeric_limits<float>::max());
            m_hltDPtCut = config.getUntrackedParameter<double>("hltDPtCut", std::numeric_limits<float>::max());
        }

        // leptons and dileptons stuff
        BRANCH(electrons, std::vector<unsigned int>);
        BRANCH(muons, std::vector<unsigned int>);
        BRANCH(leptons, std::vector<HH::Lepton>);
        BRANCH(ll, std::vector<HH::Dilepton>);
        BRANCH(met, std::vector<HH::Met>);
        BRANCH(llmet, std::vector<HH::DileptonMet>);
        BRANCH(jets, std::vector<HH::Jet>);
        BRANCH(jj, std::vector<HH::Dijet>);
        BRANCH(llmetjj, std::vector<HH::DileptonMetDijet>);
        // some few custom candidates, for convenience
        // allTight= tight lepton ID, tight lepton Iso, Tight jet ID
        BRANCH(llmetjj_allTight_btagM_ht, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_allTight_btagM_pt, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_allTight_btagM_csv, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_allTight_nobtag_ht, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_allTight_nobtag_pt, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_allTight_nobtag_csv, std::vector<HH::DileptonMetDijet>);

        // maps
        std::vector<std::vector<int>>& map_l_id_iso = tree["map_l_id_iso"].write_with_init<std::vector<std::vector<int>>>(lepID::Count * lepIso::Count, std::vector<int>(0));
        std::vector<std::vector<int>>& map_ll_id_iso = tree["map_ll_id_iso"].write_with_init<std::vector<std::vector<int>>>(lepID::Count * lepIso::Count * lepID::Count * lepIso::Count, std::vector<int>(0));
        // FIXME: add enum over met?
        std::vector<std::vector<int>>& map_llmet_id_iso = tree["map_llmet_id_iso"].write_with_init<std::vector<std::vector<int>>>(lepID::Count * lepIso::Count * lepID::Count * lepIso::Count, std::vector<int>(0));
        std::vector<std::vector<int>>& map_j_btagWP = tree["map_j_btagWP"].write_with_init<std::vector<std::vector<int>>>(jetID::Count * btagWP::Count, std::vector<int>(0));
        std::vector<std::vector<int>>& map_jj_btagWP_pair = tree["map_jj_btagWP_pair"].write_with_init<std::vector<std::vector<int>>>(jetID::Count * jetID::Count * btagWP::Count * btagWP::Count * jetPair::Count, std::vector<int>(0));
        std::vector<std::vector<int>>& map_llmetjj_id_iso_btagWP_pair = tree["map_llmetjj_id_iso_btagWP_pair"].write_with_init<std::vector<std::vector<int>>>(lepID::Count * lepIso::Count * lepID::Count * lepIso::Count * jetID::Count * jetID::Count * btagWP::Count * btagWP::Count * jetPair::Count, std::vector<int>(0));

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const AnalyzersManager&, const CategoryManager&) override;
        virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet& config) override;

        // global event stuff (selected objects multiplicity)
        BRANCH(nJets, unsigned int);
        BRANCH(nJetsL, unsigned int);
        BRANCH(nBJetsL, unsigned int);
        BRANCH(nBJetsM, unsigned int);
        BRANCH(nBJetsT, unsigned int);
        BRANCH(nMuons, unsigned int);
        BRANCH(nMuonsL, unsigned int);
        BRANCH(nMuonsT, unsigned int);
        BRANCH(nElectrons, unsigned int);
        BRANCH(nElectronsL, unsigned int);
        BRANCH(nElectronsT, unsigned int);
        BRANCH(nLeptons, unsigned int);
        BRANCH(nLeptonsL, unsigned int);
        BRANCH(nLeptonsT, unsigned int);

        float m_electronIsoCut_EB_Loose, m_electronIsoCut_EE_Loose, m_electronIsoCut_EB_Tight, m_electronIsoCut_EE_Tight, m_electronEtaCut, m_leadingElectronPtCut, m_subleadingElectronPtCut;
        float m_muonLooseIsoCut, m_muonTightIsoCut, m_muonEtaCut, m_leadingMuonPtCut, m_subleadingMuonPtCut;
        float m_jetEtaCut, m_jetPtCut, m_jet_bDiscrCut_loose, m_jet_bDiscrCut_medium, m_jet_bDiscrCut_tight;
        float m_minDR_l_j_Cut;
        float m_hltDRCut, m_hltDPtCut;
        std::string m_jet_bDiscrName;
        std::string m_electron_loose_wp_name;
        std::string m_electron_tight_wp_name;

        float getCosThetaStar_CS(const LorentzVector & h1, const LorentzVector & h2, float ebeam = 6500)
        {// cos theta star angle in the Collins Soper frame
            LorentzVector p1, p2;
            p1.SetPxPyPzE(0, 0,  ebeam, ebeam);
            p2.SetPxPyPzE(0, 0, -ebeam, ebeam);

            LorentzVector hh = h1 + h2;
            ROOT::Math::Boost boost(-hh.X() / hh.T(), -hh.Y() / hh.T(), -hh.Z() / hh.T());
            p1 = boost(p1);
            p2 = boost(p2);
            LorentzVector newh1 = boost(h1);
            ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> CSaxis(p1.Vect().Unit() - p2.Vect().Unit());

            return cos(ROOT::Math::VectorUtil::Angle(CSaxis.Unit(), newh1.Vect().Unit()));
        }
};

#endif
