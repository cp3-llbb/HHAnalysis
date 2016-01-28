#ifndef HHANALYZER_H
#define HHANALYZER_H

#include <cp3_llbb/Framework/interface/Analyzer.h>
#include <cp3_llbb/Framework/interface/Category.h>
#include <cp3_llbb/Framework/interface/ScaleFactorParser.h>
#include <cp3_llbb/HHAnalysis/interface/Types.h>

#include <Math/VectorUtil.h>

using namespace HH;
using namespace HHAnalysis;

class HHAnalyzer: public Framework::Analyzer {
    public:
        HHAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config)
        {
            // Not untracked as these parameters are mandatory
            m_electrons_producer = config.getParameter<std::string>("electronsProducer");
            m_muons_producer = config.getParameter<std::string>("muonsProducer");
            m_jets_producer = config.getParameter<std::string>("jetsProducer");
            m_met_producer = config.getParameter<std::string>("metProducer");
            m_nohf_met_producer = config.getParameter<std::string>("nohfMETProducer");
            // other parameters
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

            if (config.exists("hlt_efficiencies")){
                const edm::ParameterSet& hlt_efficiencies = config.getUntrackedParameter<edm::ParameterSet>("hlt_efficiencies");
                std::vector<std::string> hlt_efficiencies_name = hlt_efficiencies.getParameterNames();
                for (const std::string& hlt_efficiency: hlt_efficiencies_name) {
                    std::cout << "Adding hlt efficiencies : " << hlt_efficiency << std::endl;
                    ScaleFactorParser parser(hlt_efficiencies.getUntrackedParameter<edm::FileInPath>(hlt_efficiency).fullPath());
                    m_hlt_efficiencies.emplace(hlt_efficiency, std::move(parser.get_scale_factor()));
                }
            }
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
        BRANCH(llmetjj_allTight_btagL_ht, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_allTight_btagL_pt, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_allTight_btagL_csv, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_allTight_nobtag_ht, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_allTight_nobtag_pt, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_allTight_nobtag_csv, std::vector<HH::DileptonMetDijet>);

        // maps
        std::vector<std::vector<int>>& map_l = tree["map_l"].write_with_init<std::vector<std::vector<int>>>(lepID::Count * lepIso::Count, std::vector<int>(0));
        std::vector<std::vector<int>>& map_ll = tree["map_ll"].write_with_init<std::vector<std::vector<int>>>(lepID::Count * lepIso::Count * lepID::Count * lepIso::Count, std::vector<int>(0));
        // FIXME: add enum over met?
        std::vector<std::vector<int>>& map_llmet = tree["map_llmet"].write_with_init<std::vector<std::vector<int>>>(lepID::Count * lepIso::Count * lepID::Count * lepIso::Count, std::vector<int>(0));
        std::vector<std::vector<int>>& map_j = tree["map_j"].write_with_init<std::vector<std::vector<int>>>(jetID::Count * btagWP::Count, std::vector<int>(0));
        std::vector<std::vector<int>>& map_jj = tree["map_jj"].write_with_init<std::vector<std::vector<int>>>(jetID::Count * jetID::Count * btagWP::Count * btagWP::Count * jetPair::Count, std::vector<int>(0));
        std::vector<std::vector<int>>& map_llmetjj = tree["map_llmetjj"].write_with_init<std::vector<std::vector<int>>>(lepID::Count * lepIso::Count * lepID::Count * lepIso::Count * jetID::Count * jetID::Count * btagWP::Count * btagWP::Count * jetPair::Count, std::vector<int>(0));

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

        float getTriggerEfficiency(const Lepton & lep1, const Lepton & lep2)
        {
            float lep1_leg1 = 1., lep1_leg2 = 1., lep2_leg1 = 1., lep2_leg2 = 1.;
            if (lep1.isMu && lep2.isMu) {
                lep1_leg1 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
                lep1_leg2 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
                //lep1_tkLeg2 = m_hlt_efficiencies["DoubleIsoMu17Mu8_TkMu8leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
                lep2_leg1 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
                lep2_leg2 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
                //lep2_tkLeg2 = m_hlt_efficiencies["DoubleIsoMu17Mu8_TkMu8leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
            }
            else if (lep1.isMu && lep2.isEl) {
                lep1_leg1 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
                lep1_leg2 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8leg"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
                lep2_leg1 = m_hlt_efficiencies["Ele17_12Leg1"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
                lep2_leg2 = m_hlt_efficiencies["Ele17_12Leg2"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
            }
            else if (lep1.isEl && lep2.isMu) {
                lep1_leg1 = m_hlt_efficiencies["Ele17_12Leg1"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
                lep1_leg2 = m_hlt_efficiencies["Ele17_12Leg2"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
                lep2_leg1 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu17leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
                lep2_leg2 = m_hlt_efficiencies["DoubleIsoMu17Mu8_IsoMu8leg"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
            }
            else if (lep1.isEl && lep2.isEl){
                lep1_leg1 = m_hlt_efficiencies["Ele17_12Leg1"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
                lep1_leg2 = m_hlt_efficiencies["Ele17_12Leg2"].get({lep1.p4.Eta(), lep1.p4.Pt()})[0];
                lep2_leg1 = m_hlt_efficiencies["Ele17_12Leg1"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
                lep2_leg2 = m_hlt_efficiencies["Ele17_12Leg2"].get({lep2.p4.Eta(), lep2.p4.Pt()})[0];
            }
            float total_efficiency = lep1_leg1 * lep2_leg2 + lep1_leg2 * lep2_leg1 - lep1_leg1 * lep2_leg1;
            return total_efficiency;
        }

        

    private:
        // Producers name
        std::string m_electrons_producer;
        std::string m_muons_producer;
        std::string m_jets_producer;
        std::string m_met_producer;
        std::string m_nohf_met_producer;
        float m_electronIsoCut_EB_Loose, m_electronIsoCut_EE_Loose, m_electronIsoCut_EB_Tight, m_electronIsoCut_EE_Tight, m_electronEtaCut, m_leadingElectronPtCut, m_subleadingElectronPtCut;
        float m_muonLooseIsoCut, m_muonTightIsoCut, m_muonEtaCut, m_leadingMuonPtCut, m_subleadingMuonPtCut;
        float m_jetEtaCut, m_jetPtCut, m_jet_bDiscrCut_loose, m_jet_bDiscrCut_medium, m_jet_bDiscrCut_tight;
        float m_minDR_l_j_Cut;
        float m_hltDRCut, m_hltDPtCut;
        std::string m_jet_bDiscrName;
        std::string m_electron_loose_wp_name;
        std::string m_electron_tight_wp_name;
        std::map<std::string, ScaleFactor> m_hlt_efficiencies;

};

#endif
