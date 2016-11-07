#ifndef HHANALYZER_H
#define HHANALYZER_H

#include <cp3_llbb/Framework/interface/Analyzer.h>
#include <cp3_llbb/Framework/interface/Category.h>
#include <cp3_llbb/Framework/interface/BinnedValuesJSONParser.h>
#include <cp3_llbb/Framework/interface/WeightedBinnedValues.h>

#include <cp3_llbb/HHAnalysis/interface/Types.h>
#include <cp3_llbb/HHAnalysis/interface/lester_mt2_bisect.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

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
            m_electron_medium_wp_name = config.getUntrackedParameter<std::string>("electrons_medium_wp_name", "cutBasedElectronID-Spring15-50ns-V1-standalone-medium");
            m_electron_tight_wp_name = config.getUntrackedParameter<std::string>("electrons_tight_wp_name", "cutBasedElectronID-Spring15-50ns-V1-standalone-tight");

            m_jetEtaCut = config.getUntrackedParameter<double>("jetEtaCut", 2.4);
            m_jetPtCut = config.getUntrackedParameter<double>("jetPtCut", 20);
            m_jet_bDiscrName = config.getUntrackedParameter<std::string>("discr_name", "pfCombinedInclusiveSecondaryVertexV2BJetTags");
            m_jet_bDiscrCut_loose = config.getUntrackedParameter<double>("discr_cut_loose", 0.605);
            m_jet_bDiscrCut_medium = config.getUntrackedParameter<double>("discr_cut_medium", 0.89);
            m_jet_bDiscrCut_tight = config.getUntrackedParameter<double>("discr_cut_tight", 0.97);
            m_minDR_l_j_Cut = config.getUntrackedParameter<double>("minDR_l_j_Cut", 0.3);
            m_applyBJetRegression = config.getUntrackedParameter<bool>("applyBJetRegression", false);

            m_hltDRCut = config.getUntrackedParameter<double>("hltDRCut", std::numeric_limits<float>::max());
            m_hltDPtCut = config.getUntrackedParameter<double>("hltDPtCut", std::numeric_limits<float>::max());

            const edm::ParameterSet& hlt_efficiencies = config.getUntrackedParameter<edm::ParameterSet>("hlt_efficiencies");
            std::vector<std::string> hlt_efficiencies_name = hlt_efficiencies.getParameterNames();
            for (const std::string& hlt_efficiency: hlt_efficiencies_name) {
                std::cout << "    Registering new HLT efficiency: " << hlt_efficiency;
                if (hlt_efficiencies.existsAs<edm::FileInPath>(hlt_efficiency, false)) {
                    BinnedValuesJSONParser parser(hlt_efficiencies.getUntrackedParameter<edm::FileInPath>(hlt_efficiency).fullPath());
                    m_hlt_efficiencies.emplace(hlt_efficiency, std::unique_ptr<BinnedValues>(new BinnedValues(std::move(parser.get_values()))));
                    std::cout << " -> non-weighted. " << std::endl;
                } else {
                    const auto& parts = hlt_efficiencies.getUntrackedParameter<std::vector<edm::ParameterSet>>(hlt_efficiency);
                    m_hlt_efficiencies.emplace(hlt_efficiency, std::unique_ptr<BinnedValues>(new WeightedBinnedValues(parts)));
                    std::cout << " -> weighted. " << std::endl;
                }
            }

            asymm_mt2_lester_bisect::disableCopyrightMessage();
        }
        virtual void endJob(MetadataManager&) override;

        // leptons and dileptons stuff
        BRANCH(electrons, std::vector<unsigned int>);
        BRANCH(muons, std::vector<unsigned int>);
        BRANCH(leptons, std::vector<HH::Lepton>);
        BRANCH(met, std::vector<HH::Met>);
        BRANCH(jets, std::vector<HH::Jet>);
        std::vector<HH::Dilepton> ll;
        std::vector<HH::DileptonMet> llmet;
        std::vector<HH::Dijet> jj;
        std::vector<HH::DileptonMetDijet> llmetjj;
        std::vector<HH::DileptonMetDijet> llmetjj_csv;
        // some few custom candidates, for convenience
        // Januray 2016: preapproval freezing custom candidates
        BRANCH(llmetjj_HWWleptons_nobtag_csv, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_HWWleptons_btagL_csv, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_HWWleptons_btagM_csv, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_HWWleptons_btagT_csv, std::vector<HH::DileptonMetDijet>);
        // October 2016: adding some asymmetric btag candidates, for study
        BRANCH(llmetjj_HWWleptons_btagLM_csv, std::vector<HH::DileptonMetDijet>);
        BRANCH(llmetjj_HWWleptons_btagMT_csv, std::vector<HH::DileptonMetDijet>);

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const AnalyzersManager&, const CategoryManager&) override;
        virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet& config) override;

        // Various helper functions, implemented in plugins/Tools.cc
        float getCosThetaStar_CS(const LorentzVector & h1, const LorentzVector & h2, float ebeam = 6500);
        MELAAngles getMELAAngles(const LorentzVector &q1, const LorentzVector &q2, const LorentzVector &q11, const LorentzVector &q12, const LorentzVector &q21, const LorentzVector &q22, float ebeam = 6500);
        void matchOfflineLepton(const HLTProducer& hlt, Dilepton& dilepton);
        void fillTriggerEfficiencies(const Lepton & lep1, const Lepton & lep2, Dilepton & dilep);

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

        float count_has2leptons = 0.;
        float count_has2leptons_elel = 0.;
        float count_has2leptons_elmu = 0.;
        float count_has2leptons_muel = 0.;
        float count_has2leptons_mumu = 0.;
        float count_has2leptons_1llmetjj = 0.;
        float count_has2leptons_elel_1llmetjj = 0.;
        float count_has2leptons_elmu_1llmetjj = 0.;
        float count_has2leptons_muel_1llmetjj = 0.;
        float count_has2leptons_mumu_1llmetjj = 0.;
        float count_has2leptons_1llmetjj_2btagM = 0.;
        float count_has2leptons_elel_1llmetjj_2btagM = 0.;
        float count_has2leptons_elmu_1llmetjj_2btagM = 0.;
        float count_has2leptons_muel_1llmetjj_2btagM = 0.;
        float count_has2leptons_mumu_1llmetjj_2btagM = 0.;

        // ttbar system mc truth
        // Gen matching. All indexes are from the `pruned` collection
        uint16_t gen_t; // Index of the top quark
        uint16_t gen_t_beforeFSR; // Index of the top quark, before any FSR
        uint16_t gen_tbar; // Index of the anti-top quark
        uint16_t gen_tbar_beforeFSR; // Index of the anti-top quark, before any FSR

        uint16_t gen_b; // Index of the b quark coming from the top decay
        uint16_t gen_b_beforeFSR; // Index of the b quark coming from the top decay, before any FSR
        uint16_t gen_bbar; // Index of the anti-b quark coming from the anti-top decay
        uint16_t gen_bbar_beforeFSR; // Index of the anti-b quark coming from the anti-top decay, before any FSR

        uint16_t gen_jet1_t; // Index of the first jet from the top decay chain
        uint16_t gen_jet1_t_beforeFSR; // Index of the first jet from the top decay chain, before any FSR
        uint16_t gen_jet2_t; // Index of the second jet from the top decay chain
        uint16_t gen_jet2_t_beforeFSR; // Index of the second jet from the top decay chain, before any FSR

        uint16_t gen_jet1_tbar; // Index of the first jet from the anti-top decay chain
        uint16_t gen_jet1_tbar_beforeFSR; // Index of the first jet from the anti-top decay chain, before any FSR
        uint16_t gen_jet2_tbar; // Index of the second jet from the anti-top decay chain
        uint16_t gen_jet2_tbar_beforeFSR; // Index of the second jet from the anti-top decay chain, before any FSR

        uint16_t gen_lepton_t; // Index of the lepton from the top decay chain
        uint16_t gen_lepton_t_beforeFSR; // Index of the lepton from the top decay chain, before any FSR
        uint16_t gen_neutrino_t; // Index of the neutrino from the top decay chain
        uint16_t gen_neutrino_t_beforeFSR; // Index of the neutrino from the top decay chain, before any FSR

        uint16_t gen_lepton_tbar; // Index of the lepton from the anti-top decay chain
        uint16_t gen_lepton_tbar_beforeFSR; // Index of the lepton from the anti-top decay chain, before any FSR
        uint16_t gen_neutrino_tbar; // Index of the neutrino from the anti-top decay chain
        uint16_t gen_neutrino_tbar_beforeFSR; // Index of the neutrino from the anti-top decay chain, before any FSR

        BRANCH(gen_ttbar_decay_type, char); // Type of ttbar decay. Can take any values from TTDecayType enum

        // Di-higgs gen system
        BRANCH(gen_iX, char);
        BRANCH(gen_X, LorentzVector);

        BRANCH(gen_iH1, char);
        BRANCH(gen_iH2, char);
        BRANCH(gen_H1, LorentzVector);
        BRANCH(gen_H2, LorentzVector);
        BRANCH(gen_mHH, double);
        BRANCH(gen_costhetastar, double);
        BRANCH(gen_iH1_afterFSR, char);
        BRANCH(gen_iH2_afterFSR, char);
        BRANCH(gen_H1_afterFSR, LorentzVector);
        BRANCH(gen_H2_afterFSR, LorentzVector);

        BRANCH(gen_iB, char);
        BRANCH(gen_iBbar, char);
        BRANCH(gen_iB_afterFSR, char);
        BRANCH(gen_iBbar_afterFSR, char);
        BRANCH(gen_B, LorentzVector);
        BRANCH(gen_Bbar, LorentzVector);
        BRANCH(gen_B_afterFSR, LorentzVector);
        BRANCH(gen_Bbar_afterFSR, LorentzVector);

        BRANCH(gen_iV2, char);
        BRANCH(gen_iV1, char);
        BRANCH(gen_iV2_afterFSR, char);
        BRANCH(gen_iV1_afterFSR, char);
        BRANCH(gen_V2, LorentzVector);
        BRANCH(gen_V1, LorentzVector);
        BRANCH(gen_V2_afterFSR, LorentzVector);
        BRANCH(gen_V1_afterFSR, LorentzVector);

        BRANCH(gen_iLminus, char);
        BRANCH(gen_iLplus, char);
        BRANCH(gen_iLminus_afterFSR, char);
        BRANCH(gen_iLplus_afterFSR, char);
        BRANCH(gen_Lminus, LorentzVector);
        BRANCH(gen_Lplus, LorentzVector);
        BRANCH(gen_Lminus_afterFSR, LorentzVector);
        BRANCH(gen_Lplus_afterFSR, LorentzVector);

        BRANCH(gen_iNu1, char);
        BRANCH(gen_iNu2, char);
        BRANCH(gen_Nu1, LorentzVector);
        BRANCH(gen_Nu2, LorentzVector);

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
        std::string m_electron_medium_wp_name;
        std::string m_electron_tight_wp_name;
        bool m_applyBJetRegression;
        std::unordered_map<std::string, std::unique_ptr<BinnedValues>> m_hlt_efficiencies;

};

// Some macros for gen information
#define ASSIGN_HH_GEN_INFO_2(X, Y, ERROR) \
        /* Before FSR. Use isHardProcess (7) flag */ \
        if (flags.test(7)) { \
            if (gen_i##X == -1) { \
                gen_i##X = ip; \
                gen_##X = p4; \
            } else if (gen_i##Y == -1) { \
                gen_i##Y = ip; \
                gen_##Y = p4; \
            } else { \
                std::cout << "Warning: more than two " #ERROR " in the hard process" << std::endl; \
            } \
        /* After FSR. Use isLastCopy (13) flag */ \
        } else if (flags.test(13)) { \
            if (gen_i##X##_afterFSR == -1) { \
                gen_i##X##_afterFSR = ip; \
                gen_##X##_afterFSR = p4; \
            } else if (gen_i##Y##_afterFSR == -1) { \
                gen_i##Y##_afterFSR = ip; \
                gen_##Y##_afterFSR = p4; \
            } else { \
                std::cout << "Warning: more than two " #ERROR " after FSR" << std::endl; \
            } \
        }

#define ASSIGN_HH_GEN_INFO_2_NO_FSR(X, Y, ERROR) \
        /* Before FSR. Use isHardProcess (7) flag */ \
        if (flags.test(7)) { \
            if (gen_i##X == -1) { \
                gen_i##X = ip; \
                gen_##X = p4; \
            } else if (gen_i##Y == -1) { \
                gen_i##Y = ip; \
                gen_##Y = p4; \
            } else { \
                std::cout << "Warning: more than two " #ERROR " in the hard process" << std::endl; \
            } \
        }

#define ASSIGN_HH_GEN_INFO(X, ERROR) \
        /* Before FSR. Use isHardProcess (7) flag */ \
        if (flags.test(7)) { \
            if (gen_i##X == -1) { \
                gen_i##X = ip; \
                gen_##X = p4; \
            } else { \
                std::cout << "Warning: more than one " #ERROR " in the hard process" << std::endl; \
            } \
        /* After FSR. Use isLastCopy (13) flag */ \
        } else if (flags.test(13)) { \
            if (gen_i##X##_afterFSR == -1) { \
                gen_i##X##_afterFSR = ip; \
                gen_##X##_afterFSR = p4; \
            } else { \
                std::cout << "Warning: more than one " #ERROR " after FSR" << std::endl; \
            } \
        }

#define ASSIGN_HH_GEN_INFO_NO_FSR(X, ERROR) \
        /* Before FSR. Use isHardProcess (7) flag */ \
        if (flags.test(7)) { \
            if (gen_i##X == -1) { \
                gen_i##X = ip; \
                gen_##X = p4; \
            } else { \
                std::cout << "Warning: more than one " #ERROR " in the hard process" << std::endl; \
            } \
        }

#define PRINT_PARTICULE(X) \
        if (gen_i##X != -1) { \
            std::cout << "    gen_" #X ".M() = " << gen_##X.M() << std::endl; \
        }

#define PRINT_RESONANCE(X, Y) \
        if ((gen_i##X != -1) && (gen_i##Y != -1)) { \
            std::cout << "    gen_" #X ".M() = " << gen_##X.M() << std::endl; \
            std::cout << "    gen_" #Y ".M() = " << gen_##Y.M() << std::endl; \
            std::cout << "    gen_(" #X " + " #Y ").M() = " << (gen_##X + gen_##Y).M() << std::endl; \
        } \
        if ((gen_i##X##_afterFSR != -1) && (gen_i##Y##_afterFSR != -1)) { \
            std::cout << "    gen_" #X "_afterFSR.M() = " << gen_##X##_afterFSR.M() << std::endl; \
            std::cout << "    gen_" #Y "_afterFSR.M() = " << gen_##Y##_afterFSR.M() << std::endl; \
            std::cout << "    gen_(" #X " + " #Y ")_afterFSR.M() = " << (gen_##X##_afterFSR + gen_##Y##_afterFSR).M() << std::endl; \
        }

#define PRINT_RESONANCE_NO_FSR(X, Y) \
        if ((gen_i##X != -1) && (gen_i##Y != -1)) { \
            std::cout << "    gen_" #X ".M() = " << gen_##X.M() << std::endl; \
            std::cout << "    gen_" #Y ".M() = " << gen_##Y.M() << std::endl; \
            std::cout << "    gen_(" #X " + " #Y ").M() = " << (gen_##X + gen_##Y).M() << std::endl; \
        }


#endif
