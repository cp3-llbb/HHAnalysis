#pragma once

#include <vector>
#include <Math/Vector4D.h>
#include <cp3_llbb/HHAnalysis/interface/Indices.h>

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>> LorentzVector;

namespace HH {
    struct Lepton {
        LorentzVector p4;
        LorentzVector gen_p4;
        int8_t charge;
        int idx;
        int8_t hlt_idx = -1; // Index to the matched HLT object. -1 if no match. 
                             // Example : t->Draw("hh_leptons.p4.Pt() - hlt_object_p4[hh_leptons.hlt_idx].Pt()","hh_leptons.hlt_idx != -1","")
        bool hlt_already_tried_matching = false; // do the matching only once, even when the lepton is in several Dilepton
        float hlt_DR_matchedObject = std::numeric_limits<float>::max();
        float hlt_DPtOverPt_matchedObject = std::numeric_limits<float>::max();
        bool isMu;
        bool isEl;
        bool id_L; // Loose
        bool id_T; // Tight
        bool iso_L; // Loose
        bool iso_T; // Tight
        bool gen_matched;
        float gen_DR;
        float gen_DPtOverPt;
    };
    struct Dilepton {
        LorentzVector p4;
        LorentzVector gen_p4;
        std::pair<int, int> idxs; // indices in the corresponding framework collection
        int ilep1; // index in the HH::Lepton collection
        int ilep2; // index in the HH::Lepton collection
        std::pair<int8_t, int8_t> hlt_idxs = std::make_pair(-1,-1); // Stores indices of matched online objects. (-1,-1) if no match
        bool isOS; // Opposite Sign
        bool isMuMu;
        bool isElEl;
        bool isElMu;
        bool isMuEl;
        bool isSF; // Same Flavour
        bool id_LL;
        bool id_LT;
        bool id_TL;
        bool id_TT;
        bool iso_LL;
        bool iso_LT;
        bool iso_TL;
        bool iso_TT;
        float DR_l_l;
        float DPhi_l_l;
        bool gen_matched;
        float gen_DR;
        float gen_DPtOverPt;
    };
    struct Met {
        LorentzVector p4;
        LorentzVector gen_p4;
        bool isNoHF;
        bool gen_matched;
        float gen_DR;
        float gen_DPhi;
        float gen_DPtOverPt;
    };
    struct DileptonMet : public Dilepton, public Met {
        LorentzVector p4;
        LorentzVector gen_p4;
        int ill; // index in the HH::Dilepton collection
        int imet; // index in the HH::Met collection
        float DPhi_ll_met;
        float minDPhi_l_met;
        float maxDPhi_l_met;
        float MT;
        float MT_formula;
        float projectedMet;
        bool gen_matched;
        float gen_DR;
        float gen_DPhi;
        float gen_DPtOverPt;
    };
    struct Jet {
        LorentzVector p4;
        LorentzVector gen_p4;
        int idx;
        bool id_L;
        bool id_T;
        bool id_TLV;
        bool btag_L;
        bool btag_M;
        bool btag_T;
        float CSV;
        float JP;
        bool gen_matched_bParton;
        bool gen_matched_bHadron;
        bool gen_matched;
        float gen_DR;
        float gen_DPtOverPt;
    };
    struct Dijet {
        LorentzVector p4;
        LorentzVector gen_p4;
        std::pair<int, int> idxs; // indices in the framework collection
        int ijet1; // indices in the HH::Jet collection
        int ijet2;
        bool btag_LL;
        bool btag_LM;
        bool btag_LT;
        bool btag_ML;
        bool btag_MM;
        bool btag_MT;
        bool btag_TL;
        bool btag_TM;
        bool btag_TT;
        float sumCSV;
        float sumJP;
        float DR_j_j;
        float DPhi_j_j;
        bool gen_matched_bbPartons;
        bool gen_matched_bbHadrons;
        bool gen_matched;
        float gen_DR;
        float gen_DPtOverPt;
    };
    struct DileptonMetDijet : public DileptonMet, public Dijet {
        LorentzVector p4;
        LorentzVector lep1_p4;
        LorentzVector lep2_p4;
        LorentzVector jet1_p4;
        LorentzVector jet2_p4;
        LorentzVector met_p4;
        LorentzVector ll_p4;
        LorentzVector jj_p4;
        LorentzVector lljj_p4;
        LorentzVector gen_p4;
        LorentzVector gen_lep1_p4;
        LorentzVector gen_lep2_p4;
        LorentzVector gen_jet1_p4;
        LorentzVector gen_jet2_p4;
        LorentzVector gen_met_p4;
        LorentzVector gen_ll_p4;
        LorentzVector gen_jj_p4;
        LorentzVector gen_lljj_p4;
        int illmet; // index in the HH::DileptonMet collection
        int ijj; // index in the HH::Dijet collection
        float DPhi_jj_met;
        float minDPhi_j_met;
        float maxDPhi_j_met;
        float maxDR_l_j;
        float minDR_l_j;
        float DR_ll_jj;
        float DPhi_ll_jj;
        float DR_llmet_jj;
        float DPhi_llmet_jj;
        float cosThetaStar_CS;
        bool gen_matched;
        float gen_DR;
        float gen_DPhi;
        float gen_DPtOverPt;
    };
}


