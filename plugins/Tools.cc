#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>
#include <cp3_llbb/HHAnalysis/interface/Types.h>
#include <Math/Vector3D.h>

#define HH_HLT_DEBUG (false)

float HHAnalyzer::getCosThetaStar_CS(const LorentzVector & h1, const LorentzVector & h2, float ebeam /*= 6500*/) {
    // cos theta star angle in the Collins Soper frame
    LorentzVector p1, p2;
    p1.SetPxPyPzE(0, 0,  ebeam, ebeam);
    p2.SetPxPyPzE(0, 0, -ebeam, ebeam);

    LorentzVector hh = h1 + h2;
    ROOT::Math::Boost boost(-hh.X() / hh.T(), -hh.Y() / hh.T(), -hh.Z() / hh.T());
    p1 = boost(p1);
    p2 = boost(p2);
    LorentzVector newh1 = boost(h1);
    ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<float>> CSaxis(p1.Vect().Unit() - p2.Vect().Unit());

    return cos(ROOT::Math::VectorUtil::Angle(CSaxis.Unit(), newh1.Vect().Unit()));
}

MELAAngles HHAnalyzer::getMELAAngles(const LorentzVector &q1, const LorentzVector &q2, const LorentzVector &q11, const LorentzVector &q12, const LorentzVector &q21, const LorentzVector &q22, float ebeam /*= 6500*/) {
    MELAAngles angles;
    LorentzVector p1, p2;
    p1.SetPxPyPzE(0, 0,  ebeam, ebeam);
    p2.SetPxPyPzE(0, 0, -ebeam, ebeam);

    // set yourselves to the q1 + q2 rest frame, everyone! (note the prefix 'b' for 'boosted')
    LorentzVector X = q1 + q2;
    ROOT::Math::Boost boost(-X.X() / X.T(), -X.Y() / X.T(), -X.Z() / X.T());
    LorentzVector b_p1 = boost(p1);
    LorentzVector b_p2 = boost(p2);
    LorentzVector b_q1 = boost(q1);
    LorentzVector b_q11 = boost(q11);
    LorentzVector b_q12 = boost(q12);
    LorentzVector b_q21 = boost(q21);
    LorentzVector b_q22 = boost(q22);
    // Amend that, actually we also want some stuff in q1 or q2 rest frame
    ROOT::Math::Boost boost1(-q1.X() / q1.T(), -q1.Y() / q1.T(), -q1.Z() / q1.T());
    ROOT::Math::Boost boost2(-q2.X() / q2.T(), -q2.Y() / q2.T(), -q2.Z() / q2.T());
    LorentzVector b1_q2 = boost1(q2);
    LorentzVector b1_q11 = boost1(q11);
    LorentzVector b2_q1 = boost2(q1);
    LorentzVector b2_q21 = boost2(q21);
    // let's concentrate on three-momenta (note the prefix 'm' for three-'momenta')
    ROOT::Math::XYZVectorF m_q1(b_q1.Px(), b_q1.Py(), b_q1.Pz());
    ROOT::Math::XYZVectorF m_q11(b_q11.Px(), b_q11.Py(), b_q11.Pz());
    ROOT::Math::XYZVectorF m_q12(b_q12.Px(), b_q12.Py(), b_q12.Pz());
    ROOT::Math::XYZVectorF m_q21(b_q21.Px(), b_q21.Py(), b_q21.Pz());
    ROOT::Math::XYZVectorF m_q22(b_q22.Px(), b_q22.Py(), b_q22.Pz());
    // let's get as well three-momenta in q1 and q2 rest frame where appropriate;
    ROOT::Math::XYZVectorF m1_q2(b1_q2.Px(), b1_q2.Py(), b1_q2.Pz());
    ROOT::Math::XYZVectorF m1_q11(b1_q11.Px(), b1_q11.Py(), b1_q11.Pz());
    ROOT::Math::XYZVectorF m2_q1(b2_q1.Px(), b2_q1.Py(), b2_q1.Pz());
    ROOT::Math::XYZVectorF m2_q21(b2_q21.Px(), b2_q21.Py(), b2_q21.Pz());

    // Define reference vectors
    ROOT::Math::XYZVectorF n1(m_q11.Cross(m_q12).Unit());
    ROOT::Math::XYZVectorF n2(m_q21.Cross(m_q22).Unit());
    ROOT::Math::XYZVectorF nz(0., 0., 1.);
    ROOT::Math::XYZVectorF nsc(nz.Cross(m_q1).Unit());

    // MELA angles as taken from https://arxiv.org/pdf/1208.4018v3.pdf
    angles.phi = m_q1.Dot(n1.Cross(n2)) / fabs(m_q1.Dot(n1.Cross(n2))) * acos(- n1.Dot(n2));
    float phi1 = m_q1.Dot(n1.Cross(nsc)) / fabs(m_q1.Dot(n1.Cross(nsc))) * acos(n1.Dot(nsc));
    angles.psi = phi1 + angles.phi / 2.;
    angles.theta1 = acos(- m1_q2.Dot(m1_q11) / sqrt(m1_q2.Dot(m1_q2)) / sqrt(m1_q11.Dot(m1_q11)));
    angles.theta2 = acos(- m2_q1.Dot(m2_q21) / sqrt(m2_q1.Dot(m2_q1)) / sqrt(m2_q21.Dot(m2_q21)));
    // thetaStar is defined in the Collins-Soper frame
    ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<float>> CSaxis(b_p1.Vect().Unit() - b_p2.Vect().Unit());
    angles.thetaStar = ROOT::Math::VectorUtil::Angle(CSaxis.Unit(), b_q1.Vect().Unit());

    return angles;
}

void HHAnalyzer::matchOfflineLepton(const HLTProducer& hlt, HH::Dilepton& dilepton) {

    if (leptons[dilepton.ilep1].hlt_already_tried_matching && leptons[dilepton.ilep2].hlt_already_tried_matching) {
        if (HH_HLT_DEBUG) std::cout << "The HLT matching for this lepton pair has already been attempted, stopping here" << std::endl;
        return;
    }

    if (HH_HLT_DEBUG) {
        std::cout << "Trying to match offline leptons " << dilepton.ilep1 << " and " << dilepton.ilep2 << " (there is " << hlt.object_p4.size() << " candidate HLT objects): " << std::endl;
        std::cout   << "\tlepton1: " << (leptons[dilepton.ilep1].isMu ? "muon" : "electron")
            << " ; Pt: " << leptons[dilepton.ilep1].p4.Pt() 
            << " ; Eta: " << leptons[dilepton.ilep1].p4.Eta() 
            << " ; Phi: " << leptons[dilepton.ilep1].p4.Phi() 
            << " ; E: " << leptons[dilepton.ilep1].p4.E() 
            << std::endl;
        std::cout   << "\tlepton2: " << (leptons[dilepton.ilep2].isMu ? "muon" : "electron")
            << " ; Pt: " << leptons[dilepton.ilep2].p4.Pt() 
            << " ; Eta: " << leptons[dilepton.ilep2].p4.Eta() 
            << " ; Phi: " << leptons[dilepton.ilep2].p4.Phi() 
            << " ; E: " << leptons[dilepton.ilep2].p4.E() 
            << std::endl;
    }
    std::vector<int8_t> l1_all_indices;
    std::vector<int8_t> l2_all_indices;
    // Preselection
    for (size_t hlt_object = 0; hlt_object < hlt.object_p4.size(); hlt_object++) {
        float l1_dr = ROOT::Math::VectorUtil::DeltaR(leptons[dilepton.ilep1].p4, hlt.object_p4[hlt_object]);
        float l2_dr = ROOT::Math::VectorUtil::DeltaR(leptons[dilepton.ilep2].p4, hlt.object_p4[hlt_object]);
        float l1_dpt_over_pt = fabs(leptons[dilepton.ilep1].p4.Pt() - hlt.object_p4[hlt_object].Pt()) / leptons[dilepton.ilep1].p4.Pt();
        float l2_dpt_over_pt = fabs(leptons[dilepton.ilep2].p4.Pt() - hlt.object_p4[hlt_object].Pt()) / leptons[dilepton.ilep2].p4.Pt();
        if (HH_HLT_DEBUG && false) { // quite verbose even for debugging
                int8_t index = hlt_object;
                for (auto &path: hlt.object_paths[index])
                    std::cout << "\t# HLT path # " << +index << "\t" << path << std::endl;
                if (false) // extra verbose for further debugging
                    for (auto &filter: hlt.object_filters[index])
                        std::cout << "\t# HLT filter # " << +index << "\t" << filter << std::endl;
                std::cout << "\tPDG Id: " << hlt.object_pdg_id[index] 
                    << " ; Pt: " << hlt.object_p4[index].Pt() 
                    << " ; Eta: " << hlt.object_p4[index].Eta() 
                    << " ; Phi: " << hlt.object_p4[index].Phi() 
                    << " ; E: " << hlt.object_p4[index].E() 
                    << std::endl;
                std::cout << "\tΔR: " << l1_dr
                    << " ; ΔPt / Pt: " << l1_dpt_over_pt
                    << std::endl;
                std::cout << "\tΔR: " << l2_dr
                    << " ; ΔPt / Pt: " << l2_dpt_over_pt
                    << std::endl;
        }
        if (l1_dr < m_hltDRCut
            && l1_dpt_over_pt < m_hltDPtCut
            && ((fabs(hlt.object_pdg_id[hlt_object]) == 13 && leptons[dilepton.ilep1].isMu)
                || (fabs(hlt.object_pdg_id[hlt_object]) == 0 && leptons[dilepton.ilep1].isEl)) // It is unfortunate but the PDG ID is not correct in HLT objects
            ) {
            l1_all_indices.push_back(hlt_object);
        }
        if (l2_dr < m_hltDRCut
            && l2_dpt_over_pt < m_hltDPtCut
            && ((fabs(hlt.object_pdg_id[hlt_object]) == 13 && leptons[dilepton.ilep2].isMu)
                || (fabs(hlt.object_pdg_id[hlt_object]) == 0 && leptons[dilepton.ilep2].isEl)) // It is unfortunate but the PDG ID is not correct in HLT objects
            ) {
            l2_all_indices.push_back(hlt_object);
        }
    }
    if (l1_all_indices.empty()) {
        leptons[dilepton.ilep1].hlt_idx = -1;
        leptons[dilepton.ilep1].hlt_already_tried_matching = true;
        if (HH_HLT_DEBUG)
            std::cout << "\033[31mNo match found for first lepton\033[00m" << std::endl;
    }
    if (l2_all_indices.empty()) {
        leptons[dilepton.ilep2].hlt_idx = -1;
        leptons[dilepton.ilep2].hlt_already_tried_matching = true;
        if (HH_HLT_DEBUG)
            std::cout << "\033[31mNo match found for second lepton\033[00m" << std::endl;
    }
    // Check that the hlt path name is the same for both legs
    // FIXME: beware the day of adding single lepton HLT paths....
    std::vector<int8_t> l1_samepath_indices;
    std::vector<int8_t> l2_samepath_indices;
    for (auto& i1: l1_all_indices) {
        for (auto& i2: l2_all_indices) {
            if (i1 == i2)
                continue;
            for (auto& path1: hlt.object_paths[i1]) {
                for (auto& path2: hlt.object_paths[i2]) {
                    if (path1 == path2)
                    {
                        l1_samepath_indices.push_back(i1);
                        l2_samepath_indices.push_back(i2);
                    }
                }
            }
        }
    }
    if (l1_samepath_indices.empty()) {
        if (HH_HLT_DEBUG)
            std::cout << "\033[31mNo common HLT name match for the two leptons\033[00m" << std::endl;
        leptons[dilepton.ilep1].hlt_idx = -1;
        leptons[dilepton.ilep1].hlt_already_tried_matching = true;
        leptons[dilepton.ilep1].hlt_DR_matchedObject = std::numeric_limits<float>::max();
        leptons[dilepton.ilep1].hlt_DPtOverPt_matchedObject = std::numeric_limits<float>::max();
        leptons[dilepton.ilep2].hlt_idx = -1;
        leptons[dilepton.ilep2].hlt_already_tried_matching = true;
        leptons[dilepton.ilep2].hlt_DR_matchedObject = std::numeric_limits<float>::max();
        leptons[dilepton.ilep2].hlt_DPtOverPt_matchedObject = std::numeric_limits<float>::max();
        return;
    }
    // We have two hlt objects firing the same path: each lepton should be at least leg2, let's make sure of that
    leptons[dilepton.ilep1].hlt_leg1 = false;
    leptons[dilepton.ilep1].hlt_leg2 = false;
    leptons[dilepton.ilep2].hlt_leg1 = false;
    leptons[dilepton.ilep2].hlt_leg2 = false;
    // Check who is leg1 who is leg2
    if (
           // Di-muon
           (leptons[dilepton.ilep1].isMu && leptons[dilepton.ilep2].isMu) ||
           // Di-eletron
           (leptons[dilepton.ilep1].isEl && leptons[dilepton.ilep2].isEl)
       ) {

        std::vector<std::string> filter_leg1;
        std::vector<std::string> filter_leg2;

        static auto isLegMatched = [&hlt](const std::vector<int8_t> path_indices, const std::vector<std::string>& filters) -> bool {
            return
                std::any_of(path_indices.begin(), path_indices.end(), [&](int8_t index) {
                    for (const auto& filter: filters) {
                        if (std::any_of(hlt.object_filters[index].begin(), hlt.object_filters[index].end(), [&filter](const std::string& f) { return f == filter; }))
                            return true;
                    }

                    return false;
                });
        };

        if (leptons[dilepton.ilep1].isMu && leptons[dilepton.ilep2].isMu) {
            // di-muon filters: from the path name the legs can have asymetric cuts, taken filters from 
            // https://github.com/cms-analysis/MuonAnalysis-TagAndProbe/blob/fa1f8f3d469a5a78754ed4b4c43adbfad39a2544/python/common_variables_cff.py#L253-L264
            if (HH_HLT_DEBUG) std::cout << "\tfinding dilepton legs: di-muon" << std::endl;
            filter_leg1.push_back("hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17");
            filter_leg1.push_back("hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17");
            filter_leg2.push_back("hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8");
            filter_leg2.push_back("hltDiMuonGlbFiltered17TrkFiltered8");
            filter_leg2.push_back("hltL2pfL1sDoubleMu114ORDoubleMu125L1f0L2PreFiltered0");
        } else {
            if (HH_HLT_DEBUG) std::cout << "\tfinding dilepton legs: di-electron" << std::endl;
            filter_leg1.push_back("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
            filter_leg1.push_back("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter");
            filter_leg2.push_back("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
            filter_leg2.push_back("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter");
        }

        leptons[dilepton.ilep1].hlt_leg1 = isLegMatched(l1_samepath_indices, filter_leg1);
        leptons[dilepton.ilep1].hlt_leg2 = isLegMatched(l1_samepath_indices, filter_leg2);
        leptons[dilepton.ilep2].hlt_leg1 = isLegMatched(l2_samepath_indices, filter_leg1);
        leptons[dilepton.ilep2].hlt_leg2 = isLegMatched(l2_samepath_indices, filter_leg2);

    } else if ((leptons[dilepton.ilep1].isEl && leptons[dilepton.ilep2].isMu)
        || (leptons[dilepton.ilep1].isMu && leptons[dilepton.ilep2].isEl)) {
        // if the two offline objects are matching the same different flavour HLT path
        // then the leg1 and leg2 assignment is in sync with the order of the path name itself
        if (HH_HLT_DEBUG) std::cout << "\tfinding dilepton legs: different flavour" << std::endl;

        // Leg 1 is alway mu, and leg 2 always electron
        leptons[dilepton.ilep1].hlt_leg1 = leptons[dilepton.ilep1].isMu;
        leptons[dilepton.ilep1].hlt_leg2 = leptons[dilepton.ilep1].isEl;

        leptons[dilepton.ilep2].hlt_leg1 = !leptons[dilepton.ilep1].hlt_leg1;
        leptons[dilepton.ilep2].hlt_leg2 = !leptons[dilepton.ilep1].hlt_leg2;
    }

    if (HH_HLT_DEBUG) {
        std::cout << "\tLeg matching (before solving ambiguities):" << std::endl;
        std::cout << "\t    Lepton 1:" << std::endl;
        std::cout << "\t        Leg 1: " << std::boolalpha << leptons[dilepton.ilep1].hlt_leg1 << std::endl;
        std::cout << "\t        Leg 2: " << std::boolalpha << leptons[dilepton.ilep1].hlt_leg2 << std::endl;
        std::cout << "\t    Lepton 2:" << std::endl;
        std::cout << "\t        Leg 1: " << std::boolalpha << leptons[dilepton.ilep2].hlt_leg1 << std::endl;
        std::cout << "\t        Leg 2: " << std::boolalpha << leptons[dilepton.ilep2].hlt_leg2 << std::endl;
    }

    // Find best match wrt offline lepton
    float min_dr = std::numeric_limits<float>::max();
    float final_dpt_over_pt = std::numeric_limits<float>::max();
    int8_t index = -1;
    for (auto& i1: l1_samepath_indices) {
        float dr = ROOT::Math::VectorUtil::DeltaR(leptons[dilepton.ilep1].p4, hlt.object_p4[i1]);
        float dpt_over_pt = fabs(leptons[dilepton.ilep1].p4.Pt() - hlt.object_p4[i1].Pt()) / leptons[dilepton.ilep1].p4.Pt();
        if (dr < min_dr) {
            min_dr = dr;
            final_dpt_over_pt = dpt_over_pt;
            index = i1;
        }
    }
    leptons[dilepton.ilep1].hlt_idx = index;
    leptons[dilepton.ilep1].hlt_already_tried_matching = true;
    leptons[dilepton.ilep1].hlt_DR_matchedObject = min_dr;
    leptons[dilepton.ilep1].hlt_DPtOverPt_matchedObject = final_dpt_over_pt;
    min_dr = std::numeric_limits<float>::max();
    final_dpt_over_pt = std::numeric_limits<float>::max();
    index = -1;
    for (auto& i2: l2_samepath_indices) {
        float dr = ROOT::Math::VectorUtil::DeltaR(leptons[dilepton.ilep2].p4, hlt.object_p4[i2]);
        float dpt_over_pt = fabs(leptons[dilepton.ilep2].p4.Pt() - hlt.object_p4[i2].Pt()) / leptons[dilepton.ilep2].p4.Pt();
        if (dr < min_dr) {
            min_dr = dr;
            final_dpt_over_pt = dpt_over_pt;
            index = i2;
        }
    }
    leptons[dilepton.ilep2].hlt_idx = index;
    leptons[dilepton.ilep2].hlt_already_tried_matching = true;
    leptons[dilepton.ilep2].hlt_DR_matchedObject = min_dr;
    leptons[dilepton.ilep2].hlt_DPtOverPt_matchedObject = final_dpt_over_pt;

    // Solve ambiguities, if any
    if (!(leptons[dilepton.ilep1].hlt_leg1 && leptons[dilepton.ilep2].hlt_leg1)) {
        if (HH_HLT_DEBUG) std::cout << "\033[32mNo ambiguities! lucky day!" << std::endl;
        // let's set booleans clearly so that there remains only one leg1 and one leg2
        leptons[dilepton.ilep1].hlt_leg2 = leptons[dilepton.ilep1].hlt_leg1 ? false : true;
        leptons[dilepton.ilep2].hlt_leg2 = leptons[dilepton.ilep2].hlt_leg1 ? false : true;
    } else {
        if (HH_HLT_DEBUG) std::cout << "\033[33mOh, both leptons are 'leg1', so let's say the 'real' leg 1 is the one with highest hlt-object pt" << std::endl;
        if (hlt.object_p4[leptons[dilepton.ilep1].hlt_idx].Pt() > hlt.object_p4[leptons[dilepton.ilep2].hlt_idx].Pt()) {
            leptons[dilepton.ilep1].hlt_leg1 = true; leptons[dilepton.ilep1].hlt_leg2 = false;
            leptons[dilepton.ilep2].hlt_leg1 = false; leptons[dilepton.ilep2].hlt_leg2 = true;
        } else {
            leptons[dilepton.ilep1].hlt_leg1 = false; leptons[dilepton.ilep1].hlt_leg2 = true;
            leptons[dilepton.ilep2].hlt_leg1 = true; leptons[dilepton.ilep2].hlt_leg2 = false;
        }
    }
    if (HH_HLT_DEBUG) {
        index = leptons[dilepton.ilep1].hlt_idx;
        std::cout << "\033[32mLepton 1 is leg " << (leptons[dilepton.ilep1].hlt_leg1 ? 1 : 2) << " and matched with online object:\033[00m" << std::endl;
        std::cout << "\tHLT paths:" << std::endl;
        for (auto &path: hlt.object_paths[index])
            std::cout << "\t\t" << path << std::endl;
        std::cout << "\tPDG Id: " << hlt.object_pdg_id[index] 
            << " ; Pt: " << hlt.object_p4[index].Pt() 
            << " ; Eta: " << hlt.object_p4[index].Eta() 
            << " ; Phi: " << hlt.object_p4[index].Phi() 
            << " ; E: " << hlt.object_p4[index].E() 
            << std::endl;
        std::cout << "\tΔR: " << leptons[dilepton.ilep1].hlt_DR_matchedObject
            << " ; ΔPt / Pt: " << leptons[dilepton.ilep1].hlt_DPtOverPt_matchedObject
            << std::endl;
        index = leptons[dilepton.ilep2].hlt_idx;
        std::cout << "\033[32mLepton 2 is leg " << (leptons[dilepton.ilep2].hlt_leg1 ? 1 : 2) << " and matched with online object:\033[00m" << std::endl;
        std::cout << "\tHLT paths:" << std::endl;
        for (auto &path: hlt.object_paths[index])
            std::cout << "\t\t" << path << std::endl;
        std::cout << "\tPDG Id: " << hlt.object_pdg_id[index] 
            << " ; Pt: " << hlt.object_p4[index].Pt() 
            << " ; Eta: " << hlt.object_p4[index].Eta() 
            << " ; Phi: " << hlt.object_p4[index].Phi() 
            << " ; E: " << hlt.object_p4[index].E() 
            << std::endl;
        std::cout << "\tΔR: " << leptons[dilepton.ilep2].hlt_DR_matchedObject
            << " ; ΔPt / Pt: " << leptons[dilepton.ilep2].hlt_DPtOverPt_matchedObject
            << std::endl;
    }
}

void HHAnalyzer::fillTriggerEfficiencies(const Lepton & lep1, const Lepton & lep2, Dilepton & dilep) {

    float eff_lep1_leg1 = 1.;
    float eff_lep1_leg2 = 1.;
    float eff_lep2_leg1 = 1.;
    float eff_lep2_leg2 = 1.;
    Parameters p_hlt_lep1 = {{BinningVariable::Eta, lep1.p4.Eta()}, {BinningVariable::Pt, lep1.p4.Pt()}};
    Parameters p_hlt_lep2 = {{BinningVariable::Eta, lep2.p4.Eta()}, {BinningVariable::Pt, lep2.p4.Pt()}};

    if (lep1.isMu && lep2.isMu) {
        eff_lep1_leg1 = m_hlt_efficiencies.at("IsoMu17leg")->get(p_hlt_lep1)[0];
        eff_lep1_leg2 = m_hlt_efficiencies.at("IsoMu8orIsoTkMu8leg")->get(p_hlt_lep1)[0];
        eff_lep2_leg1 = m_hlt_efficiencies.at("IsoMu17leg")->get(p_hlt_lep2)[0];
        eff_lep2_leg2 = m_hlt_efficiencies.at("IsoMu8orIsoTkMu8leg")->get(p_hlt_lep2)[0];
    }
    else if (lep1.isMu && lep2.isEl) {
        eff_lep1_leg1 = m_hlt_efficiencies.at("IsoMu23leg")->get(p_hlt_lep1)[0];
        eff_lep1_leg2 = m_hlt_efficiencies.at("IsoMu8leg")->get(p_hlt_lep1)[0];
        eff_lep2_leg1 = m_hlt_efficiencies.at("EleMuHighPtleg")->get(p_hlt_lep2)[0];
        eff_lep2_leg2 = m_hlt_efficiencies.at("MuEleLowPtleg")->get(p_hlt_lep2)[0];
    }
    else if (lep1.isEl && lep2.isMu) {
        eff_lep1_leg1 = m_hlt_efficiencies.at("EleMuHighPtleg")->get(p_hlt_lep1)[0];
        eff_lep1_leg2 = m_hlt_efficiencies.at("MuEleLowPtleg")->get(p_hlt_lep1)[0];
        eff_lep2_leg1 = m_hlt_efficiencies.at("IsoMu23leg")->get(p_hlt_lep2)[0];
        eff_lep2_leg2 = m_hlt_efficiencies.at("IsoMu8leg")->get(p_hlt_lep2)[0];
    }
    else if (lep1.isEl && lep2.isEl){
        eff_lep1_leg1 = m_hlt_efficiencies.at("DoubleEleHighPtleg")->get(p_hlt_lep1)[0];
        eff_lep1_leg2 = m_hlt_efficiencies.at("DoubleEleLowPtleg")->get(p_hlt_lep1)[0];
        eff_lep2_leg1 = m_hlt_efficiencies.at("DoubleEleHighPtleg")->get(p_hlt_lep2)[0];
        eff_lep2_leg2 = m_hlt_efficiencies.at("DoubleEleLowPtleg")->get(p_hlt_lep2)[0];
    }
    else 
        std::cout << "We have something else then el or mu !!" << std::endl;

    float error_eff_lep1_leg1_up = 0.;
    float error_eff_lep1_leg2_up = 0.;
    float error_eff_lep2_leg1_up = 0.;
    float error_eff_lep2_leg2_up = 0.;

    if (lep1.isMu && lep2.isMu) {
        error_eff_lep1_leg1_up = m_hlt_efficiencies.at("IsoMu17leg")->get(p_hlt_lep1)[2];
        error_eff_lep1_leg2_up = m_hlt_efficiencies.at("IsoMu8orIsoTkMu8leg")->get(p_hlt_lep1)[2];
        error_eff_lep2_leg1_up = m_hlt_efficiencies.at("IsoMu17leg")->get(p_hlt_lep2)[2];
        error_eff_lep2_leg2_up = m_hlt_efficiencies.at("IsoMu8orIsoTkMu8leg")->get(p_hlt_lep2)[2];
    }
    else if (lep1.isMu && lep2.isEl) {
        error_eff_lep1_leg1_up = m_hlt_efficiencies.at("IsoMu23leg")->get(p_hlt_lep1)[2];
        error_eff_lep1_leg2_up = m_hlt_efficiencies.at("IsoMu8leg")->get(p_hlt_lep1)[2];
        error_eff_lep2_leg1_up = m_hlt_efficiencies.at("EleMuHighPtleg")->get(p_hlt_lep2)[2];
        error_eff_lep2_leg2_up = m_hlt_efficiencies.at("MuEleLowPtleg")->get(p_hlt_lep2)[2];
    }
    else if (lep1.isEl && lep2.isMu) {
        error_eff_lep1_leg1_up = m_hlt_efficiencies.at("EleMuHighPtleg")->get(p_hlt_lep1)[2];
        error_eff_lep1_leg2_up = m_hlt_efficiencies.at("MuEleLowPtleg")->get(p_hlt_lep1)[2];
        error_eff_lep2_leg1_up = m_hlt_efficiencies.at("IsoMu23leg")->get(p_hlt_lep2)[2];
        error_eff_lep2_leg2_up = m_hlt_efficiencies.at("IsoMu8leg")->get(p_hlt_lep2)[2];
    }
    else if (lep1.isEl && lep2.isEl){
        error_eff_lep1_leg1_up = m_hlt_efficiencies.at("DoubleEleHighPtleg")->get(p_hlt_lep1)[2];
        error_eff_lep1_leg2_up = m_hlt_efficiencies.at("DoubleEleLowPtleg")->get(p_hlt_lep1)[2];
        error_eff_lep2_leg1_up = m_hlt_efficiencies.at("DoubleEleHighPtleg")->get(p_hlt_lep2)[2];
        error_eff_lep2_leg2_up = m_hlt_efficiencies.at("DoubleEleLowPtleg")->get(p_hlt_lep2)[2];
    }

    float error_eff_lep1_leg1_down = 0.;
    float error_eff_lep1_leg2_down = 0.;
    float error_eff_lep2_leg1_down = 0.;
    float error_eff_lep2_leg2_down = 0.;

    if (lep1.isMu && lep2.isMu) {
        error_eff_lep1_leg1_down = m_hlt_efficiencies.at("IsoMu17leg")->get(p_hlt_lep1)[1];
        error_eff_lep1_leg2_down = m_hlt_efficiencies.at("IsoMu8orIsoTkMu8leg")->get(p_hlt_lep1)[1];
        error_eff_lep2_leg1_down = m_hlt_efficiencies.at("IsoMu17leg")->get(p_hlt_lep2)[1];
        error_eff_lep2_leg2_down = m_hlt_efficiencies.at("IsoMu8orIsoTkMu8leg")->get(p_hlt_lep2)[1];
    }
    else if (lep1.isMu && lep2.isEl) {
        error_eff_lep1_leg1_down = m_hlt_efficiencies.at("IsoMu23leg")->get(p_hlt_lep1)[1];
        error_eff_lep1_leg2_down = m_hlt_efficiencies.at("IsoMu8leg")->get(p_hlt_lep1)[1];
        error_eff_lep2_leg1_down = m_hlt_efficiencies.at("EleMuHighPtleg")->get(p_hlt_lep2)[1];
        error_eff_lep2_leg2_down = m_hlt_efficiencies.at("MuEleLowPtleg")->get(p_hlt_lep2)[1];
    }
    else if (lep1.isEl && lep2.isMu) {
        error_eff_lep1_leg1_down = m_hlt_efficiencies.at("EleMuHighPtleg")->get(p_hlt_lep1)[1];
        error_eff_lep1_leg2_down = m_hlt_efficiencies.at("MuEleLowPtleg")->get(p_hlt_lep1)[1];
        error_eff_lep2_leg1_down = m_hlt_efficiencies.at("IsoMu23leg")->get(p_hlt_lep2)[1];
        error_eff_lep2_leg2_down = m_hlt_efficiencies.at("IsoMu8leg")->get(p_hlt_lep2)[1];
    }
    else if (lep1.isEl && lep2.isEl){
        error_eff_lep1_leg1_down = m_hlt_efficiencies.at("DoubleEleHighPtleg")->get(p_hlt_lep1)[1];
        error_eff_lep1_leg2_down = m_hlt_efficiencies.at("DoubleEleLowPtleg")->get(p_hlt_lep1)[1];
        error_eff_lep2_leg1_down = m_hlt_efficiencies.at("DoubleEleHighPtleg")->get(p_hlt_lep2)[1];
        error_eff_lep2_leg2_down = m_hlt_efficiencies.at("DoubleEleLowPtleg")->get(p_hlt_lep2)[1];
    }


    float nominal = -(eff_lep1_leg1 * eff_lep2_leg1) +
        (1 - (1 - eff_lep1_leg2)) * eff_lep2_leg1 +
        eff_lep1_leg1 * (1 - (1 - eff_lep2_leg2));

    float error_squared_up =
        std::pow(1 - eff_lep2_leg1 - (1 - eff_lep2_leg2), 2) *
        std::pow(error_eff_lep1_leg1_up, 2) +
        std::pow(eff_lep2_leg1, 2) *
        std::pow(error_eff_lep1_leg2_up, 2) +
        std::pow(1 - eff_lep1_leg1 - (1 - eff_lep1_leg2), 2) *
        std::pow(error_eff_lep2_leg1_up, 2) +
        std::pow(eff_lep1_leg1, 2) *
        std::pow(error_eff_lep2_leg2_up, 2);

    float error_squared_down = 
        std::pow(1 - eff_lep2_leg1 - (1 - eff_lep2_leg2), 2) *
        std::pow(error_eff_lep1_leg1_down, 2) +
        std::pow(eff_lep2_leg1, 2) *
        std::pow(error_eff_lep1_leg2_down, 2) +
        std::pow(1 - eff_lep1_leg1 - (1 - eff_lep1_leg2), 2) *
        std::pow(error_eff_lep2_leg1_down, 2) +
        std::pow(eff_lep1_leg1, 2) *
        std::pow(error_eff_lep2_leg2_down, 2);

    dilep.trigger_efficiency = nominal;
    dilep.trigger_efficiency_upVariated = ((nominal + std::sqrt(error_squared_up)) > 1.)? 1. : nominal + std::sqrt(error_squared_up);
    dilep.trigger_efficiency_downVariated = ((nominal - std::sqrt(error_squared_down)) < 0.)? 0. : nominal - std::sqrt(error_squared_down);
    
    // Arun's method (not using the proper derivative formula)
    float X = eff_lep1_leg1 * eff_lep2_leg2 * std::sqrt((std::pow((error_eff_lep1_leg1_up/eff_lep1_leg1),2) + std::pow((error_eff_lep2_leg2_up/eff_lep2_leg2),2) ));
    float Y = eff_lep2_leg1 * eff_lep1_leg2 * std::sqrt((std::pow((error_eff_lep2_leg1_up/eff_lep2_leg1),2) + std::pow((error_eff_lep1_leg2_up/eff_lep1_leg2),2) ));
    float Z = eff_lep1_leg1 * eff_lep2_leg1 * std::sqrt((std::pow((error_eff_lep1_leg1_up/eff_lep1_leg1),2) + std::pow((error_eff_lep2_leg1_up/eff_lep2_leg1),2) ));
    float error_squared_up_Arun = X*X + Y*Y + Z*Z ;
    dilep.trigger_efficiency_upVariated_Arun = ((nominal + std::sqrt(error_squared_up_Arun)) > 1.)? 1. : nominal + std::sqrt(error_squared_up_Arun);

    X = eff_lep1_leg1 * eff_lep2_leg2 * std::sqrt((std::pow((error_eff_lep1_leg1_down/eff_lep1_leg1),2) + std::pow((error_eff_lep2_leg2_down/eff_lep2_leg2),2) ));
    Y = eff_lep2_leg1 * eff_lep1_leg2 * std::sqrt((std::pow((error_eff_lep2_leg1_down/eff_lep2_leg1),2) + std::pow((error_eff_lep1_leg2_down/eff_lep1_leg2),2) ));
    Z = eff_lep1_leg1 * eff_lep2_leg1 * std::sqrt((std::pow((error_eff_lep1_leg1_down/eff_lep1_leg1),2) + std::pow((error_eff_lep2_leg1_down/eff_lep2_leg1),2) ));
    float error_squared_down_Arun = X*X + Y*Y + Z*Z ;
    dilep.trigger_efficiency_downVariated_Arun = ((nominal - std::sqrt(error_squared_down_Arun)) < 0.)? 0. : nominal - std::sqrt(error_squared_down_Arun);

}


