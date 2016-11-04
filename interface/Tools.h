#pragma once

#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>
#include <cp3_llbb/HHAnalysis/interface/Types.h>
#include <Math/Vector3D.h>

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

void HHAnalyzer::fillTriggerEfficiencies(const Lepton & lep1, const Lepton & lep2, Dilepton & dilep) {

    float eff_lep1_leg1 = 1.;
    float eff_lep1_leg2 = 1.;
    float eff_lep1_tkleg2 = 0.;
    float eff_lep2_leg1 = 1.;
    float eff_lep2_leg2 = 1.;
    float eff_lep2_tkleg2 = 0.;
    Parameters p_hlt_lep1 = {{BinningVariable::Eta, lep1.p4.Eta()}, {BinningVariable::Pt, lep1.p4.Pt()}};
    Parameters p_hlt_lep2 = {{BinningVariable::Eta, lep2.p4.Eta()}, {BinningVariable::Pt, lep2.p4.Pt()}};

    if (lep1.isMu && lep2.isMu) {
        eff_lep1_leg1 = m_hlt_efficiencies.at("IsoMu17leg").get(p_hlt_lep1)[0];
        eff_lep1_leg2 = m_hlt_efficiencies.at("IsoMu8orIsoTkMu8leg").get(p_hlt_lep1)[0];
        //eff_lep1_tkleg2 = m_hlt_efficiencies.at("DoubleIsoMu17Mu8_TkMu8leg").get(p_hlt_lep1)[0];
        eff_lep2_leg1 = m_hlt_efficiencies.at("IsoMu17leg").get(p_hlt_lep2)[0];
        eff_lep2_leg2 = m_hlt_efficiencies.at("IsoMu8orIsoTkMu8leg").get(p_hlt_lep2)[0];
        //eff_lep2_tkleg2 = m_hlt_efficiencies.at("DoubleIsoMu17Mu8_TkMu8leg").get(p_hlt_lep2)[0];
    }
    else if (lep1.isMu && lep2.isEl) {
        eff_lep1_leg1 = m_hlt_efficiencies.at("IsoMu23leg").get(p_hlt_lep1)[0];
        eff_lep1_leg2 = m_hlt_efficiencies.at("IsoMu8leg").get(p_hlt_lep1)[0];
        eff_lep2_leg1 = m_hlt_efficiencies.at("EleMuHighPtleg").get(p_hlt_lep2)[0];
        eff_lep2_leg2 = m_hlt_efficiencies.at("MuEleLowPtleg").get(p_hlt_lep2)[0];
    }
    else if (lep1.isEl && lep2.isMu) {
        eff_lep1_leg1 = m_hlt_efficiencies.at("EleMuHighPtleg").get(p_hlt_lep1)[0];
        eff_lep1_leg2 = m_hlt_efficiencies.at("MuEleLowPtleg").get(p_hlt_lep1)[0];
        eff_lep2_leg1 = m_hlt_efficiencies.at("IsoMu23leg").get(p_hlt_lep2)[0];
        eff_lep2_leg2 = m_hlt_efficiencies.at("IsoMu8leg").get(p_hlt_lep2)[0];
    }
    else if (lep1.isEl && lep2.isEl){
        eff_lep1_leg1 = m_hlt_efficiencies.at("DoubleEleHighPtleg").get(p_hlt_lep1)[0];
        eff_lep1_leg2 = m_hlt_efficiencies.at("DoubleEleLowPtleg").get(p_hlt_lep1)[0];
        eff_lep2_leg1 = m_hlt_efficiencies.at("DoubleEleHighPtleg").get(p_hlt_lep2)[0];
        eff_lep2_leg2 = m_hlt_efficiencies.at("DoubleEleLowPtleg").get(p_hlt_lep2)[0];
    }
    else 
        std::cout << "We have something else then el or mu !!" << std::endl;

    float error_eff_lep1_leg1_up = 0.;
    float error_eff_lep1_leg2_up = 0.;
    float error_eff_lep1_tkleg2_up = 0.;
    float error_eff_lep2_leg1_up = 0.;
    float error_eff_lep2_leg2_up = 0.;
    float error_eff_lep2_tkleg2_up = 0.;

    if (lep1.isMu && lep2.isMu) {
        error_eff_lep1_leg1_up = m_hlt_efficiencies.at("IsoMu17leg").get(p_hlt_lep1)[2];
        error_eff_lep1_leg2_up = m_hlt_efficiencies.at("IsoMu8orIsoTkMu8leg").get(p_hlt_lep1)[2];
        //error_eff_lep1_tkleg2_up = m_hlt_efficiencies.at("DoubleIsoMu17Mu8_TkMu8leg").get(p_hlt_lep1)[2];
        error_eff_lep2_leg1_up = m_hlt_efficiencies.at("IsoMu17leg").get(p_hlt_lep2)[2];
        error_eff_lep2_leg2_up = m_hlt_efficiencies.at("IsoMu8orIsoTkMu8leg").get(p_hlt_lep2)[2];
        //error_eff_lep2_tkleg2_up = m_hlt_efficiencies.at("DoubleIsoMu17Mu8_TkMu8leg").get(p_hlt_lep2)[2];
    }
    else if (lep1.isMu && lep2.isEl) {
        error_eff_lep1_leg1_up = m_hlt_efficiencies.at("IsoMu23leg").get(p_hlt_lep1)[2];
        error_eff_lep1_leg2_up = m_hlt_efficiencies.at("IsoMu8leg").get(p_hlt_lep1)[2];
        error_eff_lep2_leg1_up = m_hlt_efficiencies.at("EleMuHighPtleg").get(p_hlt_lep2)[2];
        error_eff_lep2_leg2_up = m_hlt_efficiencies.at("MuEleLowPtleg").get(p_hlt_lep2)[2];
    }
    else if (lep1.isEl && lep2.isMu) {
        error_eff_lep1_leg1_up = m_hlt_efficiencies.at("EleMuHighPtleg").get(p_hlt_lep1)[2];
        error_eff_lep1_leg2_up = m_hlt_efficiencies.at("MuEleLowPtleg").get(p_hlt_lep1)[2];
        error_eff_lep2_leg1_up = m_hlt_efficiencies.at("IsoMu23leg").get(p_hlt_lep2)[2];
        error_eff_lep2_leg2_up = m_hlt_efficiencies.at("IsoMu8leg").get(p_hlt_lep2)[2];
    }
    else if (lep1.isEl && lep2.isEl){
        error_eff_lep1_leg1_up = m_hlt_efficiencies.at("DoubleEleHighPtleg").get(p_hlt_lep1)[2];
        error_eff_lep1_leg2_up = m_hlt_efficiencies.at("DoubleEleLowPtleg").get(p_hlt_lep1)[2];
        error_eff_lep2_leg1_up = m_hlt_efficiencies.at("DoubleEleHighPtleg").get(p_hlt_lep2)[2];
        error_eff_lep2_leg2_up = m_hlt_efficiencies.at("DoubleEleLowPtleg").get(p_hlt_lep2)[2];
    }

    float error_eff_lep1_leg1_down = 0.;
    float error_eff_lep1_leg2_down = 0.;
    float error_eff_lep1_tkleg2_down = 0.;
    float error_eff_lep2_leg1_down = 0.;
    float error_eff_lep2_leg2_down = 0.;
    float error_eff_lep2_tkleg2_down = 0.;

    if (lep1.isMu && lep2.isMu) {
        error_eff_lep1_leg1_down = m_hlt_efficiencies.at("IsoMu17leg").get(p_hlt_lep1)[1];
        error_eff_lep1_leg2_down = m_hlt_efficiencies.at("IsoMu8orIsoTkMu8leg").get(p_hlt_lep1)[1];
        //error_eff_lep1_tkleg2_down = m_hlt_efficiencies.at("DoubleIsoMu17Mu8_TkMu8leg").get(p_hlt_lep1)[1];
        error_eff_lep2_leg1_down = m_hlt_efficiencies.at("IsoMu17leg").get(p_hlt_lep2)[1];
        error_eff_lep2_leg2_down = m_hlt_efficiencies.at("IsoMu8orIsoTkMu8leg").get(p_hlt_lep2)[1];
        //error_eff_lep2_tkleg2_down = m_hlt_efficiencies.at("DoubleIsoMu17Mu8_TkMu8leg").get(p_hlt_lep2)[1];
    }
    else if (lep1.isMu && lep2.isEl) {
        error_eff_lep1_leg1_down = m_hlt_efficiencies.at("IsoMu23leg").get(p_hlt_lep1)[1];
        error_eff_lep1_leg2_down = m_hlt_efficiencies.at("IsoMu8leg").get(p_hlt_lep1)[1];
        error_eff_lep2_leg1_down = m_hlt_efficiencies.at("EleMuHighPtleg").get(p_hlt_lep2)[1];
        error_eff_lep2_leg2_down = m_hlt_efficiencies.at("MuEleLowPtleg").get(p_hlt_lep2)[1];
    }
    else if (lep1.isEl && lep2.isMu) {
        error_eff_lep1_leg1_down = m_hlt_efficiencies.at("EleMuHighPtleg").get(p_hlt_lep1)[1];
        error_eff_lep1_leg2_down = m_hlt_efficiencies.at("MuEleLowPtleg").get(p_hlt_lep1)[1];
        error_eff_lep2_leg1_down = m_hlt_efficiencies.at("IsoMu23leg").get(p_hlt_lep2)[1];
        error_eff_lep2_leg2_down = m_hlt_efficiencies.at("IsoMu8leg").get(p_hlt_lep2)[1];
    }
    else if (lep1.isEl && lep2.isEl){
        error_eff_lep1_leg1_down = m_hlt_efficiencies.at("DoubleEleHighPtleg").get(p_hlt_lep1)[1];
        error_eff_lep1_leg2_down = m_hlt_efficiencies.at("DoubleEleLowPtleg").get(p_hlt_lep1)[1];
        error_eff_lep2_leg1_down = m_hlt_efficiencies.at("DoubleEleHighPtleg").get(p_hlt_lep2)[1];
        error_eff_lep2_leg2_down = m_hlt_efficiencies.at("DoubleEleLowPtleg").get(p_hlt_lep2)[1];
    }


    float nominal = -(eff_lep1_leg1 * eff_lep2_leg1) +
        (1 - (1 - eff_lep1_leg2) * (1 - eff_lep1_tkleg2)) * eff_lep2_leg1 +
        eff_lep1_leg1 * (1 - (1 - eff_lep2_leg2) * (1 - eff_lep2_tkleg2));

    float error_squared_up =
        std::pow(1 - eff_lep2_leg1 - (1 - eff_lep2_leg2) * (1 - eff_lep2_tkleg2), 2) *
        std::pow(error_eff_lep1_leg1_up, 2) +
        std::pow(1 - eff_lep1_tkleg2, 2) * std::pow(eff_lep2_leg1, 2) *
        std::pow(error_eff_lep1_leg2_up, 2) +
        std::pow(1 - eff_lep1_leg2, 2) * std::pow(eff_lep2_leg1, 2) *
        std::pow(error_eff_lep1_tkleg2_up, 2) +
        std::pow(1 - eff_lep1_leg1 - (1 - eff_lep1_leg2) * (1 - eff_lep1_tkleg2), 2) *
        std::pow(error_eff_lep2_leg1_up, 2) +
        std::pow(eff_lep1_leg1, 2) * std::pow(1 - eff_lep2_tkleg2, 2) *
        std::pow(error_eff_lep2_leg2_up, 2) +
        std::pow(eff_lep1_leg1, 2) * std::pow(1 - eff_lep2_leg2, 2) *
        std::pow(error_eff_lep2_tkleg2_up, 2);

    float error_squared_down = 
        std::pow(1 - eff_lep2_leg1 - (1 - eff_lep2_leg2) * (1 - eff_lep2_tkleg2), 2) *
        std::pow(error_eff_lep1_leg1_down, 2) +
        std::pow(1 - eff_lep1_tkleg2, 2) * std::pow(eff_lep2_leg1, 2) *
        std::pow(error_eff_lep1_leg2_down, 2) +
        std::pow(1 - eff_lep1_leg2, 2) * std::pow(eff_lep2_leg1, 2) *
        std::pow(error_eff_lep1_tkleg2_down, 2) +
        std::pow(1 - eff_lep1_leg1 - (1 - eff_lep1_leg2) * (1 - eff_lep1_tkleg2), 2) *
        std::pow(error_eff_lep2_leg1_down, 2) +
        std::pow(eff_lep1_leg1, 2) * std::pow(1 - eff_lep2_tkleg2, 2) *
        std::pow(error_eff_lep2_leg2_down, 2) +
        std::pow(eff_lep1_leg1, 2) * std::pow(1 - eff_lep2_leg2, 2) *
        std::pow(error_eff_lep2_tkleg2_down, 2);

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


