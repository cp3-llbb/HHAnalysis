#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>
#include <cp3_llbb/Framework/interface/BTagsAnalyzer.h>
#include <cp3_llbb/HHAnalysis/interface/Categories.h>

#include <cp3_llbb/Framework/interface/GenParticlesProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/LeptonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/METProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cmath>

#define HHANADEBUG 0

void HHAnalyzer::registerCategories(CategoryManager& manager, const edm::ParameterSet& config) {
    manager.new_category<MuMuCategory>("mumu", "Category with leading leptons as two muons", config);
    manager.new_category<ElElCategory>("elel", "Category with leading leptons as two electrons", config);
    manager.new_category<ElMuCategory>("elmu", "Category with leading leptons as electron, subleading as muon", config);
    manager.new_category<MuElCategory>("muel", "Category with leading leptons as muon, subleading as electron", config);
}


void HHAnalyzer::analyze(const edm::Event& event, const edm::EventSetup&, const ProducersManager& producers, const AnalyzersManager&, const CategoryManager&) {

    float mh = event.isRealData() ? 125.02 : 125.0;
    LorentzVector null_p4(0., 0., 0., 0.);

    // ***** ***** *****
    // Trigger Matching
    // ***** ***** *****

    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    //Function that tries to match `lepton` with an online object, using a deltaR and a deltaPt cut   
    //Returns the index inside the HLTProducer collection, or -1 if no match is found.
    //(Taken from https://github.com/blinkseb/TTAnalysis/blob/c2a2d5de3e4281943c19c582afb452b8ef6457f1/plugins/TTAnalyzer.cc#L533)
    auto matchOfflineLepton = [&](HH::Lepton& lepton) {

        if (lepton.hlt_already_tried_matching)
            return lepton.hlt_idx;
        float min_dr = std::numeric_limits<float>::max();
        float final_dpt_over_pt = std::numeric_limits<float>::max();

        int8_t index = -1;
        for (size_t hlt_object = 0; hlt_object < hlt.object_p4.size(); hlt_object++) {

            float dr = ROOT::Math::VectorUtil::DeltaR(lepton.p4, hlt.object_p4[hlt_object]);
            float dpt_over_pt = fabs(lepton.p4.Pt() - hlt.object_p4[hlt_object].Pt()) / lepton.p4.Pt();

            if (dr > m_hltDRCut)
                continue;

            if (dpt_over_pt > m_hltDPtCut)
                continue;

            if (dr < min_dr) {
                min_dr = dr;
                final_dpt_over_pt = dpt_over_pt;
                index = hlt_object;
            }
        }
        lepton.hlt_idx = index;
        lepton.hlt_already_tried_matching = true;
        lepton.hlt_DR_matchedObject = min_dr;
        lepton.hlt_DPtOverPt_matchedObject = final_dpt_over_pt;
        return index;
    };


    // ********** 
    // Leptons and dileptons
    // ********** 
    const ElectronsProducer& allelectrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& allmuons = producers.get<MuonsProducer>("muons");

    leptons.clear();
    ll.clear();

    // Fill lepton structures
    for (unsigned int ielectron = 0; ielectron < allelectrons.p4.size(); ielectron++)
    {
        if (allelectrons.p4[ielectron].Pt() > m_subleadingElectronPtCut
            && fabs(allelectrons.p4[ielectron].Eta()) < m_electronEtaCut) 
        {
            electrons.push_back(ielectron);
            HH::Lepton ele;
            ele.p4 = allelectrons.p4[ielectron];
            ele.charge = allelectrons.charge[ielectron];
            ele.idx = ielectron;
            ele.isMu = false;
            ele.isEl = true;
            ele.id_L = allelectrons.ids[ielectron][m_electron_loose_wp_name];
            ele.id_T = allelectrons.ids[ielectron][m_electron_tight_wp_name];
            ele.iso_L = allelectrons.isEB[ielectron] ? (allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut_EB_Loose) : (allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut_EE_Loose);
            ele.iso_T = allelectrons.isEB[ielectron] ? (allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut_EB_Tight) : (allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut_EE_Tight);
            ele.gen_matched = allelectrons.matched[ielectron];
            ele.gen_p4 = ele.gen_matched ? allelectrons.gen_p4[ielectron] : null_p4;
            ele.gen_DR = ele.gen_matched ? ROOT::Math::VectorUtil::DeltaR(ele.p4, ele.gen_p4): -1.;
            ele.gen_DPtOverPt = ele.gen_matched ? (ele.p4.Pt() - ele.gen_p4.Pt()) / ele.p4.Pt() : -10.;
            leptons.push_back(ele);
        }
    }//end of loop on electrons

    for (unsigned int imuon = 0; imuon < allmuons.p4.size(); imuon++)
    {
        if (allmuons.p4[imuon].Pt() > m_subleadingMuonPtCut
            && fabs(allmuons.p4[imuon].Eta()) < m_muonEtaCut)
        {
            muons.push_back(imuon);
            HH::Lepton mu;
            mu.p4 = allmuons.p4[imuon];
            mu.charge = allmuons.charge[imuon];
            mu.idx = imuon;
            mu.isMu = true;
            mu.isEl = false;
            mu.id_L = allmuons.isLoose[imuon];
            mu.id_T = allmuons.isTight[imuon];
            mu.iso_L = allmuons.relativeIsoR04_deltaBeta[imuon] < m_muonLooseIsoCut;
            mu.iso_T = allmuons.relativeIsoR04_deltaBeta[imuon] < m_muonTightIsoCut;
            mu.gen_matched = allmuons.matched[imuon];
            mu.gen_p4 = mu.gen_matched ? allmuons.gen_p4[imuon] : null_p4;
            mu.gen_DR = mu.gen_matched ? ROOT::Math::VectorUtil::DeltaR(mu.p4, mu.gen_p4) : -1.;
            mu.gen_DPtOverPt = mu.gen_matched ? (mu.p4.Pt() - mu.gen_p4.Pt()) / mu.p4.Pt() : -10.;
            leptons.push_back(mu);
        }
    }//end of loop on muons

    // sort leptons by pt (ignoring flavour, id and iso)
    std::sort(leptons.begin(), leptons.end(), [](const HH::Lepton& lep1, const HH::Lepton& lep2) { return lep1.p4.Pt() > lep2.p4.Pt(); });     

    // Fill lepton maps
    for (unsigned int i_id_iso = 0; i_id_iso < map_l.size(); i_id_iso++)
        map_l[i_id_iso].clear();

    for (unsigned int ilepton = 0; ilepton < leptons.size(); ilepton++)
    {
        if (leptons[ilepton].id_L)
        {
            map_l[lepID::L * lepIso::Count + lepIso::no].push_back(ilepton);
            if (leptons[ilepton].iso_L)
                map_l[lepID::L * lepIso::Count + lepIso::L].push_back(ilepton);
            if (leptons[ilepton].iso_T)
                map_l[lepID::L * lepIso::Count + lepIso::T].push_back(ilepton);
        }
        if (leptons[ilepton].id_T)
        {
            map_l[lepID::T * lepIso::Count + lepIso::no].push_back(ilepton);
            if (leptons[ilepton].iso_L)
                map_l[lepID::T * lepIso::Count + lepIso::L].push_back(ilepton);
            if (leptons[ilepton].iso_T)
                map_l[lepID::T * lepIso::Count + lepIso::T].push_back(ilepton);
        }
    }
           
    for (unsigned int ilep1 = 0; ilep1 < leptons.size(); ilep1++)
    {
        if ((leptons[ilep1].isMu && leptons[ilep1].p4.Pt() < m_leadingMuonPtCut) || (leptons[ilep1].isEl && leptons[ilep1].p4.Pt() < m_leadingElectronPtCut)) continue;

        for (unsigned int ilep2 = ilep1+1; ilep2 < leptons.size(); ilep2++)
        {
            HH::Dilepton dilep;
            dilep.p4 = leptons[ilep1].p4 + leptons[ilep2].p4;
            dilep.idxs = std::make_pair(leptons[ilep1].idx, leptons[ilep2].idx);
            dilep.ilep1 = ilep1;
            dilep.ilep2 = ilep2;
            dilep.isOS = leptons[ilep1].charge * leptons[ilep2].charge < 0;
            dilep.isMuMu = leptons[ilep1].isMu && leptons[ilep2].isMu;
            dilep.isElEl = leptons[ilep1].isEl && leptons[ilep2].isEl;
            dilep.isElMu = leptons[ilep1].isEl && leptons[ilep2].isMu;
            dilep.isMuEl = leptons[ilep1].isMu && leptons[ilep2].isEl;
            dilep.isSF = dilep.isMuMu || dilep.isElEl;
            dilep.id_LL = leptons[ilep1].id_L && leptons[ilep2].id_L;
            dilep.id_LT = leptons[ilep1].id_L && leptons[ilep2].id_T;
            dilep.id_TL = leptons[ilep1].id_T && leptons[ilep2].id_L;
            dilep.id_TT = leptons[ilep1].id_T && leptons[ilep2].id_T;
            dilep.iso_LL = leptons[ilep1].iso_L && leptons[ilep2].iso_L;
            dilep.iso_LT = leptons[ilep1].iso_L && leptons[ilep2].iso_T;
            dilep.iso_TL = leptons[ilep1].iso_T && leptons[ilep2].iso_L;
            dilep.iso_TT = leptons[ilep1].iso_T && leptons[ilep2].iso_T;
            dilep.DR_l_l = ROOT::Math::VectorUtil::DeltaR(leptons[ilep1].p4, leptons[ilep2].p4);
            dilep.DPhi_l_l = fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ilep1].p4, leptons[ilep2].p4));
            if (!hlt.paths.empty()) dilep.hlt_idxs = std::make_pair(matchOfflineLepton(leptons[ilep1]),matchOfflineLepton(leptons[ilep2]));
            dilep.gen_matched = leptons[ilep1].gen_matched && leptons[ilep2].gen_matched;
            dilep.gen_p4 = dilep.gen_matched ? leptons[ilep1].gen_p4 + leptons[ilep2].gen_p4 : null_p4;
            dilep.gen_DR = dilep.gen_matched ? ROOT::Math::VectorUtil::DeltaR(dilep.p4, dilep.gen_p4) : -1.;
            dilep.gen_DPtOverPt = dilep.gen_matched ? (dilep.p4.Pt() - dilep.gen_p4.Pt()) / dilep.p4.Pt() : -10.;
            ll.push_back(dilep); 
        }
    }

    // Fill dilepton maps
    for (unsigned int i = 0; i < map_ll.size(); i++)
        map_ll[i].clear();

    // map is stored as lep1_ID, lep1_Iso, lep2_ID, lep2_Iso
    int bitA = lepIso::Count;
    int bitB = bitA * lepID::Count;
    int bitC = bitB * lepIso::Count;
    for (unsigned int ill = 0; ill < ll.size(); ill++)
    {
        if (ll[ill].id_LL)
        {
            map_ll[lepID::L * bitC + lepIso::no * bitB + lepID::L * bitA + lepIso::no].push_back(ill);
            if (leptons[ll[ill].ilep2].iso_L)
                map_ll[lepID::L * bitC + lepIso::no * bitB + lepID::L * bitA + lepIso::L].push_back(ill);
            if (leptons[ll[ill].ilep2].iso_L)
                map_ll[lepID::L * bitC + lepIso::no * bitB + lepID::L * bitA + lepIso::T].push_back(ill);
            if (leptons[ll[ill].ilep1].iso_L)
            {
                map_ll[lepID::L * bitC + lepIso::L * bitB + lepID::L * bitA + lepIso::no].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::L * bitC + lepIso::L * bitB + lepID::L * bitA + lepIso::L].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::L * bitC + lepIso::L * bitB + lepID::L * bitA + lepIso::T].push_back(ill);
            }
            if (leptons[ll[ill].ilep1].iso_T)
            {
                map_ll[lepID::L * bitC + lepIso::T * bitB + lepID::L * bitA + lepIso::no].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::L * bitC + lepIso::T * bitB + lepID::L * bitA + lepIso::L].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::L * bitC + lepIso::T * bitB + lepID::L * bitA + lepIso::T].push_back(ill);
            }
        }
        if (ll[ill].id_LT)
        {
            map_ll[lepID::L * bitC + lepIso::no * bitB + lepID::T * bitA + lepIso::no].push_back(ill);
            if (leptons[ll[ill].ilep2].iso_L)
                map_ll[lepID::L * bitC + lepIso::no * bitB + lepID::T * bitA + lepIso::L].push_back(ill);
            if (leptons[ll[ill].ilep2].iso_L)
                map_ll[lepID::L * bitC + lepIso::no * bitB + lepID::T * bitA + lepIso::T].push_back(ill);
            if (leptons[ll[ill].ilep1].iso_L)
            {
                map_ll[lepID::L * bitC + lepIso::L * bitB + lepID::T * bitA + lepIso::no].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::L * bitC + lepIso::L * bitB + lepID::T * bitA + lepIso::L].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::L * bitC + lepIso::L * bitB + lepID::T * bitA + lepIso::T].push_back(ill);
            }
            if (leptons[ll[ill].ilep1].iso_T)
            {
                map_ll[lepID::L * bitC + lepIso::T * bitB + lepID::T * bitA + lepIso::no].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::L * bitC + lepIso::T * bitB + lepID::T * bitA + lepIso::L].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::L * bitC + lepIso::T * bitB + lepID::T * bitA + lepIso::T].push_back(ill);
            }
        }
        if (ll[ill].id_TL)
        {
            map_ll[lepID::T * bitC + lepIso::no * bitB + lepID::L * bitA + lepIso::no].push_back(ill);
            if (leptons[ll[ill].ilep2].iso_L)
                map_ll[lepID::T * bitC + lepIso::no * bitB + lepID::L * bitA + lepIso::L].push_back(ill);
            if (leptons[ll[ill].ilep2].iso_L)
                map_ll[lepID::T * bitC + lepIso::no * bitB + lepID::L * bitA + lepIso::T].push_back(ill);
            if (leptons[ll[ill].ilep1].iso_L)
            {
                map_ll[lepID::T * bitC + lepIso::L * bitB + lepID::L * bitA + lepIso::no].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::T * bitC + lepIso::L * bitB + lepID::L * bitA + lepIso::L].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::T * bitC + lepIso::L * bitB + lepID::L * bitA + lepIso::T].push_back(ill);
            }
            if (leptons[ll[ill].ilep1].iso_T)
            {
                map_ll[lepID::T * bitC + lepIso::T * bitB + lepID::L * bitA + lepIso::no].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::T * bitC + lepIso::T * bitB + lepID::L * bitA + lepIso::L].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::T * bitC + lepIso::T * bitB + lepID::L * bitA + lepIso::T].push_back(ill);
            }
        }
        if (ll[ill].id_TT)
        {
            map_ll[lepID::T * bitC + lepIso::no * bitB + lepID::T * bitA + lepIso::no].push_back(ill);
            if (leptons[ll[ill].ilep2].iso_L)
                map_ll[lepID::T * bitC + lepIso::no * bitB + lepID::T * bitA + lepIso::L].push_back(ill);
            if (leptons[ll[ill].ilep2].iso_L)
                map_ll[lepID::T * bitC + lepIso::no * bitB + lepID::T * bitA + lepIso::T].push_back(ill);
            if (leptons[ll[ill].ilep1].iso_L)
            {
                map_ll[lepID::T * bitC + lepIso::L * bitB + lepID::T * bitA + lepIso::no].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::T * bitC + lepIso::L * bitB + lepID::T * bitA + lepIso::L].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::T * bitC + lepIso::L * bitB + lepID::T * bitA + lepIso::T].push_back(ill);
            }
            if (leptons[ll[ill].ilep1].iso_T)
            {
                map_ll[lepID::T * bitC + lepIso::T * bitB + lepID::T * bitA + lepIso::no].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::T * bitC + lepIso::T * bitB + lepID::T * bitA + lepIso::L].push_back(ill);
                if (leptons[ll[ill].ilep2].iso_L)
                    map_ll[lepID::T * bitC + lepIso::T * bitB + lepID::T * bitA + lepIso::T].push_back(ill);
            }
        }
    }

    // ***** 
    // Adding MET(s)
    // ***** 
    const METProducer& pf_met = producers.get<METProducer>("met");
    HH::Met mymet;
    mymet.p4 = pf_met.p4;
    mymet.isNoHF = false;
    mymet.gen_matched = false;
    mymet.gen_p4 = null_p4;
    mymet.gen_DR = -1.;
    mymet.gen_DPhi = -1.;
    mymet.gen_DPtOverPt = -10.;
    if (!event.isRealData())
    { // genMet is not constructed in the framework, so construct it manually out of the neutrinos hanging around the mc particles
        const GenParticlesProducer& gp = producers.get<GenParticlesProducer>("gen_particles");
        for (unsigned int ip = 0; ip < gp.pruned_p4.size(); ip++) {
            std::bitset<15> flags (gp.pruned_status_flags[ip]);
            if (!flags.test(13)) continue; // take the last copies
            if (abs(gp.pruned_pdg_id[ip]) == 12 || abs(gp.pruned_pdg_id[ip]) == 14 || abs(gp.pruned_pdg_id[ip]) == 16)
            {
                mymet.gen_matched = true;
                mymet.gen_p4 += gp.pruned_p4[ip];
            }
        }
        mymet.gen_DR = mymet.gen_matched ? ROOT::Math::VectorUtil::DeltaR(mymet.p4, mymet.gen_p4) : -1.;
        mymet.gen_DPhi = mymet.gen_matched ? fabs(ROOT::Math::VectorUtil::DeltaPhi(mymet.p4, mymet.gen_p4)) : -1.;
        mymet.gen_DPtOverPt = mymet.gen_matched ? (mymet.p4.Pt() - mymet.gen_p4.Pt()) / mymet.p4.Pt() : -10.;
    }
    met.push_back(mymet);
    //const METProducer& puppi_met = producers.get<METProducer>("puppimet");
    // TODO: adding puppi met will require changing the Met AND DileptonMet struct

    for (unsigned int imet = 0; imet < met.size(); imet++)
    {
        for (unsigned int ill = 0; ill < ll.size(); ill++)
        {
            HH::DileptonMet myllmet;
// DileptonMet inherits from Dilepton struct, initalize everything properly
// FIXME: there is very probably a cleaner way to do
            myllmet.p4 = ll[ill].p4 + met[imet].p4;
            // blind copy of the ll content
            myllmet.idxs = std::make_pair(ll[ill].idxs.first, ll[ill].idxs.second);
            myllmet.ilep1 = ll[ill].ilep1;
            myllmet.ilep2 = ll[ill].ilep2;
            myllmet.isOS = ll[ill].isOS;
            myllmet.isMuMu = ll[ill].isMuMu;
            myllmet.isElEl = ll[ill].isElEl;
            myllmet.isElMu = ll[ill].isElMu;
            myllmet.isMuEl = ll[ill].isMuEl;
            myllmet.isSF = ll[ill].isSF;
            myllmet.id_LL = ll[ill].id_LL;
            myllmet.id_LT = ll[ill].id_LT;
            myllmet.id_TL = ll[ill].id_TL;
            myllmet.id_TT = ll[ill].id_TT;
            myllmet.iso_LL = ll[ill].iso_LL;
            myllmet.iso_LT = ll[ill].iso_LT;
            myllmet.iso_TL = ll[ill].iso_TL;
            myllmet.iso_TT = ll[ill].iso_TT;
            myllmet.DR_l_l = ll[ill].DR_l_l;
            myllmet.DPhi_l_l = ll[ill].DPhi_l_l;
            // content specific to HH:DileptonMet
            myllmet.ill = ill;
            myllmet.imet = imet;
            myllmet.isNoHF = met[imet].isNoHF;
            float dphi = fabs(ROOT::Math::VectorUtil::DeltaPhi(ll[ill].p4, met[imet].p4));
            myllmet.DPhi_ll_met = dphi;
            float mindphi = std::min(fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ll[ill].idxs.first].p4, met[imet].p4)), fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ll[ill].idxs.second].p4, met[imet].p4)));
            myllmet.minDPhi_l_met = mindphi; 
            float maxdphi = std::max(fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ll[ill].idxs.first].p4, met[imet].p4)), fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ll[ill].idxs.second].p4, met[imet].p4)));
            myllmet.maxDPhi_l_met = maxdphi;
            myllmet.MT = (ll[ill].p4 + met[imet].p4).M();
            myllmet.MT_formula = std::sqrt(2 * ll[ill].p4.Pt() * met[imet].p4.Pt() * (1-std::cos(dphi)));
            myllmet.projectedMet = mindphi >= M_PI ? met[imet].p4.Pt() : met[imet].p4.Pt() * std::sin(mindphi);
            myllmet.gen_matched = ll[ill].gen_matched && met[imet].gen_matched;
            myllmet.gen_p4 = myllmet.gen_matched ? ll[ill].gen_p4 + met[imet].gen_p4 : null_p4;
            myllmet.gen_DR = myllmet.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myllmet.p4, myllmet.gen_p4) : -1.;
            myllmet.gen_DPhi = myllmet.gen_matched ? fabs(ROOT::Math::VectorUtil::DeltaPhi(myllmet.p4, myllmet.gen_p4)) : -1.;
            myllmet.gen_DPtOverPt = myllmet.gen_matched ? (myllmet.p4.Pt() - myllmet.gen_p4.Pt()) / myllmet.p4.Pt() : -10.;
            llmet.push_back(myllmet);
        }
    }

    // Fill dilepton+met maps
    // FIXME: for now only store pf_met
    // so ll and llmet structures are in sync
    for (unsigned int i = 0; i < map_llmet.size(); i++)
    {
        map_llmet[i].clear();
        for (unsigned int j = 0; j < map_ll[i].size(); j++)
        {
            map_llmet[i].push_back(map_ll[i][j]);
        }
    }

    // ***** 
    // Jets and dijets 
    // ***** 
    const JetsProducer& alljets = producers.get<JetsProducer>("jets");
    jets.clear();
    for (unsigned int imap_j = 0; imap_j < map_j.size(); imap_j++)
        map_j[imap_j].clear();

    bitA = btagWP::Count;
    // (The following include filling up the jets map
    int count = 0;
    for (unsigned int ijet = 0; ijet < alljets.p4.size(); ijet++)
    {
        if ((alljets.p4[ijet].Pt() > m_jetPtCut) 
            && (fabs(alljets.p4[ijet].Eta()) < m_jetEtaCut))
        {
            HH::Jet myjet;
            myjet.p4 = alljets.p4[ijet];
            myjet.idx = ijet;
            myjet.id_L = alljets.passLooseID[ijet];
            myjet.id_T = alljets.passTightID[ijet];
            myjet.id_TLV = alljets.passTightLeptonVetoID[ijet];
            myjet.CSV = alljets.getBTagDiscriminant(ijet, "pfCombinedInclusiveSecondaryVertexV2BJetTags");
            myjet.JP = alljets.getBTagDiscriminant(ijet, "pfJetProbabilityBJetTags");
            float mybtag = alljets.getBTagDiscriminant(ijet, m_jet_bDiscrName);
            myjet.btag_L = mybtag > m_jet_bDiscrCut_loose;
            myjet.btag_M = mybtag > m_jet_bDiscrCut_medium;
            myjet.btag_T = mybtag > m_jet_bDiscrCut_tight;
            myjet.gen_matched_bParton = (alljets.partonFlavor[ijet]) == 5;
            myjet.gen_matched_bHadron = (alljets.hadronFlavor[ijet]) == 5;
            myjet.gen_matched = alljets.matched[ijet];
            myjet.gen_p4 = myjet.gen_matched ? alljets.gen_p4[ijet] : null_p4;
            myjet.gen_DR = myjet.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myjet.p4, myjet.gen_p4) : -1.;
            myjet.gen_DPtOverPt = myjet.gen_matched ? (myjet.p4.Pt() - myjet.gen_p4.Pt()) / myjet.p4.Pt() : -10.;
            myjet.gen_b = (alljets.hadronFlavor[ijet]) == 5; // redundant with gen_matched_bHadron defined above
            myjet.gen_c = (alljets.hadronFlavor[ijet]) == 4;
            myjet.gen_l = (alljets.hadronFlavor[ijet]) < 4;
            jets.push_back(myjet);
            // filling maps
            // no jet ID
            map_j[jetID::no * bitA + btagWP::no].push_back(count);
            if (myjet.btag_L)
                map_j[jetID::no * bitA + btagWP::L].push_back(count);
            if (myjet.btag_M)
                map_j[jetID::no * bitA + btagWP::M].push_back(count);
            if (myjet.btag_T)
                map_j[jetID::no * bitA + btagWP::T].push_back(count);
            // loose jet ID
            if (myjet.id_L)
            {
                map_j[jetID::L * bitA + btagWP::no].push_back(count);
                if (myjet.btag_L)
                    map_j[jetID::L * bitA + btagWP::L].push_back(count);
                if (myjet.btag_M)
                    map_j[jetID::L * bitA + btagWP::M].push_back(count);
                if (myjet.btag_T)
                    map_j[jetID::L * bitA + btagWP::T].push_back(count);
            }
            // tight jet ID
            if (myjet.id_T)
            {
                map_j[jetID::T * bitA + btagWP::no].push_back(count);
                if (myjet.btag_L)
                    map_j[jetID::T * bitA + btagWP::L].push_back(count);
                if (myjet.btag_M)
                    map_j[jetID::T * bitA + btagWP::M].push_back(count);
                if (myjet.btag_T)
                    map_j[jetID::T * bitA + btagWP::T].push_back(count);
            }
            // tight jet ID lepton veto
            if (myjet.id_TLV)
            {
                map_j[jetID::TLV * bitA + btagWP::no].push_back(count);
                if (myjet.btag_L)
                    map_j[jetID::TLV * bitA + btagWP::L].push_back(count);
                if (myjet.btag_M)
                    map_j[jetID::TLV * bitA + btagWP::M].push_back(count);
                if (myjet.btag_T)
                    map_j[jetID::TLV * bitA + btagWP::T].push_back(count);
            }
            count++;
        }
    }

    jj.clear();
    // Do NOT change the loop logic here: we expect [0] to be made out of the leading jets
    for (unsigned int ijj = 0; ijj < map_jj.size(); ijj++)
        map_jj[ijj].clear();
    bitA = jetPair::Count;
    bitB = bitA * btagWP::Count;
    bitC = bitB * btagWP::Count;
    int bitD = bitC * jetID::Count;
    count = 0;

    for (unsigned int ijet1 = 0; ijet1 < jets.size(); ijet1++)
    {
        for (unsigned int ijet2 = ijet1 + 1; ijet2 < jets.size(); ijet2++)
        {
            HH::Dijet myjj;
            myjj.p4 = jets[ijet1].p4 + jets[ijet2].p4;
            myjj.idxs = std::make_pair(jets[ijet1].idx, jets[ijet2].idx);
            myjj.ijet1 = ijet1;
            myjj.ijet2 = ijet2;
            myjj.btag_LL = jets[ijet1].btag_L && jets[ijet2].btag_L;
            myjj.btag_LM = jets[ijet1].btag_L && jets[ijet2].btag_M;
            myjj.btag_LT = jets[ijet1].btag_L && jets[ijet2].btag_T;
            myjj.btag_ML = jets[ijet1].btag_M && jets[ijet2].btag_L;
            myjj.btag_MM = jets[ijet1].btag_M && jets[ijet2].btag_M;
            myjj.btag_MT = jets[ijet1].btag_M && jets[ijet2].btag_T;
            myjj.btag_TL = jets[ijet1].btag_T && jets[ijet2].btag_L;
            myjj.btag_TM = jets[ijet1].btag_T && jets[ijet2].btag_M;
            myjj.btag_TT = jets[ijet1].btag_T && jets[ijet2].btag_T;
            myjj.sumCSV = jets[ijet1].CSV + jets[ijet2].CSV;
            myjj.sumJP = jets[ijet1].JP + jets[ijet2].JP;
            myjj.DR_j_j = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, jets[ijet2].p4);
            myjj.DPhi_j_j = fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[ijet1].p4, jets[ijet2].p4));
            myjj.gen_matched_bbPartons = jets[ijet1].gen_matched_bParton && jets[ijet2].gen_matched_bParton; 
            myjj.gen_matched_bbHadrons = jets[ijet1].gen_matched_bHadron && jets[ijet2].gen_matched_bHadron; 
            myjj.gen_matched = jets[ijet1].gen_matched && jets[ijet2].gen_matched;
            myjj.gen_p4 = myjj.gen_matched ? jets[ijet1].gen_p4 + jets[ijet2].gen_p4 : null_p4;
            myjj.gen_DR = myjj.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myjj.p4, myjj.gen_p4) : -1.;
            myjj.gen_DPtOverPt = myjj.gen_matched ? (myjj.p4.Pt() - myjj.gen_p4.Pt()) / myjj.p4.Pt() : -10.;
            myjj.gen_bb = (jets[ijet1].gen_b && jets[ijet2].gen_b);
            myjj.gen_bc = (jets[ijet1].gen_b && jets[ijet2].gen_c) || (jets[ijet1].gen_c && jets[ijet2].gen_b);
            myjj.gen_bl = (jets[ijet1].gen_b && jets[ijet2].gen_l) || (jets[ijet1].gen_l && jets[ijet2].gen_b);
            myjj.gen_cc = (jets[ijet1].gen_c && jets[ijet2].gen_c);
            myjj.gen_cl = (jets[ijet1].gen_c && jets[ijet2].gen_l) || (jets[ijet1].gen_l && jets[ijet2].gen_c);
            myjj.gen_ll = (jets[ijet1].gen_l && jets[ijet2].gen_l);
            jj.push_back(myjj);
            // fill dijet map
            map_jj[jetID::no * bitD + jetID::no * bitC + btagWP::no * bitB + btagWP::no * bitA + jetPair::ht].push_back(count);
            if (jets[ijet1].btag_L)
                map_jj[jetID::no * bitD + jetID::no * bitC + btagWP::L * bitB + btagWP::no * bitA + jetPair::ht].push_back(count);
            if (jets[ijet2].btag_L)
                map_jj[jetID::no * bitD + jetID::no * bitC + btagWP::no * bitB + btagWP::L * bitA + jetPair::ht].push_back(count);
            if (myjj.btag_LL)
                map_jj[jetID::no * bitD + jetID::no * bitC + btagWP::L * bitB + btagWP::L * bitA + jetPair::ht].push_back(count);
            if (myjj.btag_LM)
                map_jj[jetID::no * bitD + jetID::no * bitC + btagWP::L * bitB + btagWP::M * bitA + jetPair::ht].push_back(count);
            if (myjj.btag_LT)
                map_jj[jetID::no * bitD + jetID::no * bitC + btagWP::L * bitB + btagWP::T * bitA + jetPair::ht].push_back(count);
            if (myjj.btag_ML)
                map_jj[jetID::no * bitD + jetID::no * bitC + btagWP::M * bitB + btagWP::L * bitA + jetPair::ht].push_back(count);
            if (myjj.btag_MM)
                map_jj[jetID::no * bitD + jetID::no * bitC + btagWP::M * bitB + btagWP::M * bitA + jetPair::ht].push_back(count);
            if (myjj.btag_MT)
                map_jj[jetID::no * bitD + jetID::no * bitC + btagWP::M * bitB + btagWP::T * bitA + jetPair::ht].push_back(count);
            if (myjj.btag_TL)
                map_jj[jetID::no * bitD + jetID::no * bitC + btagWP::T * bitB + btagWP::L * bitA + jetPair::ht].push_back(count);
            if (myjj.btag_TM)
                map_jj[jetID::no * bitD + jetID::no * bitC + btagWP::T * bitB + btagWP::M * bitA + jetPair::ht].push_back(count);
            if (myjj.btag_TT)
                map_jj[jetID::no * bitD + jetID::no * bitC + btagWP::T * bitB + btagWP::T * bitA + jetPair::ht].push_back(count);
            if (jets[ijet1].id_L)
            {
                for (int ibtag = 0; ibtag < btagWP::Count; ibtag++)
                    for (int jbtag = 0; jbtag < btagWP::Count; jbtag++)
                        if (map_jj[jetID::no * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].size() > 0)
                            if (map_jj[jetID::no * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].back() == count)
                                map_jj[jetID::L  * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].push_back(count);
                if (jets[ijet2].id_L)
                    for (int ibtag = 0; ibtag < btagWP::Count; ibtag++)
                        for (int jbtag = 0; jbtag < btagWP::Count; jbtag++)
                            if (map_jj[jetID::no * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].size() > 0)
                                if (map_jj[jetID::no * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].back() == count)
                                    map_jj[jetID::L  * bitD + jetID::L  * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].push_back(count);
                if (jets[ijet2].id_T)
                    for (int ibtag = 0; ibtag < btagWP::Count; ibtag++)
                        for (int jbtag = 0; jbtag < btagWP::Count; jbtag++)
                            if (map_jj[jetID::no * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].size() > 0)
                                if (map_jj[jetID::no * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].back() == count)
                                    map_jj[jetID::L  * bitD + jetID::T  * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].push_back(count);
                if (jets[ijet2].id_TLV)
                    for (int ibtag = 0; ibtag < btagWP::Count; ibtag++)
                        for (int jbtag = 0; jbtag < btagWP::Count; jbtag++)
                            if (map_jj[jetID::no * bitD + jetID::no  * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].size() > 0)
                                if (map_jj[jetID::no * bitD + jetID::no  * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].back() == count)
                                    map_jj[jetID::L  * bitD + jetID::TLV * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].push_back(count);
            }
            if (jets[ijet1].id_T)
            {
                for (int ibtag = 0; ibtag < btagWP::Count; ibtag++)
                    for (int jbtag = 0; jbtag < btagWP::Count; jbtag++)
                        if (map_jj[jetID::no * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].size() > 0)
                            if (map_jj[jetID::no * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].back() == count)
                                map_jj[jetID::T  * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].push_back(count);
                if (jets[ijet2].id_L)
                    for (int ibtag = 0; ibtag < btagWP::Count; ibtag++)
                        for (int jbtag = 0; jbtag < btagWP::Count; jbtag++)
                            if (map_jj[jetID::no * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].size() > 0)
                                if (map_jj[jetID::no * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].back() == count)
                                    map_jj[jetID::T  * bitD + jetID::L  * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].push_back(count);
                if (jets[ijet2].id_T)
                    for (int ibtag = 0; ibtag < btagWP::Count; ibtag++)
                        for (int jbtag = 0; jbtag < btagWP::Count; jbtag++)
                            if (map_jj[jetID::no * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].size() > 0)
                                if (map_jj[jetID::no * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].back() == count)
                                    map_jj[jetID::T  * bitD + jetID::T  * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].push_back(count);
                if (jets[ijet2].id_TLV)
                    for (int ibtag = 0; ibtag < btagWP::Count; ibtag++)
                        for (int jbtag = 0; jbtag < btagWP::Count; jbtag++)
                            if (map_jj[jetID::no * bitD + jetID::no  * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].size() > 0)
                                if (map_jj[jetID::no * bitD + jetID::no  * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].back() == count)
                                    map_jj[jetID::T  * bitD + jetID::TLV * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].push_back(count);
            }
            if (jets[ijet1].id_TLV)
            {
                for (int ibtag = 0; ibtag < btagWP::Count; ibtag++)
                    for (int jbtag = 0; jbtag < btagWP::Count; jbtag++)
                        if (map_jj[jetID::no  * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].size() > 0)
                            if (map_jj[jetID::no  * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].back() == count)
                                map_jj[jetID::TLV * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].push_back(count);
                if (jets[ijet2].id_L)
                    for (int ibtag = 0; ibtag < btagWP::Count; ibtag++)
                        for (int jbtag = 0; jbtag < btagWP::Count; jbtag++)
                            if (map_jj[jetID::no  * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].size() > 0)
                                if (map_jj[jetID::no  * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].back() == count)
                                    map_jj[jetID::TLV * bitD + jetID::L  * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].push_back(count);
                if (jets[ijet2].id_T)
                    for (int ibtag = 0; ibtag < btagWP::Count; ibtag++)
                        for (int jbtag = 0; jbtag < btagWP::Count; jbtag++)
                            if (map_jj[jetID::no  * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].size() > 0)
                                if (map_jj[jetID::no  * bitD + jetID::no * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].back() == count)
                                    map_jj[jetID::TLV * bitD + jetID::T  * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].push_back(count);
                if (jets[ijet2].id_TLV)
                    for (int ibtag = 0; ibtag < btagWP::Count; ibtag++)
                        for (int jbtag = 0; jbtag < btagWP::Count; jbtag++)
                            if (map_jj[jetID::no  * bitD + jetID::no  * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].size() > 0)
                                if (map_jj[jetID::no  * bitD + jetID::no  * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].back() == count)
                                    map_jj[jetID::TLV * bitD + jetID::TLV * bitC + ibtag * bitB + jbtag * bitA + jetPair::ht].push_back(count);
            }
            count++;
        }
    }

    // Fill the rest of dijet combinatoric maps
    bitA = jetPair::Count;
    bitB = bitA * btagWP::Count;
    bitC = bitB * btagWP::Count;
    bitD = bitC * jetID::Count;
    for (int ijetid1 = 0; ijetid1 < jetID::Count; ijetid1++)
    {
        for (int ijetid2 = 0; ijetid2 < jetID::Count; ijetid2++)
        {
            for (int ibtag1 = 0; ibtag1 < btagWP::Count; ibtag1++)
            {
                for (int ibtag2 = 0; ibtag2 < btagWP::Count; ibtag2++)
                {
                    int iht = ijetid1 * bitD + ijetid2 * bitC + ibtag1 * bitB + ibtag2 * bitA + jetPair::ht;
                    int jpt = ijetid1 * bitD + ijetid2 * bitC + ibtag1 * bitB + ibtag2 * bitA + jetPair::pt;
                    int jmh = ijetid1 * bitD + ijetid2 * bitC + ibtag1 * bitB + ibtag2 * bitA + jetPair::mh;
                    int jcsv = ijetid1 * bitD + ijetid2 * bitC + ibtag1 * bitB + ibtag2 * bitA + jetPair::csv;
                    int jjp = ijetid1 * bitD + ijetid2 * bitC + ibtag1 * bitB + ibtag2 * bitA + jetPair::jp;
                    int jptOverM = ijetid1 * bitD + ijetid2 * bitC + ibtag1 * bitB + ibtag2 * bitA + jetPair::ptOverM;
                    std::vector<int> tmp = map_jj[iht];
                    // do the ptjj sorted maps
                    std::sort(tmp.begin(), tmp.end(), [&](const int& a, const int& b){return jj[a].p4.Pt() > jj[b].p4.Pt();});
                    map_jj[jpt] = tmp;
                    // do the closest to mh sorted maps
                    tmp = map_jj[iht];
                    std::sort(tmp.begin(), tmp.end(), [&](const int& a, const int& b){return fabs(jj[a].p4.M() - mh) < fabs(jj[b].p4.M() - mh);});
                    map_jj[jmh] = tmp;
                    // do the sum(btag) sorted maps
                    tmp = map_jj[iht];
                    std::sort(tmp.begin(), tmp.end(), [&](const int& a, const int& b){return jj[a].sumCSV > jj[b].sumCSV;});
                    map_jj[jcsv] = tmp;
                    tmp = map_jj[iht];
                    std::sort(tmp.begin(), tmp.end(), [&](const int& a, const int& b){return jj[a].sumJP > jj[b].sumJP;});
                    map_jj[jjp] = tmp;
                    // do the ptjj / mjj sorted maps
                    tmp = map_jj[iht];
                    std::sort(tmp.begin(), tmp.end(), [&](const int& a, const int& b){return (jj[a].p4.Pt() / jj[a].p4.M()) > (jj[b].p4.Pt() / jj[b].p4.M());});
                    map_jj[jptOverM] = tmp;
                }
            }
        }
    }
 
    // ********** 
    // lljj, llbb, +pf_met
    // ********** 
    llmetjj.clear();
    for (unsigned int illmet = 0; illmet < llmet.size(); illmet++)
    {
        for (unsigned int ijj = 0; ijj < jj.size(); ijj++)
        {
            unsigned int imet = llmet[illmet].imet;
            unsigned int ill = llmet[illmet].ill;
            unsigned int ijet1 = jj[ijj].ijet1;
            unsigned int ijet2 = jj[ijj].ijet2;
            unsigned int ilep1 = ll[ill].ilep1;
            unsigned int ilep2 = ll[ill].ilep2;
            HH::DileptonMetDijet myllmetjj;
            myllmetjj.p4 = ll[ill].p4 + jj[ijj].p4 + met[imet].p4;
            myllmetjj.lep1_p4 = leptons[ilep1].p4;
            myllmetjj.lep2_p4 = leptons[ilep2].p4;
            myllmetjj.jet1_p4 = leptons[ijet1].p4;
            myllmetjj.jet2_p4 = leptons[ijet2].p4;
            myllmetjj.met_p4 = met[imet].p4;
            myllmetjj.ll_p4 = ll[ill].p4;
            myllmetjj.jj_p4 = jj[ijj].p4;
            myllmetjj.lljj_p4 = ll[ill].p4 + jj[ijj].p4;
            // gen info
            myllmetjj.gen_matched = ll[ill].gen_matched && jj[ijj].gen_matched && met[imet].gen_matched;
            myllmetjj.gen_p4 = myllmetjj.gen_matched ? ll[ill].gen_p4 + jj[ijj].gen_p4 + met[imet].gen_p4 : null_p4;
            myllmetjj.gen_DR = myllmetjj.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myllmetjj.p4, myllmetjj.gen_p4) : -1.;
            myllmetjj.gen_DPhi = myllmetjj.gen_matched ? fabs(ROOT::Math::VectorUtil::DeltaPhi(myllmetjj.p4, myllmetjj.gen_p4)) : -1.;
            myllmetjj.gen_DPtOverPt = myllmetjj.gen_matched ? (myllmetjj.p4.Pt() - myllmetjj.gen_p4.Pt()) / myllmetjj.p4.Pt() : -10.;
            myllmetjj.gen_lep1_p4 = leptons[ilep1].gen_p4;
            myllmetjj.gen_lep2_p4 = leptons[ilep2].gen_p4;
            myllmetjj.gen_jet1_p4 = leptons[ijet1].gen_p4;
            myllmetjj.gen_jet2_p4 = leptons[ijet2].gen_p4;
            myllmetjj.gen_met_p4 = met[imet].gen_p4;
            myllmetjj.gen_ll_p4 = ll[ill].gen_p4;
            myllmetjj.gen_jj_p4 = jj[ijj].gen_p4;
            myllmetjj.gen_lljj_p4 = ll[ill].gen_p4 + jj[ijj].gen_p4;
            // blind copy of the jj content
            myllmetjj.ijet1 = jj[ijj].ijet1;
            myllmetjj.ijet2 = jj[ijj].ijet2;
            myllmetjj.btag_LL = jj[ijj].btag_LL;
            myllmetjj.btag_LM = jj[ijj].btag_LM;
            myllmetjj.btag_LT = jj[ijj].btag_LT;
            myllmetjj.btag_ML = jj[ijj].btag_ML;
            myllmetjj.btag_MM = jj[ijj].btag_MM;
            myllmetjj.btag_MT = jj[ijj].btag_MT;
            myllmetjj.btag_TL = jj[ijj].btag_TL;
            myllmetjj.btag_TM = jj[ijj].btag_TM;
            myllmetjj.btag_TT = jj[ijj].btag_TT;
            myllmetjj.sumCSV = jj[ijj].sumCSV;
            myllmetjj.sumJP = jj[ijj].sumJP;
            myllmetjj.DR_j_j = jj[ijj].DR_j_j;
            myllmetjj.DPhi_j_j = jj[ijj].DPhi_j_j;
            myllmetjj.gen_matched_bbPartons = jj[ijj].gen_matched_bbPartons;
            myllmetjj.gen_matched_bbHadrons = jj[ijj].gen_matched_bbHadrons;
            myllmetjj.gen_bb = jj[ijj].gen_bb;
            myllmetjj.gen_bc = jj[ijj].gen_bc;
            myllmetjj.gen_bl = jj[ijj].gen_bl;
            myllmetjj.gen_cc = jj[ijj].gen_cc;
            myllmetjj.gen_cl = jj[ijj].gen_cl;
            myllmetjj.gen_ll = jj[ijj].gen_ll;
            // blind copy of the llmet content
            myllmetjj.ilep1 = ll[ill].ilep1;
            myllmetjj.ilep2 = ll[ill].ilep2;
            myllmetjj.isOS = ll[ill].isOS;
            myllmetjj.isMuMu = ll[ill].isMuMu;
            myllmetjj.isElEl = ll[ill].isElEl;
            myllmetjj.isElMu = ll[ill].isElMu;
            myllmetjj.isMuEl = ll[ill].isMuEl;
            myllmetjj.isSF = ll[ill].isSF;
            myllmetjj.id_LL = ll[ill].id_LL;
            myllmetjj.id_LT = ll[ill].id_LT;
            myllmetjj.id_TL = ll[ill].id_TL;
            myllmetjj.id_TT = ll[ill].id_TT;
            myllmetjj.iso_LL = ll[ill].iso_LL;
            myllmetjj.iso_LT = ll[ill].iso_LT;
            myllmetjj.iso_TL = ll[ill].iso_TL;
            myllmetjj.iso_TT = ll[ill].iso_TT;
            myllmetjj.DR_l_l = ll[ill].DR_l_l;
            myllmetjj.DPhi_l_l = ll[ill].DPhi_l_l;
            myllmetjj.ill = ill;
            myllmetjj.imet = imet;
            myllmetjj.isNoHF = met[imet].isNoHF;
            myllmetjj.DPhi_ll_met = llmet[illmet].DPhi_ll_met;
            myllmetjj.minDPhi_l_met = llmet[illmet].minDPhi_l_met; 
            myllmetjj.maxDPhi_l_met = llmet[illmet].maxDPhi_l_met;
            myllmetjj.MT = llmet[illmet].MT;
            myllmetjj.MT_formula = llmet[illmet].MT_formula;
            myllmetjj.projectedMet = llmet[illmet].projectedMet;
            // content specific to HH::DijetMet
            // NB: computed for the first time here, no intermediate jjmet collection
            myllmetjj.DPhi_jj_met = fabs(ROOT::Math::VectorUtil::DeltaPhi(jj[ijj].p4, met[imet].p4));
            myllmetjj.minDPhi_j_met = std::min(fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[jj[ijj].ijet1].p4, met[imet].p4)), fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[jj[ijj].ijet2].p4, met[imet].p4)));
            myllmetjj.maxDPhi_j_met = std::max(fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[jj[ijj].ijet1].p4, met[imet].p4)), fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[jj[ijj].ijet2].p4, met[imet].p4)));
            // content specific to HH::DileptonMetDijet
            myllmetjj.illmet = illmet;
            myllmetjj.ijj = ijj;
            float DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2;
            DR_j1l1 = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, leptons[ilep1].p4);
            DR_j1l2 = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, leptons[ilep2].p4);
            DR_j2l1 = ROOT::Math::VectorUtil::DeltaR(jets[ijet2].p4, leptons[ilep1].p4);
            DR_j2l2 = ROOT::Math::VectorUtil::DeltaR(jets[ijet2].p4, leptons[ilep2].p4);
            myllmetjj.maxDR_l_j = std::max({DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2});
            myllmetjj.minDR_l_j = std::min({DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2});
            myllmetjj.DR_ll_jj = ROOT::Math::VectorUtil::DeltaR(ll[ill].p4, jj[ijj].p4);
            myllmetjj.DPhi_ll_jj = fabs(ROOT::Math::VectorUtil::DeltaPhi(ll[ill].p4, jj[ijj].p4));
            myllmetjj.DR_llmet_jj = ROOT::Math::VectorUtil::DeltaR(llmet[illmet].p4, jj[ijj].p4);
            myllmetjj.DPhi_llmet_jj = fabs(ROOT::Math::VectorUtil::DeltaPhi(llmet[illmet].p4, jj[ijj].p4));
            myllmetjj.cosThetaStar_CS = fabs(getCosThetaStar_CS(llmet[illmet].p4, jj[ijj].p4));
            myllmetjj.MT_fullsystem = myllmetjj.p4.Mt();
            if (myllmetjj.minDR_l_j < m_minDR_l_j_Cut)
                continue;
            llmetjj.push_back(myllmetjj);
        }
    }

    if (HHANADEBUG) {std::cout << "##### Now debugging the llmetjj maps" << std::endl;}
    // llmetjj maps: cross product of llmet and jj maps
    for (unsigned int i = 0; i < map_llmetjj.size(); i++)
        map_llmetjj[i].clear();

    bitA = jetPair::Count;
    bitB = bitA * btagWP::Count;
    bitC = bitB * btagWP::Count;
    bitD = bitC * jetID::Count;
    int bitE = bitD * jetID::Count;
    int bitF = bitE * lepIso::Count;
    int bitG = bitF * lepID::Count;
    int bitH = bitG * lepIso::Count;

    int bit0 = lepIso::Count;
    int bit1 = bit0 * lepID::Count;
    int bit2 = bit1 * lepIso::Count;

    int bitX = jetPair::Count;
    int bitY = bitX * btagWP::Count;
    int bitZ = bitY * btagWP::Count;
    int bitXX = bitZ * jetID::Count;
    for (int il1id = 0; il1id < lepID::Count; il1id++)
    {
        for (int il1iso = 0; il1iso < lepIso::Count; il1iso++)
        {
            for (int il2id = 0; il2id < lepID::Count; il2id++)
            {
                for (int il2iso = 0; il2iso < lepIso::Count; il2iso++)
                {
                    for (int ij1id = 0; ij1id < jetID::Count; ij1id++)
                    {
                        for (int ij2id = 0; ij2id < jetID::Count; ij2id++)
                        {
                            for(int ibtag1 = 0; ibtag1 < btagWP::Count; ibtag1++)
                            {
                                for(int ibtag2 = 0; ibtag2 < btagWP::Count; ibtag2++)
                                {
                                    for(int ipair = 0; ipair < jetPair::Count; ipair++)
                                    {
                                        int illmetjj = il1id * bitH
                                            + il1iso * bitG
                                            + il2id * bitF
                                            + il2iso * bitE
                                            + ij1id * bitD
                                            + ij2id * bitC
                                            + ibtag1 * bitB
                                            + ibtag2 * bitA
                                            + ipair;
                                        int illmet = il1id * bit2
                                            + il1iso * bit1
                                            + il2id * bit0
                                            + il2iso;
                                        int ijj = ij1id * bitXX
                                            + ij2id * bitZ
                                            + ibtag1 * bitY
                                            + ibtag2 * bitX
                                            + ipair;
                                        if (map_llmet[illmet].size() == 0 || map_jj[ijj].size() == 0)
                                            continue;
                                        if (HHANADEBUG)
                                        {
                                            std::cout << "Now treating illmetjj= " << illmetjj << " (illmet, ijj)= (" << illmet << ", " << ijj << ")" << std::endl;
                                            std::cout << "map_llmet[" << illmet << "].size()= " << map_llmet[illmet].size() << std::endl;
                                            for (unsigned int j = 0; j < map_llmet[illmet].size(); j++)
                                                std::cout << "\tmap_llmet[" << illmet << "][" << j << "]= " << map_llmet[illmet][j] << std::endl;
                                            std::cout << "map_jj[" << ijj << "].size()= " << map_jj[ijj].size() << std::endl;
                                            for (unsigned int j = 0; j < map_jj[ijj].size(); j++)
                                                std::cout << "\tmap_jj[" << ijj << "][" << j << "]= " << map_jj[ijj][j] << std::endl;
                                        }
                                        for (unsigned int jjj = 0; jjj < map_jj[ijj].size(); jjj++)
                                        {
                                            for (unsigned int i = 0; i < llmetjj.size(); i++)
                                            {
                                                if (std::find(map_llmet[illmet].begin(), map_llmet[illmet].end(), llmetjj[i].illmet) == map_llmet[illmet].end())
                                                    continue;
                                                if (llmetjj[i].ijj == map_jj[ijj][jjj])
                                                {
                                                    map_llmetjj[illmetjj].push_back(i);
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // Adding some few custom candidates, for convenience
    // llmetjj_allTight_btagL
    int icustom = lepID::T * bitH
        + lepIso::T * bitG
        + lepID::T  * bitF
        + lepIso::T * bitE
        + jetID::T  * bitD
        + jetID::T  * bitC
        + btagWP::L * bitB
        + btagWP::L * bitA
        + jetPair::ht;
    llmetjj_allTight_btagL_ht.clear();
    for (unsigned int icandidate = 0; icandidate < map_llmetjj[icustom].size(); icandidate++)
        llmetjj_allTight_btagL_ht.push_back(llmetjj[map_llmetjj[icustom][icandidate]]);

    icustom = icustom - jetPair::ht + jetPair::pt;
    llmetjj_allTight_btagL_pt.clear();
    for (unsigned int icandidate = 0; icandidate < map_llmetjj[icustom].size(); icandidate++)
        llmetjj_allTight_btagL_pt.push_back(llmetjj[map_llmetjj[icustom][icandidate]]);

    icustom = icustom - jetPair::pt + jetPair::csv;
    llmetjj_allTight_btagL_csv.clear();
    for (unsigned int icandidate = 0; icandidate < map_llmetjj[icustom].size(); icandidate++)
        llmetjj_allTight_btagL_csv.push_back(llmetjj[map_llmetjj[icustom][icandidate]]);

    // llmetjj_allTight_nobtag
    icustom = lepID::T * bitH
        + lepIso::T  * bitG
        + lepID::T   * bitF
        + lepIso::T  * bitE
        + jetID::T   * bitD
        + jetID::T   * bitC
        + btagWP::no * bitB
        + btagWP::no * bitA
        + jetPair::ht;
    llmetjj_allTight_nobtag_ht.clear();
    for (unsigned int icandidate = 0; icandidate < map_llmetjj[icustom].size(); icandidate++)
         llmetjj_allTight_nobtag_ht.push_back(llmetjj[map_llmetjj[icustom][icandidate]]);

    icustom = icustom - jetPair::ht + jetPair::pt;
    llmetjj_allTight_nobtag_pt.clear();
    for (unsigned int icandidate = 0; icandidate < map_llmetjj[icustom].size(); icandidate++)
        llmetjj_allTight_nobtag_pt.push_back(llmetjj[map_llmetjj[icustom][icandidate]]);

    icustom = icustom - jetPair::pt + jetPair::csv;
    llmetjj_allTight_nobtag_csv.clear();
    for (unsigned int icandidate = 0; icandidate < map_llmetjj[icustom].size(); icandidate++)
        llmetjj_allTight_nobtag_csv.push_back(llmetjj[map_llmetjj[icustom][icandidate]]);


    // ***** ***** *****
    // Event variables
    // ***** ***** *****
    nJets = jets.size();
    nJetsL = 0;
    for (unsigned int ijet = 0; ijet < jets.size(); ijet++)
        if (jets[ijet].id_L)
            nJetsL++;
    nBJetsL = 0;
    nBJetsM = 0;
    nBJetsT = 0;
    for (unsigned int ijet = 0; ijet < jets.size(); ijet++)
    {
        if (!jets[ijet].id_L) continue;
        if (jets[ijet].btag_L)
            nBJetsL++;
        if (jets[ijet].btag_M)
            nBJetsM++;
        if (jets[ijet].btag_T)
            nBJetsT++;
    }
    nMuons = muons.size();
    nMuonsL = 0;
    nMuonsT = 0;
    nElectrons = electrons.size();
    nElectronsL = 0;
    nElectronsT = 0;
    nLeptons = leptons.size();
    nLeptonsL = 0;
    nLeptonsT = 0;
    for (unsigned int ilepton = 0; ilepton < leptons.size(); ilepton++)
    {
        if (leptons[ilepton].id_L && leptons[ilepton].iso_L)
        {
            nLeptonsL++;
            if (leptons[ilepton].isMu)
                nMuonsL++;
            if (leptons[ilepton].isEl)
                nElectronsL++;
        }
        if (leptons[ilepton].id_T && leptons[ilepton].iso_T)
        {
            nLeptonsT++;
            if (leptons[ilepton].isMu)
                nMuonsT++;
            if (leptons[ilepton].isEl)
                nElectronsT++;
        }
    }



    if (!event.isRealData())
    {
// ***** ***** *****
// Get the MC truth information on the hard process
// ***** ***** *****
// from https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HepMCCandidate/interface/GenStatusFlags.h
//    enum StatusBits {
//0      kIsPrompt = 0,
//1      kIsDecayedLeptonHadron,
//2      kIsTauDecayProduct,
//3      kIsPromptTauDecayProduct,
//4      kIsDirectTauDecayProduct,
//5      kIsDirectPromptTauDecayProduct,
//6      kIsDirectHadronDecayProduct,
//7      kIsHardProcess,
//8      kFromHardProcess,
//9      kIsHardProcessTauDecayProduct,
//10      kIsDirectHardProcessTauDecayProduct,
//11      kFromHardProcessBeforeFSR,
//12      kIsFirstCopy,
//13      kIsLastCopy,
//14      kIsLastCopyBeforeFSR
//    };

        const GenParticlesProducer& gp = producers.get<GenParticlesProducer>("gen_particles");
    
        BRANCH(gen_B1, LorentzVector);
        BRANCH(gen_B2, LorentzVector);
        BRANCH(gen_B1FSR, LorentzVector);
        BRANCH(gen_B2FSR, LorentzVector);
        BRANCH(gen_Nu1, LorentzVector);
        BRANCH(gen_Nu2, LorentzVector);
        BRANCH(gen_L1, LorentzVector);
        BRANCH(gen_L2, LorentzVector);
        BRANCH(gen_LL, LorentzVector);
        BRANCH(gen_BB, LorentzVector);
        BRANCH(gen_BBFSR, LorentzVector);
        BRANCH(gen_L1FSRNu, LorentzVector);
        BRANCH(gen_L2FSRNu, LorentzVector);
        BRANCH(gen_L1FSR, LorentzVector);
        BRANCH(gen_L2FSR, LorentzVector);
        BRANCH(gen_LLFSR, LorentzVector);
        BRANCH(gen_NuNu, LorentzVector);
        BRANCH(gen_LLNuNu, LorentzVector);
        BRANCH(gen_LLFSRNuNu, LorentzVector);
        BRANCH(gen_LLFSRNuNuBB, LorentzVector);
        gen_B1.SetPxPyPzE(0.,0.,0.,0.);
        gen_B2.SetPxPyPzE(0.,0.,0.,0.);
        gen_B1FSR.SetPxPyPzE(0.,0.,0.,0.);
        gen_B2FSR.SetPxPyPzE(0.,0.,0.,0.);
        gen_Nu1.SetPxPyPzE(0.,0.,0.,0.);
        gen_Nu2.SetPxPyPzE(0.,0.,0.,0.);
        gen_L1.SetPxPyPzE(0.,0.,0.,0.);
        gen_L2.SetPxPyPzE(0.,0.,0.,0.);
        gen_LL.SetPxPyPzE(0.,0.,0.,0.);
        gen_BB.SetPxPyPzE(0.,0.,0.,0.);
        gen_BBFSR.SetPxPyPzE(0.,0.,0.,0.);
        gen_L1FSRNu.SetPxPyPzE(0.,0.,0.,0.);
        gen_L2FSRNu.SetPxPyPzE(0.,0.,0.,0.);
        gen_L1FSR.SetPxPyPzE(0.,0.,0.,0.);
        gen_L2FSR.SetPxPyPzE(0.,0.,0.,0.);
        gen_LLFSR.SetPxPyPzE(0.,0.,0.,0.);
        gen_NuNu.SetPxPyPzE(0.,0.,0.,0.);
        gen_LLNuNu.SetPxPyPzE(0.,0.,0.,0.);
        gen_LLFSRNuNu.SetPxPyPzE(0.,0.,0.,0.);
        gen_LLFSRNuNuBB.SetPxPyPzE(0.,0.,0.,0.);
    
        BRANCH(gen_iX, int);
        BRANCH(gen_iH1, int);
        BRANCH(gen_iH2, int);
        BRANCH(gen_iV1, int);
        BRANCH(gen_iV2, int);
        BRANCH(gen_iB1, int);
        BRANCH(gen_iB2, int);
        BRANCH(gen_iL1, int);
        BRANCH(gen_iL2, int);
        BRANCH(gen_iNu1, int);
        BRANCH(gen_iNu2, int);
        BRANCH(gen_iG1, std::vector<int>);
        BRANCH(gen_iG2, std::vector<int>);
        BRANCH(gen_iGlu1, std::vector<int>);
        BRANCH(gen_iGlu2, std::vector<int>);
        gen_iX = -1;
        gen_iH1 = gen_iH2 = -1;
        gen_iV1 = gen_iV2 = -1;
        gen_iB1 = gen_iB2 = -1;
        gen_iL1 = gen_iL2 = -1;
        gen_iNu1 = gen_iNu2 = -1;
        gen_iG1.clear();
        gen_iG2.clear();
        gen_iGlu1.clear();
        gen_iGlu2.clear();
        bool isThereFSRforL1 = false;
        bool isThereFSRforL2 = false;
        bool isThereFSRforB1 = false;
        bool isThereFSRforB2 = false;
    
        // Construct the signal gen-level info
        for (unsigned int ip = 0; ip < gp.pruned_p4.size(); ip++) {
            std::bitset<15> flags (gp.pruned_status_flags[ip]);
            if (!flags.test(14)) continue; // take the last copies before FSR
            if (flags.test(8)) // first look at the hard process
            {
                if (HHANADEBUG)
                    std::cout << "\tip= " << ip << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[ip] << "\tflags= " << flags << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[ip].Pt() << ", " << gp.pruned_p4[ip].Eta() << ", " << gp.pruned_p4[ip].Phi() << ", " << gp.pruned_p4[ip].E() << ", " << gp.pruned_p4[ip].M() << ")" << std::endl;
                if (abs(gp.pruned_pdg_id[ip]) == 25)
                {
                    if (gen_iH1 == -1)
                        gen_iH1 = ip;
                    else if (gen_iH2 == -1)
                        gen_iH2 = ip;
                }
                else if (abs(gp.pruned_pdg_id[ip]) == 5)
                {
                    if (gen_iB1 == -1)
                    {
                        if (flags.test(8) && !flags.test(7))
                            isThereFSRforB1 = true;
                        gen_iB1 = ip;
                    }
                    else if (gen_iB2 == -1)
                    {
                        if (flags.test(8) && !flags.test(7))
                            isThereFSRforB2 = true;
                        gen_iB2 = ip;
                    }
                }
                else if (abs(gp.pruned_pdg_id[ip]) == 11 || abs(gp.pruned_pdg_id[ip]) == 13)
                {
                // if the lepton is not in the hard process itself, then look for FSR photons later on
                    if (gen_iL1 == -1)
                    {
                        if (flags.test(8) && !flags.test(7))
                            isThereFSRforL1 = true;
                        gen_iL1 = ip;
                    }
                    else if (gen_iL2 == -1)
                    {
                        if (flags.test(8) && !flags.test(7))
                            isThereFSRforL2 = true;
                        gen_iL2 = ip;
                    }
                }
                else if (abs(gp.pruned_pdg_id[ip]) == 12 || abs(gp.pruned_pdg_id[ip]) == 14 || abs(gp.pruned_pdg_id[ip]) == 16)
                {
                    if (gen_iNu1 == -1)
                        gen_iNu1 = ip;
                    else if (gen_iNu2 == -1)
                        gen_iNu2 = ip;
                }
                else if (abs(gp.pruned_pdg_id[ip]) == 35 || abs(gp.pruned_pdg_id[ip]) == 39)
                {
                    if (gen_iX == -1)
                        gen_iX = ip;
                }
                else if (abs(gp.pruned_pdg_id[ip]) == 23 || abs(gp.pruned_pdg_id[ip]) == 24)
                {
                    if (gen_iV1 == -1)
                        gen_iV1 = ip;
                    else if (gen_iV2 == -1)
                        gen_iV2 = ip;
                }
            } // end if coming from hard process
        } // end of loop over pruned gen particles
        // new loop to find FSR photons
        if (isThereFSRforL1 || isThereFSRforL2)
        {
            for (unsigned int ip = 0; ip < gp.pruned_p4.size(); ip++) {
                if (gp.pruned_pdg_id[ip] != 22) continue;
                std::bitset<15> flags (gp.pruned_status_flags[ip]);
                if (!flags.test(14)) continue; // take the last copies
                for (unsigned int imother = 0; imother < gp.pruned_mothers_index[ip].size(); imother++)
                {
                    if (isThereFSRforL1)
                    {
                        for (unsigned int jmother = 0; jmother < gp.pruned_mothers_index[gen_iL1].size(); jmother++)
                        {
                            if (gp.pruned_mothers_index[ip][imother] == gp.pruned_mothers_index[gen_iL1][jmother])
                                gen_iG1.push_back(ip);
                        }
                    }
                    if (isThereFSRforL2) 
                    {
                        for (unsigned int jmother = 0; jmother < gp.pruned_mothers_index[gen_iL2].size(); jmother++)
                        {
                            if (gp.pruned_mothers_index[ip][imother] == gp.pruned_mothers_index[gen_iL2][jmother])
                                gen_iG2.push_back(ip);
                        }
                    }
                }
            } // end of loop over pruned gen particles
        } // end of FSR loop
        // new loop to find FSR gluons
        if (isThereFSRforB1 || isThereFSRforB2)
        {
            for (unsigned int ip = 0; ip < gp.pruned_p4.size(); ip++) {
                if (gp.pruned_pdg_id[ip] != 21) continue;
                std::bitset<15> flags (gp.pruned_status_flags[ip]);
                if (!flags.test(14)) continue; // take the last copies
                for (unsigned int imother = 0; imother < gp.pruned_mothers_index[ip].size(); imother++)
                {
                    if (isThereFSRforB1)
                    {
                        for (unsigned int jmother = 0; jmother < gp.pruned_mothers_index[gen_iB1].size(); jmother++)
                        {
                            if (gp.pruned_mothers_index[ip][imother] == gp.pruned_mothers_index[gen_iB1][jmother])
                                gen_iGlu1.push_back(ip);
                        }
                    }
                    if (isThereFSRforB2) 
                    {
                        for (unsigned int jmother = 0; jmother < gp.pruned_mothers_index[gen_iB2].size(); jmother++)
                        {
                            if (gp.pruned_mothers_index[ip][imother] == gp.pruned_mothers_index[gen_iB2][jmother])
                                gen_iGlu2.push_back(ip);
                        }
                    }
                }
            } // end of loop over pruned gen particles
        } // end of FSR loop
    
        if (HHANADEBUG)
        {
            if (gen_iX != -1)
                std::cout << "\tiX= " << gen_iX << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iX] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iX]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iX].Pt() << ", " << gp.pruned_p4[gen_iX].Eta() << ", " << gp.pruned_p4[gen_iX].Phi() << ", " << gp.pruned_p4[gen_iX].E() << ", " << gp.pruned_p4[gen_iX].M() << ")" << std::endl;
            std::cout << "\tiH1= " << gen_iH1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iH1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iH1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iH1].Pt() << ", " << gp.pruned_p4[gen_iH1].Eta() << ", " << gp.pruned_p4[gen_iH1].Phi() << ", " << gp.pruned_p4[gen_iH1].E() << ", " << gp.pruned_p4[gen_iH1].M() << ")" << std::endl;
            std::cout << "\tiH2= " << gen_iH2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iH2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iH2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iH2].Pt() << ", " << gp.pruned_p4[gen_iH2].Eta() << ", " << gp.pruned_p4[gen_iH2].Phi() << ", " << gp.pruned_p4[gen_iH2].E() << ", " << gp.pruned_p4[gen_iH2].M() << ")" << std::endl;
            std::cout << "\tiB1= " << gen_iB1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iB1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iB1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iB1].Pt() << ", " << gp.pruned_p4[gen_iB1].Eta() << ", " << gp.pruned_p4[gen_iB1].Phi() << ", " << gp.pruned_p4[gen_iB1].E() << ", " << gp.pruned_p4[gen_iB1].M() << ")" << std::endl;
            std::cout << "\tiB2= " << gen_iB2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iB2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iB2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iB2].Pt() << ", " << gp.pruned_p4[gen_iB2].Eta() << ", " << gp.pruned_p4[gen_iB2].Phi() << ", " << gp.pruned_p4[gen_iB2].E() << ", " << gp.pruned_p4[gen_iB2].M() << ")" << std::endl;
            if (gen_iV1 != -1)
                std::cout << "\tiV1= " << gen_iV1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iV1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iV1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iV1].Pt() << ", " << gp.pruned_p4[gen_iV1].Eta() << ", " << gp.pruned_p4[gen_iV1].Phi() << ", " << gp.pruned_p4[gen_iV1].E() << ", " << gp.pruned_p4[gen_iV1].M() << ")" << std::endl;
            if (gen_iV2 != -1)
                std::cout << "\tiV2= " << gen_iV2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iV2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iV2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iV2].Pt() << ", " << gp.pruned_p4[gen_iV2].Eta() << ", " << gp.pruned_p4[gen_iV2].Phi() << ", " << gp.pruned_p4[gen_iV2].E() << ", " << gp.pruned_p4[gen_iV2].M() << ")" << std::endl;
            std::cout << "\tiL1= " << gen_iL1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iL1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iL1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iL1].Pt() << ", " << gp.pruned_p4[gen_iL1].Eta() << ", " << gp.pruned_p4[gen_iL1].Phi() << ", " << gp.pruned_p4[gen_iL1].E() << ", " << gp.pruned_p4[gen_iL1].M() << ")" << std::endl;
            std::cout << "\tiL2= " << gen_iL2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iL2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iL2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iL2].Pt() << ", " << gp.pruned_p4[gen_iL2].Eta() << ", " << gp.pruned_p4[gen_iL2].Phi() << ", " << gp.pruned_p4[gen_iL2].E() << ", " << gp.pruned_p4[gen_iL2].M() << ")" << std::endl;
            std::cout << "\tiNu1= " << gen_iNu1 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iNu1] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iNu1]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iNu1].Pt() << ", " << gp.pruned_p4[gen_iNu1].Eta() << ", " << gp.pruned_p4[gen_iNu1].Phi() << ", " << gp.pruned_p4[gen_iNu1].E() << ", " << gp.pruned_p4[gen_iNu1].M() << ")" << std::endl;
            std::cout << "\tiNu2= " << gen_iNu2 << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iNu2] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iNu2]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iNu2].Pt() << ", " << gp.pruned_p4[gen_iNu2].Eta() << ", " << gp.pruned_p4[gen_iNu2].Phi() << ", " << gp.pruned_p4[gen_iNu2].E() << ", " << gp.pruned_p4[gen_iNu2].M() << ")" << std::endl;
        
            for (unsigned int iG = 0; iG < gen_iG1.size(); iG++)
            {
                std::cout << "\tiG1[" << iG << "]= " << gen_iG1[iG] << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iG1[iG]] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iG1[iG]]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iG1[iG]].Pt() << ", " << gp.pruned_p4[gen_iG1[iG]].Eta() << ", " << gp.pruned_p4[gen_iG1[iG]].Phi() << ", " << gp.pruned_p4[gen_iG1[iG]].E() << ", " << gp.pruned_p4[gen_iG1[iG]].M() << ")" << std::endl;
            }
            for (unsigned int iG = 0; iG < gen_iG2.size(); iG++)
            {
                std::cout << "\tiG2[" << iG << "]= " << gen_iG2[iG] << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iG2[iG]] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iG2[iG]]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iG2[iG]].Pt() << ", " << gp.pruned_p4[gen_iG2[iG]].Eta() << ", " << gp.pruned_p4[gen_iG2[iG]].Phi() << ", " << gp.pruned_p4[gen_iG2[iG]].E() << ", " << gp.pruned_p4[gen_iG2[iG]].M() << ")" << std::endl;
            }
            for (unsigned int iGlu = 0; iGlu < gen_iGlu1.size(); iGlu++)
            {
                std::cout << "\tiGlu1[" << iGlu << "]= " << gen_iGlu1[iGlu] << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iGlu1[iGlu]] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iGlu1[iGlu]]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iGlu1[iGlu]].Pt() << ", " << gp.pruned_p4[gen_iGlu1[iGlu]].Eta() << ", " << gp.pruned_p4[gen_iGlu1[iGlu]].Phi() << ", " << gp.pruned_p4[gen_iGlu1[iGlu]].E() << ", " << gp.pruned_p4[gen_iGlu1[iGlu]].M() << ")" << std::endl;
            }
            for (unsigned int iGlu = 0; iGlu < gen_iGlu2.size(); iGlu++)
            {
                std::cout << "\tiGlu2[" << iGlu << "]= " << gen_iGlu2[iGlu] << "\tgp.pruned_pdg_id= " << gp.pruned_pdg_id[gen_iGlu2[iGlu]] << "\tflags= " << std::bitset<15>(gp.pruned_status_flags[gen_iGlu2[iGlu]]) << "\t(pt, eta, phi, e, mass)= (" << gp.pruned_p4[gen_iGlu2[iGlu]].Pt() << ", " << gp.pruned_p4[gen_iGlu2[iGlu]].Eta() << ", " << gp.pruned_p4[gen_iGlu2[iGlu]].Phi() << ", " << gp.pruned_p4[gen_iGlu2[iGlu]].E() << ", " << gp.pruned_p4[gen_iGlu2[iGlu]].M() << ")" << std::endl;
            }
        } // end if HHANADEBUG
    
    
        gen_B1 += gp.pruned_p4[gen_iB1];
        gen_B2 += gp.pruned_p4[gen_iB2];
        gen_Nu1 += gp.pruned_p4[gen_iNu1];
        gen_Nu2 += gp.pruned_p4[gen_iNu2];
        gen_L1 += gp.pruned_p4[gen_iL1];
        gen_L2 += gp.pruned_p4[gen_iL2];
        gen_LL += gen_L1 + gen_L2;
        gen_NuNu += gen_Nu1 + gen_Nu2;
        gen_B1FSR += gen_B1;
        gen_B2FSR += gen_B2;
        gen_L1FSR += gen_L1;
        gen_L2FSR += gen_L2;
        for (unsigned int iG = 0; iG < gen_iG1.size(); iG++)
            gen_L1FSR += gp.pruned_p4[gen_iG1[iG]];
        for (unsigned int iG = 0; iG < gen_iG2.size(); iG++)
            gen_L2FSR += gp.pruned_p4[gen_iG2[iG]];
        for (unsigned int iGlu = 0; iGlu < gen_iGlu1.size(); iGlu++)
            gen_B1FSR += gp.pruned_p4[gen_iGlu1[iGlu]];
        for (unsigned int iGlu = 0; iGlu < gen_iGlu2.size(); iGlu++)
            gen_B2FSR += gp.pruned_p4[gen_iGlu2[iGlu]];
        gen_BB += gen_B1 + gen_B2;
        gen_BBFSR += gen_B1FSR + gen_B2FSR;
        gen_L1FSRNu += gen_L1FSR + gen_Nu1;
        gen_L2FSRNu += gen_L2FSR + gen_Nu2;
        gen_LLFSR += gen_L1FSR + gen_L2FSR;
        gen_LLNuNu += gen_LL + gen_NuNu;
        gen_LLFSRNuNu += gen_LLFSR + gen_NuNu;
        gen_LLFSRNuNuBB += gen_LLFSRNuNu + gen_BBFSR;
    
        if (HHANADEBUG)
        {
            std::cout << "\tgen_B1.M()= " << gen_B1.M() << std::endl;
            std::cout << "\tgen_B2.M()= " << gen_B2.M() << std::endl;
            std::cout << "\tgen_B1FSR.M()= " << gen_B1FSR.M() << std::endl;
            std::cout << "\tgen_B2FSR.M()= " << gen_B2FSR.M() << std::endl;
            std::cout << "\tgen_Nu1.M()= " << gen_Nu1.M() << std::endl;
            std::cout << "\tgen_Nu2.M()= " << gen_Nu2.M() << std::endl;
            std::cout << "\tgen_L1.M()= " << gen_L1.M() << std::endl;
            std::cout << "\tgen_L2.M()= " << gen_L2.M() << std::endl;
            std::cout << "\tgen_LL.M()= " << gen_LL.M() << std::endl;
            std::cout << "\tgen_L1FSRNu.M()= " << gen_L1FSRNu.M() << std::endl;
            std::cout << "\tgen_L2FSRNu.M()= " << gen_L2FSRNu.M() << std::endl;
            std::cout << "\tgen_L1FSR.M()= " << gen_L1FSR.M() << std::endl;
            std::cout << "\tgen_L2FSR.M()= " << gen_L2FSR.M() << std::endl;
            std::cout << "\tgen_LLFSR.M()= " << gen_LLFSR.M() << std::endl;
            std::cout << "\tgen_NuNu.M()= " << gen_NuNu.M() << std::endl;
            std::cout << "\tgen_LLNuNu.M()= " << gen_LLNuNu.M() << std::endl;
            std::cout << "\tgen_LLFSRNuNu.M()= " << gen_LLFSRNuNu.M() << std::endl;
            std::cout << "\tgen_BB.M()= " << gen_BB.M() << std::endl;
            std::cout << "\tgen_BBFSR.M()= " << gen_BBFSR.M() << std::endl;
            std::cout << "\tgen_LLFSRNuNuBB.M()= " << gen_LLFSRNuNuBB.M() << std::endl;
        }

        // ***** ***** *****
        // Matching
        // ***** ***** *****
        BRANCH(gen_deltaR_jet_B1, std::vector<float>);    
        BRANCH(gen_deltaR_jet_B2, std::vector<float>);    
        BRANCH(gen_deltaR_jet_B1FSR, std::vector<float>);    
        BRANCH(gen_deltaR_jet_B2FSR, std::vector<float>);    
        BRANCH(gen_deltaR_electron_L1, std::vector<float>);    
        BRANCH(gen_deltaR_electron_L2, std::vector<float>);    
        BRANCH(gen_deltaR_electron_L1FSR, std::vector<float>);    
        BRANCH(gen_deltaR_electron_L2FSR, std::vector<float>);    
        BRANCH(gen_deltaR_muon_L1, std::vector<float>);    
        BRANCH(gen_deltaR_muon_L2, std::vector<float>);    
        BRANCH(gen_deltaR_muon_L1FSR, std::vector<float>);    
        BRANCH(gen_deltaR_muon_L2FSR, std::vector<float>);    
    
        for (auto p4: alljets.gen_p4) {
            gen_deltaR_jet_B1.push_back(deltaR(p4, gen_B1));
            gen_deltaR_jet_B2.push_back(deltaR(p4, gen_B2));
            gen_deltaR_jet_B1FSR.push_back(deltaR(p4, gen_B1FSR));
            gen_deltaR_jet_B2FSR.push_back(deltaR(p4, gen_B2FSR));
        }
        for (auto p4: allelectrons.gen_p4) {
            gen_deltaR_electron_L1.push_back(deltaR(p4, gen_L1));
            gen_deltaR_electron_L2.push_back(deltaR(p4, gen_L2));
            gen_deltaR_electron_L1FSR.push_back(deltaR(p4, gen_L1FSR));
            gen_deltaR_electron_L2FSR.push_back(deltaR(p4, gen_L2FSR));
        }
        for (auto p4: allmuons.gen_p4) {
            gen_deltaR_muon_L1.push_back(deltaR(p4, gen_L1));
            gen_deltaR_muon_L2.push_back(deltaR(p4, gen_L2));
            gen_deltaR_muon_L1FSR.push_back(deltaR(p4, gen_L1FSR));
            gen_deltaR_muon_L2FSR.push_back(deltaR(p4, gen_L2FSR));
        }
    } // end of if !event.isRealData()

}

