#include <cp3_llbb/HHAnalysis/interface/HHAnalyzer.h>
#include <cp3_llbb/Framework/interface/BTagsAnalyzer.h>
#include <cp3_llbb/HHAnalysis/interface/Categories.h>
#include <cp3_llbb/HHAnalysis/interface/GenStatusFlags.h>

#include <cp3_llbb/Framework/interface/EventProducer.h>
#include <cp3_llbb/Framework/interface/GenParticlesProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/LeptonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/METProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cmath>

#define HH_GEN_DEBUG (false)
#define TT_GEN_DEBUG (false)

void HHAnalyzer::registerCategories(CategoryManager& manager, const edm::ParameterSet& config) {
    edm::ParameterSet newconfig = edm::ParameterSet(config);
    newconfig.addUntrackedParameter("m_analyzer_name", this->m_name);
    manager.new_category<MuMuCategory>("mumu", "Category with leading leptons as two muons", newconfig);
    manager.new_category<ElElCategory>("elel", "Category with leading leptons as two electrons", newconfig);
    manager.new_category<ElMuCategory>("elmu", "Category with leading leptons as electron, subleading as muon", newconfig);
    manager.new_category<MuElCategory>("muel", "Category with leading leptons as muon, subleading as electron", newconfig);
}


void HHAnalyzer::analyze(const edm::Event& event, const edm::EventSetup&, const ProducersManager& producers, const AnalyzersManager&, const CategoryManager&) {

    //float mh = event.isRealData() ? 125.09 : 125.0;
    LorentzVector null_p4(0., 0., 0., 0.);
    const EventProducer& fwevent = producers.get<EventProducer>("event");
    float event_weight = fwevent.weight;
    float tmp_count_has2leptons = 0.;
    float tmp_count_has2leptons_elel = 0.;
    float tmp_count_has2leptons_elmu = 0.;
    float tmp_count_has2leptons_muel = 0.;
    float tmp_count_has2leptons_mumu = 0.;
    float tmp_count_has2leptons_1llmetjj = 0.;
    float tmp_count_has2leptons_elel_1llmetjj = 0.;
    float tmp_count_has2leptons_elmu_1llmetjj = 0.;
    float tmp_count_has2leptons_muel_1llmetjj = 0.;
    float tmp_count_has2leptons_mumu_1llmetjj = 0.;
    float tmp_count_has2leptons_1llmetjj_2btagM = 0.;
    float tmp_count_has2leptons_elel_1llmetjj_2btagM = 0.;
    float tmp_count_has2leptons_elmu_1llmetjj_2btagM = 0.;
    float tmp_count_has2leptons_muel_1llmetjj_2btagM = 0.;
    float tmp_count_has2leptons_mumu_1llmetjj_2btagM = 0.;

    // ***** ***** *****
    // Trigger Matching
    // ***** ***** *****

    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    // the actual trigger matching to dilepton HLT paths happens only once we have a dilepton candidate to consider in the event

    // ********** 
    // Leptons and dileptons
    // ********** 
    const ElectronsProducer& allelectrons = producers.get<ElectronsProducer>(m_electrons_producer);
    const MuonsProducer& allmuons = producers.get<MuonsProducer>(m_muons_producer);

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
            ele.id_M = allelectrons.ids[ielectron][m_electron_medium_wp_name];
            ele.id_T = allelectrons.ids[ielectron][m_electron_tight_wp_name];
            ele.id_HWW = ele.id_T;
            ele.iso_L = allelectrons.isEB[ielectron] ? (allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut_EB_Loose) : (allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut_EE_Loose);
            ele.iso_T = allelectrons.isEB[ielectron] ? (allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut_EB_Tight) : (allelectrons.relativeIsoR03_withEA[ielectron] < m_electronIsoCut_EE_Tight);
            ele.iso_HWW = ele.iso_T; // FIXME
            ele.gen_matched = allelectrons.matched[ielectron];
            ele.gen_p4 = ele.gen_matched ? allelectrons.gen_p4[ielectron] : null_p4;
            ele.gen_DR = ele.gen_matched ? ROOT::Math::VectorUtil::DeltaR(ele.p4, ele.gen_p4): -1.;
            ele.gen_DPtOverPt = ele.gen_matched ? (ele.p4.Pt() - ele.gen_p4.Pt()) / ele.p4.Pt() : -10.;
            ele.hlt_leg1 = false;
            ele.hlt_leg2 = false;
            // some selection
            if (!ele.id_HWW || !ele.iso_HWW)
                continue;
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
            mu.id_M = allmuons.isMedium[imuon];
            mu.id_T = allmuons.isTight[imuon];
            mu.id_HWW = mu.id_M && (mu.p4.Pt() < 20. ? fabs(allmuons.dxy[imuon]) < 0.01 : fabs(allmuons.dxy[imuon]) < 0.02) && (fabs(allmuons.dz[imuon]) < 0.1);
            mu.iso_L = allmuons.relativeIsoR04_deltaBeta[imuon] < m_muonLooseIsoCut;
            mu.iso_T = allmuons.relativeIsoR04_deltaBeta[imuon] < m_muonTightIsoCut;
            mu.iso_HWW = mu.iso_T; // For the isolation use relative PF isolation (cone size = 0.4) with deltaBeta PU corrections (a.l.a. Run I) and WP < 0.15
            mu.gen_matched = allmuons.matched[imuon];
            mu.gen_p4 = mu.gen_matched ? allmuons.gen_p4[imuon] : null_p4;
            mu.gen_DR = mu.gen_matched ? ROOT::Math::VectorUtil::DeltaR(mu.p4, mu.gen_p4) : -1.;
            mu.gen_DPtOverPt = mu.gen_matched ? (mu.p4.Pt() - mu.gen_p4.Pt()) / mu.p4.Pt() : -10.;
            mu.hlt_leg1 = false;
            mu.hlt_leg2 = false;
            // some selection
            if (!mu.id_HWW || !mu.iso_HWW)
                continue;
            leptons.push_back(mu);
        }
    }//end of loop on muons

    // sort leptons by pt (ignoring flavour, id and iso)
    std::sort(leptons.begin(), leptons.end(), [](const HH::Lepton& lep1, const HH::Lepton& lep2) { return lep1.p4.Pt() > lep2.p4.Pt(); });     

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
            dilep.isPlusMinus = leptons[ilep1].charge > 0 && leptons[ilep2].charge < 0;
            dilep.isMinusPlus = leptons[ilep1].charge < 0 && leptons[ilep2].charge > 0;
            dilep.isMuMu = leptons[ilep1].isMu && leptons[ilep2].isMu;
            dilep.isElEl = leptons[ilep1].isEl && leptons[ilep2].isEl;
            dilep.isElMu = leptons[ilep1].isEl && leptons[ilep2].isMu;
            dilep.isMuEl = leptons[ilep1].isMu && leptons[ilep2].isEl;
            dilep.isSF = dilep.isMuMu || dilep.isElEl;
            dilep.id_LL = leptons[ilep1].id_L && leptons[ilep2].id_L;
            dilep.id_LM = (leptons[ilep1].id_L && leptons[ilep2].id_M) || (leptons[ilep2].id_L && leptons[ilep1].id_M);
            dilep.id_LT = (leptons[ilep1].id_L && leptons[ilep2].id_T) || (leptons[ilep2].id_L && leptons[ilep1].id_T);
            dilep.id_LHWW = (leptons[ilep1].id_L && leptons[ilep2].id_HWW) || (leptons[ilep2].id_L && leptons[ilep1].id_HWW);
            dilep.id_ML = (leptons[ilep1].id_M && leptons[ilep2].id_L) || (leptons[ilep2].id_M && leptons[ilep1].id_L);
            dilep.id_MM = leptons[ilep1].id_M && leptons[ilep2].id_M;
            dilep.id_MT = (leptons[ilep1].id_T && leptons[ilep2].id_M) || (leptons[ilep2].id_T && leptons[ilep1].id_M);
            dilep.id_MHWW = (leptons[ilep1].id_M && leptons[ilep2].id_HWW) || (leptons[ilep2].id_M && leptons[ilep1].id_HWW);
            dilep.id_TL = (leptons[ilep1].id_T && leptons[ilep2].id_L) || (leptons[ilep2].id_T && leptons[ilep1].id_L);
            dilep.id_TM = (leptons[ilep1].id_T && leptons[ilep2].id_M) || (leptons[ilep2].id_T && leptons[ilep1].id_M);
            dilep.id_TT = leptons[ilep1].id_T && leptons[ilep2].id_T;
            dilep.id_THWW = (leptons[ilep1].id_T && leptons[ilep2].id_HWW) || (leptons[ilep2].id_T && leptons[ilep1].id_HWW);
            dilep.id_HWWL = (leptons[ilep1].id_HWW && leptons[ilep2].id_L) || (leptons[ilep2].id_HWW && leptons[ilep1].id_L);
            dilep.id_HWWM = (leptons[ilep1].id_HWW && leptons[ilep2].id_M) || (leptons[ilep2].id_HWW && leptons[ilep1].id_M);
            dilep.id_HWWT = (leptons[ilep1].id_HWW && leptons[ilep2].id_T) || (leptons[ilep2].id_HWW && leptons[ilep1].id_T);
            dilep.id_HWWHWW = leptons[ilep1].id_HWW && leptons[ilep2].id_HWW;
            dilep.iso_LL = leptons[ilep1].iso_L && leptons[ilep2].iso_L;
            dilep.iso_LT = (leptons[ilep1].iso_L && leptons[ilep2].iso_T) || (leptons[ilep2].iso_L && leptons[ilep1].iso_T);
            dilep.iso_LHWW = (leptons[ilep1].iso_L && leptons[ilep2].iso_HWW) || (leptons[ilep2].iso_L && leptons[ilep1].iso_HWW);
            dilep.iso_TL = (leptons[ilep1].iso_T && leptons[ilep2].iso_L) || (leptons[ilep2].iso_T && leptons[ilep1].iso_L);
            dilep.iso_TT = leptons[ilep1].iso_T && leptons[ilep2].iso_T;
            dilep.iso_THWW = (leptons[ilep1].iso_T && leptons[ilep2].iso_HWW) || (leptons[ilep2].iso_T && leptons[ilep1].iso_HWW);
            dilep.iso_HWWL = (leptons[ilep1].iso_HWW && leptons[ilep2].iso_L) || (leptons[ilep2].iso_HWW && leptons[ilep1].iso_L);
            dilep.iso_HWWT = (leptons[ilep1].iso_HWW && leptons[ilep2].iso_T) || (leptons[ilep2].iso_HWW && leptons[ilep1].iso_T);
            dilep.iso_HWWHWW = leptons[ilep1].iso_HWW && leptons[ilep2].iso_HWW;
            dilep.DR_l_l = ROOT::Math::VectorUtil::DeltaR(leptons[ilep1].p4, leptons[ilep2].p4);
            dilep.DPhi_l_l = fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ilep1].p4, leptons[ilep2].p4));
            dilep.ht_l_l = leptons[ilep1].p4.Pt() + leptons[ilep2].p4.Pt();
            if (!hlt.paths.empty()) {
                matchOfflineLepton(hlt, dilep);
                dilep.hlt_idxs = std::make_pair(leptons[dilep.ilep1].hlt_idx, leptons[dilep.ilep2].hlt_idx);
            }
            dilep.gen_matched = leptons[ilep1].gen_matched && leptons[ilep2].gen_matched;
            dilep.gen_p4 = dilep.gen_matched ? leptons[ilep1].gen_p4 + leptons[ilep2].gen_p4 : null_p4;
            dilep.gen_DR = dilep.gen_matched ? ROOT::Math::VectorUtil::DeltaR(dilep.p4, dilep.gen_p4) : -1.;
            dilep.gen_DPtOverPt = dilep.gen_matched ? (dilep.p4.Pt() - dilep.gen_p4.Pt()) / dilep.p4.Pt() : -10.;
            if (event.isRealData()) {
               dilep.trigger_efficiency = 1.;
               dilep.trigger_efficiency_downVariated = 1.;
               dilep.trigger_efficiency_downVariated_Arun = 1.;
               dilep.trigger_efficiency_upVariated = 1.;
               dilep.trigger_efficiency_upVariated_Arun = 1.;
            }
            else {
               fillTriggerEfficiencies(leptons[ilep1], leptons[ilep2], dilep);
            }
            // Some selection
            // Note that ID and isolation criteria are in both electron and muon loops
            if (!dilep.isOS)
                continue;

            // Throw event if there is no matched dilepton trigger path (only on data)
            if (event.isRealData()
                && !((leptons[dilep.ilep1].hlt_leg1 && leptons[dilep.ilep2].hlt_leg2)
                || (leptons[dilep.ilep1].hlt_leg2 && leptons[dilep.ilep2].hlt_leg1))) {
                continue;
            }

            // Counters
            tmp_count_has2leptons = event_weight;
            if (dilep.isElEl)
                tmp_count_has2leptons_elel = event_weight;
            if (dilep.isElMu)
                tmp_count_has2leptons_elmu = event_weight;
            if (dilep.isElEl)
                tmp_count_has2leptons_muel = event_weight;
            if (dilep.isElEl)
                tmp_count_has2leptons_mumu = event_weight;

            // Fill
            ll.push_back(dilep); 
        }
    }
    // have the ll collection sorted by ht
    std::sort(ll.begin(), ll.end(), [&](HH::Dilepton& a, HH::Dilepton& b){return a.ht_l_l > b.ht_l_l;});

    // ***** 
    // Adding MET(s)
    // ***** 
    met.clear();
    llmet.clear();
    const METProducer& pf_met = producers.get<METProducer>(m_met_producer);
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
    //const METProducer& nohf_met = producers.get<METProducer>(m_nohf_met_producer);  // so that nohfmet is available in the tree
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
            myllmet.isPlusMinus = ll[ill].isPlusMinus;
            myllmet.isMinusPlus = ll[ill].isMinusPlus;
            myllmet.isMuMu = ll[ill].isMuMu;
            myllmet.isElEl = ll[ill].isElEl;
            myllmet.isElMu = ll[ill].isElMu;
            myllmet.isMuEl = ll[ill].isMuEl;
            myllmet.isSF = ll[ill].isSF;
            myllmet.id_LL = ll[ill].id_LL;
            myllmet.id_LM = ll[ill].id_LM;
            myllmet.id_LT = ll[ill].id_LT;
            myllmet.id_LHWW = ll[ill].id_LHWW;
            myllmet.id_ML = ll[ill].id_ML;
            myllmet.id_MM = ll[ill].id_MM;
            myllmet.id_MT = ll[ill].id_MT;
            myllmet.id_MHWW = ll[ill].id_MHWW;
            myllmet.id_TL = ll[ill].id_TL;
            myllmet.id_TM = ll[ill].id_TM;
            myllmet.id_TT = ll[ill].id_TT;
            myllmet.id_THWW = ll[ill].id_THWW;
            myllmet.id_HWWL = ll[ill].id_HWWL;
            myllmet.id_HWWM = ll[ill].id_HWWM;
            myllmet.id_HWWT = ll[ill].id_HWWT;
            myllmet.id_HWWHWW = ll[ill].id_HWWHWW;
            myllmet.iso_LL = ll[ill].iso_LL;
            myllmet.iso_LT = ll[ill].iso_LT;
            myllmet.iso_LHWW = ll[ill].iso_LHWW;
            myllmet.iso_TL = ll[ill].iso_TL;
            myllmet.iso_TT = ll[ill].iso_TT;
            myllmet.iso_THWW = ll[ill].iso_THWW;
            myllmet.iso_HWWL = ll[ill].iso_HWWL;
            myllmet.iso_HWWT = ll[ill].iso_HWWT;
            myllmet.iso_HWWHWW = ll[ill].iso_HWWHWW;
            myllmet.DR_l_l = ll[ill].DR_l_l;
            myllmet.DPhi_l_l = ll[ill].DPhi_l_l;
            myllmet.ht_l_l = ll[ill].ht_l_l;
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

    // ***** 
    // Jets and dijets 
    // ***** 
    const JetsProducer& alljets = producers.get<JetsProducer>(m_jets_producer);

    for (unsigned int ijet = 0; ijet < alljets.p4.size(); ijet++)
    {
        float correctionFactor = m_applyBJetRegression ? alljets.regPt[ijet] / alljets.p4[ijet].Pt() : 1.;
/*
        std::cout << "m_jets_producer= " << m_jets_producer
            << "\tm_applyBJetRegression= " << m_applyBJetRegression
            << "\talljets.p4[" << ijet << "].Pt()= " << alljets.p4[ijet].Pt()
            << "\talljets.regPt[" << ijet << "]= " << alljets.regPt[ijet]
            << "\tcorrectionFactor= " << correctionFactor
            << std::endl;
*/
        if ((alljets.p4[ijet].Pt() * correctionFactor > m_jetPtCut) 
            && (fabs(alljets.p4[ijet].Eta()) < m_jetEtaCut))
        {
            HH::Jet myjet;
            myjet.p4 = alljets.p4[ijet] * correctionFactor;
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
            myjet.gen_matched_bParton = (std::abs(alljets.partonFlavor[ijet]) == 5);
            myjet.gen_matched_bHadron = (alljets.hadronFlavor[ijet]) == 5;
            myjet.gen_matched = alljets.matched[ijet];
            myjet.gen_p4 = myjet.gen_matched ? alljets.gen_p4[ijet] : null_p4;
            myjet.gen_DR = myjet.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myjet.p4, myjet.gen_p4) : -1.;
            myjet.gen_DPtOverPt = myjet.gen_matched ? (myjet.p4.Pt() - myjet.gen_p4.Pt()) / myjet.p4.Pt() : -10.;
            myjet.gen_b = (alljets.hadronFlavor[ijet]) == 5; // redundant with gen_matched_bHadron defined above
            myjet.gen_c = (alljets.hadronFlavor[ijet]) == 4;
            myjet.gen_l = (alljets.hadronFlavor[ijet]) < 4;
            // Some selection
            if (!myjet.id_L)
                continue;
            bool isThereACloseSelectedLepton = false;
            for (auto& mylepton: leptons) {
                if (ROOT::Math::VectorUtil::DeltaR(myjet.p4, mylepton.p4) < m_minDR_l_j_Cut) {
                    isThereACloseSelectedLepton = true;
                    break;
                }
            }
            if (isThereACloseSelectedLepton)
                continue;
            jets.push_back(myjet);
        }
    }

    jj.clear();
    // Do NOT change the loop logic here: we expect [0] to be made out of the leading jets
    for (unsigned int ijet1 = 0; ijet1 < jets.size(); ijet1++)
    {
        for (unsigned int ijet2 = ijet1 + 1; ijet2 < jets.size(); ijet2++)
        {
            HH::Dijet myjj;
            myjj.p4 = jets[ijet1].p4 + jets[ijet2].p4;
            myjj.idxs = std::make_pair(jets[ijet1].idx, jets[ijet2].idx);
            myjj.ijet1 = ijet1;
            myjj.ijet2 = ijet2;
            myjj.jid_LL = jets[ijet1].id_L && jets[ijet2].id_L;
            myjj.jid_TT = jets[ijet1].id_T && jets[ijet2].id_T;
            myjj.jid_TLVTLV = jets[ijet1].id_TLV && jets[ijet2].id_TLV;
            myjj.btag_LL = jets[ijet1].btag_L && jets[ijet2].btag_L;
            myjj.btag_LM = (jets[ijet1].btag_L && jets[ijet2].btag_M) || (jets[ijet2].btag_L && jets[ijet1].btag_M);
            myjj.btag_LT = (jets[ijet1].btag_L && jets[ijet2].btag_T) || (jets[ijet2].btag_L && jets[ijet1].btag_T);
            myjj.btag_ML = (jets[ijet1].btag_M && jets[ijet2].btag_L) || (jets[ijet2].btag_M && jets[ijet1].btag_L);
            myjj.btag_MM = jets[ijet1].btag_M && jets[ijet2].btag_M;
            myjj.btag_MT = (jets[ijet1].btag_M && jets[ijet2].btag_T) || (jets[ijet2].btag_M && jets[ijet1].btag_T);
            myjj.btag_TL = (jets[ijet1].btag_T && jets[ijet2].btag_L) || (jets[ijet2].btag_T && jets[ijet1].btag_L);
            myjj.btag_TM = (jets[ijet1].btag_T && jets[ijet2].btag_M) || (jets[ijet2].btag_T && jets[ijet1].btag_M);
            myjj.btag_TT = jets[ijet1].btag_T && jets[ijet2].btag_T;
            myjj.sumCSV = jets[ijet1].CSV + jets[ijet2].CSV;
            myjj.sumJP = jets[ijet1].JP + jets[ijet2].JP;
            myjj.DR_j_j = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, jets[ijet2].p4);
            myjj.DPhi_j_j = fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[ijet1].p4, jets[ijet2].p4));
            myjj.ht_j_j = jets[ijet1].p4.Pt() + jets[ijet2].p4.Pt();
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
        }
    }
    // have the jj collection sorted by ht
    std::sort(jj.begin(), jj.end(), [&](HH::Dijet& a, HH::Dijet& b){return a.p4.Pt() > b.p4.Pt();});

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
            myllmetjj.jet1_p4 = jets[ijet1].p4;
            myllmetjj.jet2_p4 = jets[ijet2].p4;
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
            myllmetjj.jid_LL = jj[ijj].jid_LL;
            myllmetjj.jid_TT = jj[ijj].jid_TT;
            myllmetjj.jid_TLVTLV = jj[ijj].jid_TLVTLV;
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
            myllmetjj.ht_j_j = jj[ijj].ht_j_j;
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
            myllmetjj.isPlusMinus = ll[ill].isPlusMinus;
            myllmetjj.isMinusPlus = ll[ill].isMinusPlus;
            myllmetjj.isMuMu = ll[ill].isMuMu;
            myllmetjj.isElEl = ll[ill].isElEl;
            myllmetjj.isElMu = ll[ill].isElMu;
            myllmetjj.isMuEl = ll[ill].isMuEl;
            myllmetjj.isSF = ll[ill].isSF;
            myllmetjj.id_LL = ll[ill].id_LL;
            myllmetjj.id_LM = ll[ill].id_LM;
            myllmetjj.id_LT = ll[ill].id_LT;
            myllmetjj.id_LHWW = ll[ill].id_LHWW;
            myllmetjj.id_ML = ll[ill].id_ML;
            myllmetjj.id_MM = ll[ill].id_MM;
            myllmetjj.id_MT = ll[ill].id_MT;
            myllmetjj.id_MHWW = ll[ill].id_MHWW;
            myllmetjj.id_TL = ll[ill].id_TL;
            myllmetjj.id_TM = ll[ill].id_TM;
            myllmetjj.id_TT = ll[ill].id_TT;
            myllmetjj.id_THWW = ll[ill].id_THWW;
            myllmetjj.id_HWWL = ll[ill].id_HWWL;
            myllmetjj.id_HWWM = ll[ill].id_HWWM;
            myllmetjj.id_HWWT = ll[ill].id_HWWT;
            myllmetjj.id_HWWHWW = ll[ill].id_HWWHWW;
            myllmetjj.iso_LL = ll[ill].iso_LL;
            myllmetjj.iso_LT = ll[ill].iso_LT;
            myllmetjj.iso_LHWW = ll[ill].iso_LHWW;
            myllmetjj.iso_TL = ll[ill].iso_TL;
            myllmetjj.iso_TT = ll[ill].iso_TT;
            myllmetjj.iso_THWW = ll[ill].iso_THWW;
            myllmetjj.iso_HWWL = ll[ill].iso_HWWL;
            myllmetjj.iso_HWWT = ll[ill].iso_HWWT;
            myllmetjj.iso_HWWHWW = ll[ill].iso_HWWHWW;
            myllmetjj.DR_l_l = ll[ill].DR_l_l;
            myllmetjj.DPhi_l_l = ll[ill].DPhi_l_l;
            myllmetjj.ht_l_l = ll[ill].ht_l_l;
            myllmetjj.trigger_efficiency = ll[ill].trigger_efficiency;
            myllmetjj.trigger_efficiency_downVariated = ll[ill].trigger_efficiency_downVariated;
            myllmetjj.trigger_efficiency_downVariated_Arun = ll[ill].trigger_efficiency_downVariated_Arun;
            myllmetjj.trigger_efficiency_upVariated = ll[ill].trigger_efficiency_upVariated;
            myllmetjj.trigger_efficiency_upVariated_Arun = ll[ill].trigger_efficiency_upVariated_Arun;
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
            myllmetjj.melaAngles = getMELAAngles(llmet[illmet].p4, jj[ijj].p4, leptons[ilep1].p4, leptons[ilep2].p4, jets[ijet1].p4, jets[ijet2].p4);
            myllmetjj.visMelaAngles = getMELAAngles(ll[ill].p4, jj[ijj].p4, leptons[ilep1].p4, leptons[ilep2].p4, jets[ijet1].p4, jets[ijet2].p4); // only take the visible part of the H(ww) candidate

            // Compute MT2. See https://arxiv.org/pdf/1309.6318v1.pdf and https://arxiv.org/pdf/1411.4312v5.pdf
            double px_invisible = myllmetjj.lep1_p4.px() + myllmetjj.lep2_p4.px() + myllmetjj.met_p4.px();
            double py_invisible = myllmetjj.lep1_p4.py() + myllmetjj.lep2_p4.py() + myllmetjj.met_p4.py();

            myllmetjj.MT2 = asymm_mt2_lester_bisect::get_mT2(
                    myllmetjj.jet1_p4.M(), myllmetjj.jet1_p4.px(), myllmetjj.jet1_p4.py(),
                    myllmetjj.jet2_p4.M(), myllmetjj.jet2_p4.px(), myllmetjj.jet2_p4.py(),
                    px_invisible, py_invisible,
                    myllmetjj.lep1_p4.M(), myllmetjj.lep2_p4.M(),
                    0.5 // Absolute precision
                    );


            // Counters
            tmp_count_has2leptons_1llmetjj = event_weight;
            if (myllmetjj.isElEl)
                tmp_count_has2leptons_elel_1llmetjj = event_weight;
            if (myllmetjj.isElMu)
                tmp_count_has2leptons_elmu_1llmetjj = event_weight;
            if (myllmetjj.isElEl)
                tmp_count_has2leptons_muel_1llmetjj = event_weight;
            if (myllmetjj.isElEl)
                tmp_count_has2leptons_mumu_1llmetjj = event_weight;
            if (myllmetjj.btag_MM)
            {
                tmp_count_has2leptons_1llmetjj_2btagM = event_weight;
                if (myllmetjj.isElEl)
                    tmp_count_has2leptons_elel_1llmetjj_2btagM = event_weight;
                if (myllmetjj.isElMu)
                    tmp_count_has2leptons_elmu_1llmetjj_2btagM = event_weight;
                if (myllmetjj.isElEl)
                    tmp_count_has2leptons_muel_1llmetjj_2btagM = event_weight;
                if (myllmetjj.isElEl)
                    tmp_count_has2leptons_mumu_1llmetjj_2btagM = event_weight;
            }
            // Fill
            llmetjj.push_back(myllmetjj);
        }
    }



    // Sort the collections
//    llmetjj_ht.clear(); llmetjj_ht = llmetjj;
//    llmetjj_pt.clear(); llmetjj_pt = llmetjj;
//    llmetjj_mh.clear(); llmetjj_mh = llmetjj;
    llmetjj_csv.clear(); llmetjj_csv = llmetjj;
//    llmetjj_ptOverM.clear(); llmetjj_ptOverM = llmetjj;
//    std::sort(llmetjj_ht.begin(), llmetjj_ht.end(), [&](HH::DileptonMetDijet& a, const HH::DileptonMetDijet& b){return (a.jet1_p4.Pt() + a.jet2_p4.Pt()) > (b.jet1_p4.Pt() + bjet2_p4.Pt());});
//    std::sort(llmetjj_pt.begin(), llmetjj_pt.end(), [&](HH::DileptonMetDijet& a, const HH::DileptonMetDijet& b){return a.jj_p4.Pt() > b.jj_p4.Pt();});
//    std::sort(llmetjj_mh.begin(), llmetjj_mh.end(), [&](HH::DileptonMetDijet& a, const HH::DileptonMetDijet& b){return fabs(a.jj_p4.M() - mh) > fabs(b.jj_p4.M() - mh);});
    std::sort(llmetjj_csv.begin(), llmetjj_csv.end(), [&](HH::DileptonMetDijet& a, const HH::DileptonMetDijet& b){return a.sumCSV > b.sumCSV;});
//    std::sort(llmetjj_ptOverM.begin(), llmetjj_ptOverM.end(), [&](HH::DileptonMetDijet& a, const HH::DileptonMetDijet& b){return (a.jj_p4.Pt() / a.jj_p4.M()) > (b.jj_p4.Pt() / b.jj_p4.M());});

    // Adding some few custom candidates, for convenience
    for (auto &myllmetjj_csv: llmetjj_csv) {
        if (!myllmetjj_csv.id_HWWHWW) continue;
        if (!myllmetjj_csv.iso_HWWHWW) continue;
        // jetID::L is enforced while filling the jet collection
        llmetjj_HWWleptons_nobtag_csv.push_back(myllmetjj_csv);
        if (myllmetjj_csv.btag_LL)
            llmetjj_HWWleptons_btagL_csv.push_back(myllmetjj_csv);
        if (myllmetjj_csv.btag_MM)
            llmetjj_HWWleptons_btagM_csv.push_back(myllmetjj_csv);
        if (myllmetjj_csv.btag_TT)
            llmetjj_HWWleptons_btagT_csv.push_back(myllmetjj_csv);
        // October 2016: asymmetric btag candidates
        if (myllmetjj_csv.btag_LM || myllmetjj_csv.btag_ML)
            llmetjj_HWWleptons_btagLM_csv.push_back(myllmetjj_csv);
        if (myllmetjj_csv.btag_MT || myllmetjj_csv.btag_TM)
            llmetjj_HWWleptons_btagMT_csv.push_back(myllmetjj_csv);
    }


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

    count_has2leptons += tmp_count_has2leptons;
    count_has2leptons_elel += tmp_count_has2leptons_elel;
    count_has2leptons_elmu += tmp_count_has2leptons_elmu;
    count_has2leptons_muel += tmp_count_has2leptons_muel;
    count_has2leptons_mumu += tmp_count_has2leptons_mumu;
    count_has2leptons_1llmetjj += tmp_count_has2leptons_1llmetjj;
    count_has2leptons_elel_1llmetjj += tmp_count_has2leptons_elel_1llmetjj;
    count_has2leptons_elmu_1llmetjj += tmp_count_has2leptons_elmu_1llmetjj;
    count_has2leptons_muel_1llmetjj += tmp_count_has2leptons_muel_1llmetjj;
    count_has2leptons_mumu_1llmetjj += tmp_count_has2leptons_mumu_1llmetjj;
    count_has2leptons_1llmetjj_2btagM += tmp_count_has2leptons_1llmetjj_2btagM;
    count_has2leptons_elel_1llmetjj_2btagM += tmp_count_has2leptons_elel_1llmetjj_2btagM;
    count_has2leptons_elmu_1llmetjj_2btagM += tmp_count_has2leptons_elmu_1llmetjj_2btagM;
    count_has2leptons_muel_1llmetjj_2btagM += tmp_count_has2leptons_muel_1llmetjj_2btagM;
    count_has2leptons_mumu_1llmetjj_2btagM += tmp_count_has2leptons_mumu_1llmetjj_2btagM;


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

#if HH_GEN_DEBUG
    std::function<void(size_t)> print_mother_chain = [&gp, &print_mother_chain](size_t p) {

        if (gp.pruned_mothers_index[p].empty()) {
            std::cout << std::endl;
            return;
        }

        size_t index = gp.pruned_mothers_index[p][0];
            std::cout << " <- #" << index << "(" << gp.pruned_pdg_id[index] << ")";
            print_mother_chain(index);
    };
#endif

        std::function<bool(size_t, size_t)> pruned_decays_from = [&pruned_decays_from, &gp](size_t particle_index, size_t mother_index) -> bool {
            // Iterator over all pruned particles to find if the particle `particle_index` has `mother_index` in its decay history
            if (gp.pruned_mothers_index[particle_index].empty())
                return false;

            size_t index = gp.pruned_mothers_index[particle_index][0];

            if (index == mother_index) {
                return true;
            }

            if (pruned_decays_from(index, mother_index))
                return true;

            return false;
        };

        std::function<bool(size_t, size_t)> pruned_decays_from_pdg_id = [&pruned_decays_from_pdg_id, &gp](size_t particle_index, uint64_t pdg_id) -> bool {
            // Iterator over all pruned particles to find if the particle `particle_index` decays from a particle with pdg id == pdg_id
            if (gp.pruned_mothers_index[particle_index].empty())
                return false;

            size_t index = gp.pruned_mothers_index[particle_index][0];

            if (std::abs(gp.pruned_pdg_id[index]) == pdg_id) {
                return true;
            }

            if (pruned_decays_from_pdg_id(index, pdg_id))
                return true;

            return false;
        };

        // Construct signal gen info

        gen_iX = -1;
        gen_iH1 = gen_iH2 = -1;
        gen_iH1_afterFSR = gen_iH2_afterFSR = -1;
        gen_iB = gen_iBbar = -1;
        gen_iB_afterFSR = gen_iBbar_afterFSR = -1;
        gen_iV2 = gen_iV1 = -1;
        gen_iV2_afterFSR = gen_iV1_afterFSR = -1;
        gen_iLminus = gen_iLplus = -1;
        gen_iLminus_afterFSR = gen_iLplus_afterFSR = -1;
        gen_iNu1 = gen_iNu2 = -1;

        for (unsigned int ip = 0; ip < gp.pruned_p4.size(); ip++) {
            std::bitset<15> flags (gp.pruned_status_flags[ip]);

            if (!flags.test(8))
                continue;

            int64_t pdg_id = gp.pruned_pdg_id[ip];

#if HH_GEN_DEBUG
            std::cout << "[" << ip << "] pdg id: " << pdg_id << "  flags: " << flags << "  p = " << gp.pruned_p4[ip] << std::endl;
            print_mother_chain(ip);
#endif

            auto p4 = gp.pruned_p4[ip];

            if (std::abs(pdg_id) == 35 || std::abs(pdg_id) == 39) {
                ASSIGN_HH_GEN_INFO_NO_FSR(X, "X");
            } else if (pdg_id == 25) {
                ASSIGN_HH_GEN_INFO_2(H1, H2, "Higgs");
            }

            // Only look for Higgs decays if we have found the two Higgs
            if ((gen_iH1 == -1) || (gen_iH2 == -1))
                continue;

            // And if the particle actually come directly from a Higgs
            bool from_h1_decay = pruned_decays_from(ip, gen_iH1);
            bool from_h2_decay = pruned_decays_from(ip, gen_iH2);

            // Only keep particles coming from the tops decay
            if (! from_h1_decay && ! from_h2_decay)
                continue;

            if (pdg_id == 5) {
                ASSIGN_HH_GEN_INFO(B, "B");
            } else if (pdg_id == -5) {
                ASSIGN_HH_GEN_INFO(Bbar, "Bbar");
            }

            // Ignore B decays
            if (pruned_decays_from_pdg_id(ip, 5))
                continue;

            if ((pdg_id == 11) || (pdg_id == 13) || (pdg_id == 15)) {
                ASSIGN_HH_GEN_INFO(Lminus, "L-");
            } else if ((pdg_id == -11) || (pdg_id == -13) || (pdg_id == -15)) {
                ASSIGN_HH_GEN_INFO(Lplus, "L+");
            } else if ((pdg_id == 23) || (std::abs(pdg_id) == 24)) {
                ASSIGN_HH_GEN_INFO_2(V1, V2, "W/Z bosons");
            } else if ((std::abs(pdg_id) == 12) || (std::abs(pdg_id) == 14) || (std::abs(pdg_id) == 16)) {
                ASSIGN_HH_GEN_INFO_2_NO_FSR(Nu1, Nu2, "neutrinos");
            }
        }

        // Swap neutrinos if needed
        if ((gen_iNu1 != -1) && (gen_iNu2 != -1)) {
            if (gp.pruned_pdg_id[gen_iNu1] > 0) {
                std::swap(gen_iNu1, gen_iNu2);
                std::swap(gen_Nu1, gen_Nu2);
            }
        }

        if ((gen_iH1 != -1) && (gen_iH2 != -1)) {
            gen_mHH = (gen_H1 + gen_H2).M();
            gen_costhetastar = getCosThetaStar_CS(gen_H1, gen_H2);
        }

#if HH_GEN_DEBUG
        PRINT_PARTICULE(X);
        PRINT_RESONANCE(H1, H2);
        PRINT_RESONANCE(B, Bbar);
        PRINT_RESONANCE(V1, V2);
        PRINT_RESONANCE(Lminus, Lplus);
        PRINT_RESONANCE_NO_FSR(Nu1, Nu2);

        // Rebuild resonances for consistency checks
        auto LminusNu1 = gen_Lminus + gen_Nu1;
        std::cout << "    gen_(L- Nu1).M() = " << LminusNu1.M() << std::endl;

        auto LplusNu2 = gen_Lplus + gen_Nu2;
        std::cout << "    gen_(L+ Nu2).M() = " << LplusNu2.M() << std::endl;

        auto LminusNu1_afterFSR = gen_Lminus_afterFSR + gen_Nu1;
        std::cout << "    gen_(L- Nu1)_afterFSR.M() = " << LminusNu1_afterFSR.M() << std::endl;

        auto LplusNu2_afterFSR = gen_Lplus_afterFSR + gen_Nu2;
        std::cout << "    gen_(L+ Nu2)_afterFSR.M() = " << LplusNu2_afterFSR.M() << std::endl;
        
        auto LLNuNu = gen_Lplus + gen_Lminus + gen_Nu1 + gen_Nu2;
        std::cout << "    gen_(LL NuNu).M() = " << LLNuNu.M() << std::endl;

        auto LLNuNu_afterFSR = gen_Lplus_afterFSR + gen_Lminus_afterFSR + gen_Nu1 + gen_Nu2;
        std::cout << "    gen_(LL NuNu)_afterFSR.M() = " << LLNuNu_afterFSR.M() << std::endl;

        auto LLNuNuBB = gen_Lplus + gen_Lminus + gen_Nu1 + gen_Nu2 + gen_B + gen_Bbar;
        std::cout << "    gen_(LL NuNu BB).M() = " << LLNuNuBB.M() << std::endl;

        auto LLNuNuBB_afterFSR = gen_Lplus_afterFSR + gen_Lminus_afterFSR + gen_Nu1 + gen_Nu2 + gen_B_afterFSR + gen_Bbar_afterFSR;
        std::cout << "    gen_(LL NuNu BB)_afterFSR.M() = " << LLNuNuBB_afterFSR.M() << std::endl;
#endif

        // ***** ***** *****
        // Matching
        // ***** ***** *****
        BRANCH(gen_deltaR_jet_B, std::vector<float>);    
        BRANCH(gen_deltaR_jet_Bbar, std::vector<float>);    
        BRANCH(gen_deltaR_jet_B_afterFSR, std::vector<float>);    
        BRANCH(gen_deltaR_jet_Bbar_afterFSR, std::vector<float>);    
        BRANCH(gen_deltaR_electron_L1, std::vector<float>);    
        BRANCH(gen_deltaR_electron_L2, std::vector<float>);    
        BRANCH(gen_deltaR_electron_L1_afterFSR, std::vector<float>);    
        BRANCH(gen_deltaR_electron_L2_afterFSR, std::vector<float>);    
        BRANCH(gen_deltaR_muon_L1, std::vector<float>);    
        BRANCH(gen_deltaR_muon_L2, std::vector<float>);    
        BRANCH(gen_deltaR_muon_L1_afterFSR, std::vector<float>);    
        BRANCH(gen_deltaR_muon_L2_afterFSR, std::vector<float>);    
    
        for (auto p4: alljets.gen_p4) {
            gen_deltaR_jet_B.push_back(deltaR(p4, gen_B));
            gen_deltaR_jet_Bbar.push_back(deltaR(p4, gen_Bbar));
            gen_deltaR_jet_B_afterFSR.push_back(deltaR(p4, gen_B_afterFSR));
            gen_deltaR_jet_Bbar_afterFSR.push_back(deltaR(p4, gen_Bbar_afterFSR));
        }
        for (auto p4: allelectrons.gen_p4) {
            gen_deltaR_electron_L1.push_back(deltaR(p4, gen_Lminus));
            gen_deltaR_electron_L2.push_back(deltaR(p4, gen_Lplus));
            gen_deltaR_electron_L1_afterFSR.push_back(deltaR(p4, gen_Lminus_afterFSR));
            gen_deltaR_electron_L2_afterFSR.push_back(deltaR(p4, gen_Lplus_afterFSR));
        }
        for (auto p4: allmuons.gen_p4) {
            gen_deltaR_muon_L1.push_back(deltaR(p4, gen_Lminus));
            gen_deltaR_muon_L2.push_back(deltaR(p4, gen_Lplus));
            gen_deltaR_muon_L1_afterFSR.push_back(deltaR(p4, gen_Lminus_afterFSR));
            gen_deltaR_muon_L2_afterFSR.push_back(deltaR(p4, gen_Lplus_afterFSR));
        }

    // TTBAR MC TRUTH
    const GenParticlesProducer& gen_particles = gp;

    // 'Pruned' particles are from the hard process
    // 'Packed' particles are stable particles

#if TT_GEN_DEBUG
    std::function<void(size_t)> print_mother_chain = [&gen_particles, &print_mother_chain](size_t p) {

        if (gen_particles.pruned_mothers_index[p].empty()) {
            std::cout << std::endl;
            return;
        }

        size_t index = gen_particles.pruned_mothers_index[p][0];
            std::cout << " <- #" << index << "(" << gen_particles.pruned_pdg_id[index] << ")";
            print_mother_chain(index);
    };
#endif

#define ASSIGN_INDEX( X ) \
    if (flags.isLastCopy()) { \
        gen_##X = i; \
    }\
    if (flags.isFirstCopy()) { \
        gen_##X##_beforeFSR = i; \
    }

// Assign index to X if it's empty, or Y if not
#define ASSIGN_INDEX2(X, Y, ERROR) \
    if (flags.isLastCopy()) { \
        if (gen_##X == 0) \
            gen_##X = i; \
        else if (gen_##Y == 0)\
            gen_##Y = i; \
        else \
            std::cout << ERROR << std::endl; \
    } \
    if (flags.isFirstCopy()) { \
        if (gen_##X##_beforeFSR == 0) \
            gen_##X##_beforeFSR = i; \
        else if (gen_##Y##_beforeFSR == 0)\
            gen_##Y##_beforeFSR = i; \
        else \
            std::cout << ERROR << std::endl; \
    }

    gen_t = 0; // Index of the top quark
    gen_t_beforeFSR = 0; // Index of the top quark, before any FSR
    gen_tbar = 0; // Index of the anti-top quark
    gen_tbar_beforeFSR = 0; // Index of the anti-top quark, before any FSR

    gen_b = 0; // Index of the b quark coming from the top decay
    gen_b_beforeFSR = 0; // Index of the b quark coming from the top decay, before any FSR
    gen_bbar = 0; // Index of the anti-b quark coming from the anti-top decay
    gen_bbar_beforeFSR = 0; // Index of the anti-b quark coming from the anti-top decay, before any FSR

    gen_jet1_t = 0; // Index of the first jet from the top decay chain
    gen_jet1_t_beforeFSR = 0; // Index of the first jet from the top decay chain, before any FSR
    gen_jet2_t = 0; // Index of the second jet from the top decay chain
    gen_jet2_t_beforeFSR = 0; // Index of the second jet from the top decay chain, before any FSR

    gen_jet1_tbar = 0; // Index of the first jet from the anti-top decay chain
    gen_jet1_tbar_beforeFSR = 0; // Index of the first jet from the anti-top decay chain, before any FSR
    gen_jet2_tbar = 0; // Index of the second jet from the anti-top decay chain
    gen_jet2_tbar_beforeFSR = 0; // Index of the second jet from the anti-top decay chain, before any FSR

    gen_lepton_t = 0; // Index of the lepton from the top decay chain
    gen_lepton_t_beforeFSR = 0; // Index of the lepton from the top decay chain, before any FSR
    gen_neutrino_t = 0; // Index of the neutrino from the top decay chain
    gen_neutrino_t_beforeFSR = 0; // Index of the neutrino from the top decay chain, before any FSR

    gen_lepton_tbar = 0; // Index of the lepton from the anti-top decay chain
    gen_lepton_tbar_beforeFSR = 0; // Index of the lepton from the anti-top decay chain, before any FSR
    gen_neutrino_tbar = 0; // Index of the neutrino from the anti-top decay chain
    gen_neutrino_tbar_beforeFSR = 0; // Index of the neutrino from the anti-top decay chain, before any FSR
    for (size_t i = 0; i < gen_particles.pruned_pdg_id.size(); i++) {

        int16_t pdg_id = gen_particles.pruned_pdg_id[i];
        uint16_t a_pdg_id = std::abs(pdg_id);

        // We only care of particles with PDG id <= 16 (16 is neutrino tau)
        if (a_pdg_id > 16)
            continue;

        GenStatusFlags flags(gen_particles.pruned_status_flags[i]);

        if (! flags.isLastCopy() && ! flags.isFirstCopy())
            continue;

        if (! flags.fromHardProcess())
            continue;

#if TT_GEN_DEBUG
        std::cout << "---" << std::endl;
        std::cout << "Gen particle #" << i << ": PDG id: " << gen_particles.pruned_pdg_id[i];
        print_mother_chain(i);
        flags.dump();
#endif

        if (pdg_id == 6) {
            ASSIGN_INDEX(t);
            continue;
        } else if (pdg_id == -6) {
            ASSIGN_INDEX(tbar);
            continue;
        }

        if (gen_t == 0 || gen_tbar == 0) {
            // Don't bother if we don't have found the tops
            continue;
        }

        bool from_t_decay = pruned_decays_from(i, gen_t);
        bool from_tbar_decay = pruned_decays_from(i, gen_tbar);

        // Only keep particles coming from the tops decay
        if (! from_t_decay && ! from_tbar_decay)
            continue;

        if (pdg_id == 5) {
            // Maybe it's a b coming from the W decay
            if (!flags.isFirstCopy() && flags.isLastCopy() && gen_b == 0) {

                // This can be a B decaying from a W
                // However, we can't rely on the presence of the W in the decay chain, as it may be generator specific
                // Since it's the last copy (ie, after FSR), we can check if this B comes from the B assigned to the W decay (ie, gen_jet1_t_beforeFSR, gen_jet2_t_beforeFSR)
                // If yes, then it's not the B coming directly from the top decay
                if ((gen_jet1_t_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet1_t_beforeFSR]) == 5) ||
                    (gen_jet2_t_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet2_t_beforeFSR]) == 5) ||
                    (gen_jet1_tbar_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet1_tbar_beforeFSR]) == 5) ||
                    (gen_jet2_tbar_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet2_tbar_beforeFSR]) == 5)) {

#if TT_GEN_DEBUG
                    std::cout << "A quark coming from W decay is a b" << std::endl;
#endif

                    if (! (gen_jet1_tbar_beforeFSR != 0 && pruned_decays_from(i, gen_jet1_tbar_beforeFSR)) &&
                        ! (gen_jet2_tbar_beforeFSR != 0 && pruned_decays_from(i, gen_jet2_tbar_beforeFSR)) &&
                        ! (gen_jet1_t_beforeFSR != 0 && pruned_decays_from(i, gen_jet1_t_beforeFSR)) &&
                        ! (gen_jet2_t_beforeFSR != 0 && pruned_decays_from(i, gen_jet2_t_beforeFSR))) {
#if TT_GEN_DEBUG
                        std::cout << "This after-FSR b quark is not coming from a W decay" << std::endl;
#endif
                        gen_b = i;
                        continue;
                    }
#if TT_GEN_DEBUG
                    else {
                        std::cout << "This after-FSR b quark comes from a W decay" << std::endl;
                    }
#endif
                } else {
#if TT_GEN_DEBUG
                    std::cout << "Assigning gen_b" << std::endl;
#endif
                    gen_b = i;
                    continue;
                }
            } else if (flags.isFirstCopy() && gen_b_beforeFSR == 0) {
                gen_b_beforeFSR = i;
                continue;
            } else {
#if TT_GEN_DEBUG
                std::cout << "This should not happen!" << std::endl;
#endif
            }
        } else if (pdg_id == -5) {
            if (!flags.isFirstCopy() && flags.isLastCopy() && gen_bbar == 0) {

                // This can be a B decaying from a W
                // However, we can't rely on the presence of the W in the decay chain, as it may be generator specific
                // Since it's the last copy (ie, after FSR), we can check if this B comes from the B assigned to the W decay (ie, gen_jet1_t_beforeFSR, gen_jet2_t_beforeFSR)
                // If yes, then it's not the B coming directly from the top decay
                if ((gen_jet1_t_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet1_t_beforeFSR]) == 5) ||
                    (gen_jet2_t_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet2_t_beforeFSR]) == 5) ||
                    (gen_jet1_tbar_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet1_tbar_beforeFSR]) == 5) ||
                    (gen_jet2_tbar_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet2_tbar_beforeFSR]) == 5)) {

#if TT_GEN_DEBUG
                    std::cout << "A quark coming from W decay is a bbar" << std::endl;
#endif

                    if (! (gen_jet1_tbar_beforeFSR != 0 && pruned_decays_from(i, gen_jet1_tbar_beforeFSR)) &&
                        ! (gen_jet2_tbar_beforeFSR != 0 && pruned_decays_from(i, gen_jet2_tbar_beforeFSR)) &&
                        ! (gen_jet1_t_beforeFSR != 0 && pruned_decays_from(i, gen_jet1_t_beforeFSR)) &&
                        ! (gen_jet2_t_beforeFSR != 0 && pruned_decays_from(i, gen_jet2_t_beforeFSR))) {
#if TT_GEN_DEBUG
                        std::cout << "This after-fsr b anti-quark is not coming from a W decay" << std::endl;
#endif
                        gen_bbar = i;
                        continue;
                    }
#if TT_GEN_DEBUG
                    else {
                        std::cout << "This after-fsr b anti-quark comes from a W decay" << std::endl;
                    }
#endif
                } else {
#if TT_GEN_DEBUG
                    std::cout << "Assigning gen_bbar" << std::endl;
#endif
                    gen_bbar = i;
                    continue;
                }
            } else if (flags.isFirstCopy() && gen_bbar_beforeFSR == 0) {
                gen_bbar_beforeFSR = i;
                continue;
            }
        }

        if ((gen_tbar == 0) || (gen_t == 0))
            continue;

        if (gen_t != 0 && from_t_decay) {
#if TT_GEN_DEBUG
        std::cout << "Coming from the top chain decay" << std::endl;
#endif
            if (a_pdg_id >= 1 && a_pdg_id <= 5) {
                ASSIGN_INDEX2(jet1_t, jet2_t, "Error: more than two quarks coming from top decay");
            } else if (a_pdg_id == 11 || a_pdg_id == 13 || a_pdg_id == 15) {
                ASSIGN_INDEX(lepton_t);
            } else if (a_pdg_id == 12 || a_pdg_id == 14 || a_pdg_id == 16) {
                ASSIGN_INDEX(neutrino_t);
            } else {
                std::cout << "Error: unknown particle coming from top decay - #" << i << " ; PDG Id: " << pdg_id << std::endl;
            }
        } else if (gen_tbar != 0 && from_tbar_decay) {
#if TT_GEN_DEBUG
        std::cout << "Coming from the anti-top chain decay" << std::endl;
#endif
            if (a_pdg_id >= 1 && a_pdg_id <= 5) {
                ASSIGN_INDEX2(jet1_tbar, jet2_tbar, "Error: more than two quarks coming from anti-top decay");
            } else if (a_pdg_id == 11 || a_pdg_id == 13 || a_pdg_id == 15) {
                ASSIGN_INDEX(lepton_tbar);
            } else if (a_pdg_id == 12 || a_pdg_id == 14 || a_pdg_id == 16) {
                ASSIGN_INDEX(neutrino_tbar);
            } else {
                std::cout << "Error: unknown particle coming from anti-top decay - #" << i << " ; PDG Id: " << pdg_id << std::endl;
            }
        }
    }

    if (!gen_t || !gen_tbar) {
#if TT_GEN_DEBUG
        std::cout << "This is not a ttbar event" << std::endl;
#endif
        gen_ttbar_decay_type = NotTT;
        return;
    }

    if ((gen_jet1_t != 0) && (gen_jet2_t != 0) && (gen_jet1_tbar != 0) && (gen_jet2_tbar != 0)) {
#if TT_GEN_DEBUG
        std::cout << "Hadronic ttbar decay" << std::endl;
#endif
        gen_ttbar_decay_type = Hadronic;
    } else if (
            ((gen_lepton_t != 0) && (gen_lepton_tbar == 0)) ||
            ((gen_lepton_t == 0) && (gen_lepton_tbar != 0))
            ) {

#if TT_GEN_DEBUG
        std::cout << "Semileptonic ttbar decay" << std::endl;
#endif

        uint16_t lepton_pdg_id;
        if (gen_lepton_t != 0)
            lepton_pdg_id = std::abs(gen_particles.pruned_pdg_id[gen_lepton_t]);
        else
            lepton_pdg_id = std::abs(gen_particles.pruned_pdg_id[gen_lepton_tbar]);

        if (lepton_pdg_id == 11)
            gen_ttbar_decay_type = Semileptonic_e;
        else if (lepton_pdg_id == 13)
            gen_ttbar_decay_type = Semileptonic_mu;
        else
            gen_ttbar_decay_type = Semileptonic_tau;
    } else if (gen_lepton_t != 0 && gen_lepton_tbar != 0) {
        uint16_t lepton_t_pdg_id = std::abs(gen_particles.pruned_pdg_id[gen_lepton_t]);
        uint16_t lepton_tbar_pdg_id = std::abs(gen_particles.pruned_pdg_id[gen_lepton_tbar]);

#if TT_GEN_DEBUG
        std::cout << "Dileptonic ttbar decay" << std::endl;
#endif

        if (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 11)
            gen_ttbar_decay_type = Dileptonic_ee;
        else if (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 13)
            gen_ttbar_decay_type = Dileptonic_mumu;
        else if (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 15)
            gen_ttbar_decay_type = Dileptonic_tautau;
        else if (
                (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 13) ||
                (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 11)
                ) {
            gen_ttbar_decay_type = Dileptonic_mue;
        }
        else if (
                (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 15) ||
                (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 11)
                ) {
            gen_ttbar_decay_type = Dileptonic_etau;
        }
        else if (
                (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 15) ||
                (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 13)
                ) {
            gen_ttbar_decay_type = Dileptonic_mutau;
        } else {
            std::cout << "Error: unknown dileptonic ttbar decay." << std::endl;
            gen_ttbar_decay_type = NotTT;
            return;
        }
    } else {
        std::cout << "Error: unknown ttbar decay." << std::endl;
        gen_ttbar_decay_type = UnknownTT;
    }
    } // end of if !event.isRealData()

}

void HHAnalyzer::endJob(MetadataManager& metadata)
{
    metadata.add(this->m_name + "count_has2leptons", count_has2leptons);
    metadata.add(this->m_name + "count_has2leptons_elel", count_has2leptons_elel);
    metadata.add(this->m_name + "count_has2leptons_elmu", count_has2leptons_elmu);
    metadata.add(this->m_name + "count_has2leptons_muel", count_has2leptons_muel);
    metadata.add(this->m_name + "count_has2leptons_mumu", count_has2leptons_mumu);
    metadata.add(this->m_name + "count_has2leptons_1llmetjj", count_has2leptons_1llmetjj);
    metadata.add(this->m_name + "count_has2leptons_elel_1llmetjj", count_has2leptons_elel_1llmetjj);
    metadata.add(this->m_name + "count_has2leptons_elmu_1llmetjj", count_has2leptons_elmu_1llmetjj);
    metadata.add(this->m_name + "count_has2leptons_muel_1llmetjj", count_has2leptons_muel_1llmetjj);
    metadata.add(this->m_name + "count_has2leptons_mumu_1llmetjj", count_has2leptons_mumu_1llmetjj);
    metadata.add(this->m_name + "count_has2leptons_1llmetjj_2btagM", count_has2leptons_1llmetjj_2btagM);
    metadata.add(this->m_name + "count_has2leptons_elel_1llmetjj_2btagM", count_has2leptons_elel_1llmetjj_2btagM);
    metadata.add(this->m_name + "count_has2leptons_elmu_1llmetjj_2btagM", count_has2leptons_elmu_1llmetjj_2btagM);
    metadata.add(this->m_name + "count_has2leptons_muel_1llmetjj_2btagM", count_has2leptons_muel_1llmetjj_2btagM);
    metadata.add(this->m_name + "count_has2leptons_mumu_1llmetjj_2btagM", count_has2leptons_mumu_1llmetjj_2btagM);
}
