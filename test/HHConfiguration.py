
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework
from cp3_llbb.Framework import METProducer
from cp3_llbb.Framework.CmdLine import CmdLine

options = CmdLine()
runOnData = options.runOnData == 1

globalTag_ = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
processName_ = 'PAT'
if runOnData :
    globalTag_ = '80X_dataRun2_2016SeptRepro_v7'
    processName_ = 'RECO'

framework = Framework.Framework(runOnData, eras.Run2_25ns, globalTag=globalTag_, processName=processName_)

framework.addAnalyzer('hh_analyzer', cms.PSet(
        type = cms.string('hh_analyzer'),
        prefix = cms.string('hh_'),
        enable = cms.bool(True),
        categories_parameters = cms.PSet(
            # Per-category lepton pt cuts
            mumu_leadingLeptonPtCut = cms.untracked.double(20), # muon
            mumu_subleadingLeptonPtCut = cms.untracked.double(10), # muon
            elel_leadingLeptonPtCut = cms.untracked.double(25), # electron
            elel_subleadingLeptonPtCut = cms.untracked.double(15), # electron
            muel_leadingLeptonPtCut = cms.untracked.double(25), # muon
            muel_subleadingLeptonPtCut = cms.untracked.double(15), # electron
            elmu_leadingLeptonPtCut = cms.untracked.double(25), # electron
            elmu_subleadingLeptonPtCut = cms.untracked.double(10), # muon
        ),
        parameters = cms.PSet(
            # Producers
            electronsProducer = cms.string('electrons'),
            muonsProducer = cms.string('muons'),
            jetsProducer = cms.string('jets'),
            metProducer = cms.string('met'),
            nohfMETProducer = cms.string('nohf_met'),

            # Pre-selection pt cut, applied to all leptons
            leadingElectronPtCut = cms.untracked.double(25),
            subleadingElectronPtCut = cms.untracked.double(15),
            leadingMuonPtCut = cms.untracked.double(20),
            subleadingMuonPtCut = cms.untracked.double(10),

            electronEtaCut = cms.untracked.double(2.5),
            muonLooseIsoCut = cms.untracked.double(.25), # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO 
            muonTightIsoCut = cms.untracked.double(.15), # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO 
            muonEtaCut = cms.untracked.double(2.4),
            electrons_loose_wp_name = cms.untracked.string("cutBasedElectronID-Summer16-80X-V1-loose"),
            electrons_medium_wp_name = cms.untracked.string("cutBasedElectronID-Summer16-80X-V1-medium"),
            electrons_tight_wp_name = cms.untracked.string("cutBasedElectronID-Summer16-80X-V1-tight"),
            electrons_hlt_safe_wp_name = cms.untracked.string("cutBasedElectronHLTPreselection-Summer16-V1"),
            jetEtaCut = cms.untracked.double(2.4),
            jetPtCut = cms.untracked.double(20),

            # BTAG INFO
            # Working points from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
            discr_name =  cms.untracked.string("pfCombinedMVAV2BJetTags"),
            discr_cut_loose =  cms.untracked.double(-0.5884),
            discr_cut_medium =  cms.untracked.double(0.4432),
            discr_cut_tight =  cms.untracked.double(0.9432),

            minDR_l_j_Cut = cms.untracked.double(0.3),
            hltDRCut = cms.untracked.double(0.1),
            hltDPtCut = cms.untracked.double(0.5),  # cut will be DPt/Pt < hltDPtCut
            applyBJetRegression = cms.untracked.bool(False), # BE SURE TO ACTIVATE computeRegression FLAG BELOW

            hlt_efficiencies = cms.untracked.PSet(

                    IsoMu17leg = cms.untracked.FileInPath('cp3_llbb/HHAnalysis/data/Efficiencies/Muon_DoubleIsoMu17Mu8_IsoMu17leg.json'),
                    IsoMu8orIsoTkMu8leg = cms.untracked.FileInPath('cp3_llbb/HHAnalysis/data/Efficiencies/Muon_DoubleIsoMu17TkMu8_IsoMu8legORTkMu8leg.json'),

                    DoubleEleHighPtleg = cms.untracked.FileInPath('cp3_llbb/HHAnalysis/data/Efficiencies/Electron_IsoEle23Leg.json'),
                    DoubleEleLowPtleg = cms.untracked.FileInPath('cp3_llbb/HHAnalysis/data/Efficiencies/Electron_IsoEle12Leg.json'),

                    EleMuHighPtleg = cms.untracked.FileInPath('cp3_llbb/HHAnalysis/data/Efficiencies/Electron_IsoEle23Leg.json'),
                    MuEleLowPtleg = cms.untracked.FileInPath('cp3_llbb/HHAnalysis/data/Efficiencies/Electron_IsoEle12Leg.json'),

                    IsoMu8leg = cms.untracked.FileInPath('cp3_llbb/HHAnalysis/data/Efficiencies/Muon_XPathIsoMu8leg.json'),
                    IsoMu23leg = cms.untracked.FileInPath('cp3_llbb/HHAnalysis/data/Efficiencies/Muon_XPathIsoMu23leg.json'),
            )
        )
    )
)

# Remove fat jets
framework.removeProducer('fat_jets')

framework.getProducer('hlt').parameters.triggers = cms.untracked.FileInPath('cp3_llbb/HHAnalysis/data/triggers.xml')
# framework.getProducer('jets').parameters.cut = cms.untracked.string("pt > 20")
#framework.getProducer('jets').parameters.computeRegression = cms.untracked.bool(True)

if runOnData:
    framework.redoJEC()

framework.applyMuonCorrection('rochester')

framework.applyElectronRegression()
framework.applyElectronSmearing()

if not runOnData:
    framework.smearJets(resolutionFile='cp3_llbb/Framework/data/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt', scaleFactorFile='cp3_llbb/Framework/data/Spring16_25nsV10_MC_SF_AK4PFchs.txt')
    framework.doSystematics(['jec', 'jer'], jec={'uncertaintiesFile': 'cp3_llbb/HHAnalysis/data/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt', 'splitBySources': True})

process = framework.create()

process.framework.treeFlushSize = cms.untracked.uint64(5 * 1024 * 1024)

if runOnData : 
    process.source.fileNames = cms.untracked.vstring(
            '/store/data/Run2016F/DoubleMuon/MINIAOD/23Sep2016-v1/50000/040EDEBA-0490-E611-A424-008CFA110C68.root'
        )
else : 
    process.source.fileNames = cms.untracked.vstring(
            # Signal
            '/store/mc/RunIISummer16MiniAODv2/GluGluToHHTo2B2VTo2L2Nu_node_SM_13TeV-madgraph-v2/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/2E1015E2-71D9-E611-911E-02163E019E19.root'

            # TT
            # '/store/mc/RunIISummer16MiniAODv2/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/00ED79D3-CFC1-E611-B748-3417EBE64483.root'

            # DY
            # '/store/mc/RunIISummer16MiniAODv2/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/00CEFB4F-C1D2-E611-BBF4-7845C4FC3C11.root<Paste>'
        )

#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.source.skipEvents = cms.untracked.uint32(10)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
