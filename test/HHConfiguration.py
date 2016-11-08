
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework
from cp3_llbb.Framework import METProducer
from cp3_llbb.Framework.CmdLine import CmdLine

options = CmdLine()
runOnData = options.runOnData == 1

globalTag_ = '80X_mcRun2_asymptotic_2016_miniAODv2_v1'
processName_ = 'PAT'
if runOnData :
    globalTag_ = '80X_dataRun2_Prompt_ICHEP16JEC_v0'
    processName_ = 'RECO'

framework = Framework.Framework(runOnData, eras.Run2_25ns, globalTag=globalTag_, processName=processName_)

framework.addAnalyzer('hh_analyzer', cms.PSet(
        type = cms.string('hh_analyzer'),
        prefix = cms.string('hh_'),
        enable = cms.bool(True),
        parameters = cms.PSet(
            # Producers
            electronsProducer = cms.string('electrons'),
            muonsProducer = cms.string('muons'),
            jetsProducer = cms.string('jets'),
            metProducer = cms.string('met'),
            nohfMETProducer = cms.string('nohf_met'),
            # Here are the default value (just to show what is configurable)
            # FIXME: isolation has not been updated for 80X
            electronIsoCut_EB_Loose = cms.untracked.double(0.0893), # https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
            electronIsoCut_EE_Loose = cms.untracked.double(0.121),
            electronIsoCut_EB_Tight = cms.untracked.double(0.0354),
            electronIsoCut_EE_Tight = cms.untracked.double(0.0646),
            leadingElectronPtCut = cms.untracked.double(20),
            subleadingElectronPtCut = cms.untracked.double(15),
            electronEtaCut = cms.untracked.double(2.5),
            muonLooseIsoCut = cms.untracked.double(.25), # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO 
            muonTightIsoCut = cms.untracked.double(.15), # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO 
            leadingMuonPtCut = cms.untracked.double(20),
            subleadingMuonPtCut = cms.untracked.double(10),
            muonEtaCut = cms.untracked.double(2.4),
            electrons_loose_wp_name = cms.untracked.string("cutBasedElectronID-Summer16-80X-V1-loose"),
            electrons_medium_wp_name = cms.untracked.string("cutBasedElectronID-Summer16-80X-V1-medium"),
            electrons_tight_wp_name = cms.untracked.string("cutBasedElectronID-Summer16-80X-V1-tight"),
            jetEtaCut = cms.untracked.double(2.4),
            jetPtCut = cms.untracked.double(20),
            # BTAG INFO
            discr_name =  cms.untracked.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
            discr_cut_loose =  cms.untracked.double(0.460),
            discr_cut_medium =  cms.untracked.double(0.800),
            discr_cut_tight =  cms.untracked.double(0.935),
            minDR_l_j_Cut = cms.untracked.double(0.3),
            hltDRCut = cms.untracked.double(0.1),
            hltDPtCut = cms.untracked.double(0.5),  # cut will be DPt/Pt < hltDPtCut
            applyBJetRegression = cms.untracked.bool(False), # BE SURE TO ACTIVATE computeRegression FLAG BELOW

            hlt_efficiencies = cms.untracked.PSet(

                    IsoMu17leg = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Muon_DoubleMu_IsoMu17leg_Run2016_PTvsETA_Run271036to276811_HWW_weighted.json'),
                    IsoMu8orIsoTkMu8leg = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Muon_DoubleMu_IsoMu8orIsoTkMu8leg_Run2016_PTvsETA_Run271036to276811_HWW_weighted.json'),

                    DoubleEleHighPtleg = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Electron_HLT_DoubleEleLegHigPt_HWW.json'),
                    DoubleEleLowPtleg = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Electron_HLT_DoubleEleLegLowPt_HWW.json'),

                    EleMuHighPtleg = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Electron_HLT_EleMuLegHigPt_HWW.json'),
                    MuEleLowPtleg = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Electron_HLT_MuEleLegLowPt_HWW.json'),

                    IsoMu8leg = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Muon_DoubleMu_IsoMu8leg_Run2016_PTvsETA_Run271036to276811_HWW_weighted.json'),
                    IsoMu23leg = cms.untracked.FileInPath('cp3_llbb/Framework/data/Efficiencies/Muon_DoubleMu_IsoMu23leg_Run2016_PTvsETA_Run271036to276811_HWW_weighted.json'),
            )
        )
    )
)
# Add PUPPI MET
puppiCfg = cms.PSet(METProducer.default_configuration.clone())
puppiCfg.prefix = cms.string('puppimet_')
puppiCfg.parameters.met = cms.untracked.InputTag('slimmedMETsPuppi')
framework.addProducer('puppimet', puppiCfg)

# Remove fat jets
framework.removeProducer('fat_jets')

framework.getProducer('hlt').parameters.triggers = cms.untracked.FileInPath('cp3_llbb/HHAnalysis/data/triggers.xml')
#framework.getProducer('jets').parameters.cut = cms.untracked.string("pt > 20")
#framework.getProducer('jets').parameters.computeRegression = cms.untracked.bool(True)

framework.redoJEC()
framework.smearJets()
framework.doSystematics(['jec', 'jer'])
process = framework.create()

if runOnData : 
    process.source.fileNames = cms.untracked.vstring(
        '/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/A6AC80E5-121A-E611-A689-02163E01439E.root'
        )
else : 
    process.source.fileNames = cms.untracked.vstring(
        '/store/mc/RunIISpring16MiniAODv2/GluGluToRadionToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/10000/702A12BA-1B3B-E611-9A14-0025907B4EC0.root'
#        '/store/mc/RunIISpring16MiniAODv2/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/60000/022B4AD3-7A1B-E611-812B-28924A33B9AA.root'
        )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
