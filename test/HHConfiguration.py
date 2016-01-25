
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework
from cp3_llbb.Framework import METProducer

runOnData = False

if runOnData :
    globalTag = '74X_dataRun2_v2'
else : 
    globalTag = '74X_mcRun2_asymptotic_v2'

framework = Framework.Framework(runOnData, eras.Run2_25ns, globalTag=globalTag)

framework.addAnalyzer('hh_analyzer', cms.PSet(
        type = cms.string('hh_analyzer'),
        prefix = cms.string('hh_'),
        enable = cms.bool(True),
        parameters = cms.PSet(
            # Here are the default value (just to show what is configurable)
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
            electrons_loose_wp_name = cms.untracked.string("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
            electrons_tight_wp_name = cms.untracked.string("cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
            jetEtaCut = cms.untracked.double(2.4),
            jetPtCut = cms.untracked.double(20),
            discr_name =  cms.untracked.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
            discr_cut_loose =  cms.untracked.double(0.605),
            discr_cut_medium =  cms.untracked.double(0.89),
            discr_cut_tight =  cms.untracked.double(0.97),
            minDR_l_j_Cut = cms.untracked.double(0.3),
            hltDRCut = cms.untracked.double(0.3),
            hltDPtCut = cms.untracked.double(0.5)  # cut will be DPt/Pt < hltDPtCut
            )
        )
    ) 

#framework.redoJEC()
#framework.smearJets()
framework.doSystematics(['jec']) #, 'jer'])

process = framework.create()


## Add PUPPI MET
#process.framework.producers.puppimet = cms.PSet(METProducer.default_configuration.clone())
#process.framework.producers.puppimet.prefix = cms.string('puppimet_')
#process.framework.producers.puppimet.parameters.met = cms.untracked.InputTag('slimmedMETsPuppi')

if runOnData : 
    process.source.fileNames = cms.untracked.vstring(
        '/store/data/Run2015D/MuonEG/MINIAOD/PromptReco-v4/000/258/159/00000/64914E6C-F26B-E511-B0C8-02163E0142D1.root'
        )
else : 
    process.source.fileNames = cms.untracked.vstring(
        'file:////storage/data/cms/store/user/brfranco/testFiles/TTTo2L2Nu_13TeV-powheg_RunIISpring15MiniAODv2_74X_mcRun2_asymptotic_v2-v1.root'
        )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
