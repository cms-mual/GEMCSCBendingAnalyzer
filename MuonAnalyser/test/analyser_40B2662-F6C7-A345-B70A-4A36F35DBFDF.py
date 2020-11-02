import FWCore.ParameterSet.Config as cms 
from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('analyser',Run3)
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load("RecoTracker.TrackProducer.TrackRefitters_cff") 
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 5000
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('nEvents',
                        -1, #Max number of events
                        VarParsing.multiplicity.singleton,
                        VarParsing.varType.int,
                        "Number of events")
options.parseArguments()
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(options.nEvents)
)
process.maxEvents.input = cms.untracked.int32(5000)
process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(options.inputFiles),
                                inputCommands = cms.untracked.vstring(
                        "keep *",
                        "drop TotemTimingDigiedmDetSetVector_totemTimingRawToDigi_TotemTiming_reRECO",
                        "drop TotemTimingRecHitedmDetSetVector_totemTimingRecHits__reRECO"
                        )
                                )
process.source.fileNames.append('file:/eos/cms/store/data/Commissioning2020/Cosmics/RAW-RECO/CosmicSP-PromptReco-v1/000/337/962/00000/640B2662-F6C7-A345-B70A-4A36F35DBFDF.root')
#process.source.fileNames.append('file:/eos/cms/store/data/Run2018D/SingleMuon/RAW-RECO/ZMu-12Nov2019_UL2018-v4/00000/8AAA8B57-D308-1543-B660-BDC6B10E25C7.root')
#process.source.fileNames.append('root://xrootd-cms.infn.it//store/data/Run2018D/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/322/430/00000/B44676BC-D8B4-E811-9F1D-FA163E987C59.root')
process.options = cms.untracked.PSet()
process.TFileService = cms.Service("TFileService", fileName = cms.string("out_ana_40B2662-F6C7-A345-B70A-4A36F35DBFDF.root"))
process.analyser = cms.EDAnalyzer('analyser',
        process.MuonServiceProxy,
        gemRecHits = cms.InputTag("gemRecHits"),
        #muons = cms.InputTag("muons"),
        muons = cms.InputTag("muons1Leg"),
        vertexCollection = cms.InputTag("offlinePrimaryVerticies")
)
process.p = cms.EndPath(process.TrackRefitter*process.analyser)
