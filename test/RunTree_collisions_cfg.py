import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("DTNT",eras.Run2_2017)

# Input
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",

  fileNames = cms.untracked.vstring
  (

    '/store/data/Run2017F/SingleMuon/RAW-RECO/ZMu-17Nov2017-v1/00000/06FB9E1B-BEF1-E711-9122-0025905A608E.root',

  ),
  secondaryFileNames = cms.untracked.vstring(
  )
)


##process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/Services_cff')

process.load('Configuration/StandardSequences/GeometryRecoDB_cff')
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.DTGeometryBuilder.dtGeometryDB_cfi")

process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v2'

process.load('EventFilter.ScalersRawToDigi.ScalersRawToDigi_cfi')
process.load('RecoLuminosity.LumiProducer.lumiProducer_cff')

# Unpackers, DT readout , TwinMux, BMTF, RPC readout  
process.load("EventFilter.DTRawToDigi.dtunpackerDDUGlobal_cfi")
process.dtunpacker.readOutParameters.debug = False
process.dtunpacker.readOutParameters.rosParameters.debug = False

process.load("EventFilter.L1TXRawToDigi.twinMuxStage2Digis_cfi")

process.load('EventFilter.L1TRawToDigi.bmtfDigis_cfi')

process.load("EventFilter.RPCRawToDigi.rpcUnpackingModule_cfi")

# RPC RecHit
process.load("RecoLocalMuon.RPCRecHit.rpcRecHits_cfi")
process.rpcRecHits.rpcDigiLabel = cms.InputTag('rpcUnpackingModule')

# Emulator
process.load('Configuration.StandardSequences.SimL1Emulator_cff')

process.twinMuxEmul = process.simTwinMuxDigis.clone()

process.twinMuxEmul.DTDigi_Source = "twinMuxStage2Digis:PhIn:"
process.twinMuxEmul.DTThetaDigi_Source = "twinMuxStage2Digis:ThIn:"
process.twinMuxEmul.RPC_Source = "rpcUnpackingModule"

process.DTMuonSelection = cms.EDFilter("DTMuonSelection",
                                 src = cms.InputTag('muons'),
                                 Muons = cms.InputTag('muons'),
                                 SAMuons = cms.InputTag('standAloneMuons'),
                                 dtSegmentLabel = cms.InputTag('dt4DSegments'),
                                 etaMin = cms.double(-1.25),
                                 etaMax = cms.double(1.25),
                                 ptMin = cms.double(0.),#3.),
                                 tightness = cms.int32(1) # 0 = loose (e.g. unstable collisions, minimum bias, requires a DT segment)
                                                          # 1 = medium (e.g. cosmics, requires a stand alone muon)
                                                          # 2 = tight (good collisions, requires a global muon)
)


# DT Ntuple
process.load("UserCode/DTDPGAnalysis/DTTTreGenerator_cfi")

process.myDTNtuple.localDTmuons = cms.untracked.bool(False)
process.myDTNtuple.outputFile = "DTNtuple.root"

process.myDTNtuple.dtDigiLabel = cms.InputTag("dtunpacker")

process.myDTNtuple.dtTrigTwinMuxOutEmulLabel = cms.InputTag("twinMuxEmul")

# The path
process.p = cms.Path(  process.DTMuonSelection 
                       + process.dtunpacker 
                       + process.rpcUnpackingModule 
                       + process.twinMuxStage2Digis 
                       + process.bmtfDigis 
                       + process.scalersRawToDigi 
                       + process.lumiProducer 
                       + process.twinMuxEmul
                       # + process.rpcRecHits 
                       + process.myDTNtuple
                    )


