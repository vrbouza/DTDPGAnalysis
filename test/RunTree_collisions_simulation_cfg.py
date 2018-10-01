import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("DTNT",eras.Run2_2018)

process.load('Configuration/StandardSequences/Services_cff')

process.load('Configuration/StandardSequences/GeometryRecoDB_cff')  ##  solve STA problem
process.load('Configuration/EventContent/EventContent_cff')
process.load("Geometry.DTGeometryBuilder.dtGeometryDB_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

#for RAW
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")  #DB v2, at least since GR_E_V42

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic')

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",

  fileNames = cms.untracked.vstring
  (
    'file://./10808.0_SingleMuPt100+SingleMuPt100_pythia8_2018_GenSimFull+DigiFull_2018+RecoFull_2018+HARVESTFull_2018+ALCAFull_2018+NanoFull_2018/step3.root',
  ),
  secondaryFileNames = cms.untracked.vstring
  (
  )
)

process.load("UserCode/DTDPGAnalysis/DTTTreGenerator_cfi")

process.myDTNtuple.localDTmuons = cms.untracked.bool(False)
process.myDTNtuple.outputFile = "DTNtuple.root"
process.myDTNtuple.bmtfOutputDigis = cms.InputTag("bmtfDigis:BMTF")
process.myDTNtuple.rpcLabel = cms.InputTag("muonRPCDigis")
process.myDTNtuple.runOnSimulation = cms.bool(True)


process.p = cms.Path(process.myDTNtuple)
