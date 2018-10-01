import FWCore.ParameterSet.Config as cms
import sys

def customiseEventContentForDTNtuple(process):

    if hasattr(process,'RECOSIMoutput') :
        
        sys.stderr.write("[customiseEventContentForDTNtuple] : Found RECOSIMoutput, customising event content \n")

        process.RECOSIMoutput.outputCommands.append('keep *_muonDTDigis_*_*')
        process.RECOSIMoutput.outputCommands.append('keep *_twinMuxStage2Digis_*_*')
        process.RECOSIMoutput.outputCommands.append('keep *_bmtfDigis_*_*')

        process.RECOSIMoutput.outputCommands.append('keep *_rpcRecHits_*_*')
        process.RECOSIMoutput.outputCommands.append('keep *_muonRPCDigis_*_*')

    return process
