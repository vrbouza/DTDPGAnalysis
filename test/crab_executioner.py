from WMCore.Configuration import Configuration
config = Configuration()


config.section_("General")
config.General.requestName      = 'VicTuples'
config.General.workArea         = './'
config.General.transferOutputs  = True
config.General.transferLogs     = True

config.section_("JobType")
config.JobType.pluginName       = 'Analysis'
config.JobType.psetName         = 'RunTree_collisions_cfg.py'

config.section_("Data")
config.Data.inputDataset        = '/RelValSingleMuPt100/CMSSW_9_3_7-PU25ns_93X_upgrade2023_realistic_v5_2023D17PU70-v1/GEN-SIM-RECO'
config.Data.inputDBS            = 'global'
config.Data.splitting           = 'FileBased' #LumiBased, Automatic y otros
config.Data.unitsPerJob         = 2
config.Data.outLFNDirBase       = '/store/user/rodrigvi/'
config.Data.publication         = False

config.section_("Site")
config.Site.storageSite         = 'T2_ES_IFCA'
