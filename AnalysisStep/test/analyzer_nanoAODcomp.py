
#DATA_TAG = "ReReco" # Change to PromptReco for Run2016 period H
LEPTON_SETUP = 2018  # current default = 2018 
#ELECORRTYPE = "None" # "None" to switch off
#ELEREGRESSION = "None" # "None" to switch off
APPLYMUCORR = False  # Switch off muon scale corrections (To compare with nanoAOD)
#APPLYJEC = False     #
#APPLYJER = False     #
#RECORRECTMET = False #
KINREFIT = False    # control KinZFitter (very slow)
PROCESS_CR = False   # Uncomment to run CR paths and trees
#ADDLOOSEELE = True  # Run paths for loose electrons
#APPLYTRIG = False    # Skip events failing required triggers. They are stored with sel<0 if set to False
#KEEPLOOSECOMB = True # Do not skip loose lepton ZZ combinations (for debugging)
ADDZTREE = False # Add tree for Z analysis
ADDLHEKINEMATICS = True  #

PD = ""
MCFILTER = ""

#For DATA: 
#IsMC = False
#PD = "DoubleMu"
#DATA_TAG = "ReReco" # Change to "PromptReco" for Run2018 period D

# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "analyzer.py")
execfile(PyFilePath + "prod/pyFragments/RecoProbabilities.py")

if not IsMC:
	process.source.inputCommands = cms.untracked.vstring("keep *", "drop LHERunInfoProduct_*_*_*", "drop LHEEventProduct_*_*_*") ###FIXME In 9X this removes all collections for MC

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

process.source.fileNames = cms.untracked.vstring(

#This file includes 18000 events from /GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM
# '/store/mc/RunIIAutumn18MiniAOD/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/270000/476D26E9-498D-1849-8DB0-30E26269609D.root'
# and matches the nanoAOD file: '/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/9B09F47A-6958-9E41-90BB-155D34BC6F01.root'

#This file includes 6000 events from /GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM
 '/store/mc/RunIIAutumn18MiniAOD/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/270000/0FDDC278-35AA-5840-A23E-1B5082080832.root'
# and matches the nanoAOD file: '/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root'
        
###Run2018A-ReReco-v2                                                                                                                 
        #'/store/data/Run2018A/DoubleMuon/MINIAOD/17Sep2018-v2/120000/39DE1F78-583A-1948-8E09-E47E33DCCBED.root'                                      

###Run2018D-PromptReco-v2                                                                                                 
        #'/store/data/Run2018D/DoubleMuon/MINIAOD/PromptReco-v2/000/325/172/00000/DB7F6A17-F4E1-B844-B8A4-7F53C1E17E8C.root'
    )

#process.calibratedPatElectrons.isSynchronization = cms.bool(True) #not needed anymore since new EGamma smearing is event deterministic
#process.calibratedMuons.isSynchronization = cms.bool(True)

process.maxEvents.input = -1
#process.source.skipEvents = cms.untracked.uint32(5750)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


#process.source.eventsToProcess = cms.untracked.VEventRange("1:12:66503")
process.maxEvents.input = 100

# Debug
process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
     muonSrcs = cms.PSet(
#       slimmedMuons = cms.InputTag("slimmedMuons"),
        muons = cms.InputTag("appendPhotons:muons"),
     ),
     electronSrcs = cms.PSet(
##       slimmedElectron = cms.InputTag("slimmedElectrons"),
#        electrons = cms.InputTag("appendPhotons:electrons"),
     ),
     candidateSrcs = cms.PSet(
        Z     = cms.InputTag("ZCand"),
#        ZZ  = cms.InputTag("ZZCand"),
#        ZLL  = cms.InputTag("ZLLCand"),
#        ZL  = cms.InputTag("ZlCand"),
     ),
#     jetSrc = cms.InputTag("cleanJets"),
)

# Create lepton sync file
#process.PlotsZZ.dumpForSync = True;
#process.p = cms.EndPath( process.PlotsZZ)

# Keep all events in the tree, even if no candidate is selected
#process.ZZTree.skipEmptyEvents = False

# replace the paths in analyzer.py
#process.trees = cms.EndPath(process.ZZTree)

#Dump reconstructed variables
#process.appendPhotons.debug = cms.untracked.bool(True)
#process.fsrPhotons.debug = cms.untracked.bool(True)
process.dump = cms.Path(process.dumpUserData)

#Print MC history
#process.mch = cms.EndPath(process.printTree)


#Monitor memory usage
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
