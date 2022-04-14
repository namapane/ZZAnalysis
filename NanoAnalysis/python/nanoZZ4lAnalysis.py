####
# Steering file for ZZ analysis starting from nanoAODs
# Example for customization and running: test/runLocal.py
####

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import *

from ZZAnalysis.NanoAnalysis.tools import setConf, getConf
from ZZAnalysis.NanoAnalysis.triggerAndSkim import * # Trigger requirements are defined here
from ZZAnalysis.NanoAnalysis.lepFiller import *
from ZZAnalysis.NanoAnalysis.ZZFiller import *

#import os
#import sys

### Definition of analysis cuts
cuts = dict(
    ### lepton IDs
    muPt = 5.,
    elePt = 7.,
    relIso = 0.35,
    sip3d = 4.,
    dxy =  0.5,
    dz = 1.,
    fsr_dRET2 = 0.012,
    fsr_Iso = 1.8,
    # Relaxed ID (used for CRs for fake rate method) includes SIP, but a version without SIP will be required for relaxed-SIP CRs
    muRelaxedId = (lambda l : l.pt > cuts["muPt"] and abs(l.eta) < 2.4 and abs(l.dxy) < cuts["dxy"] and abs(l.dz) < cuts["dz"] and abs(l.sip3d) < cuts["sip3d"] and (l.isGlobal or (l.isTracker and l.nStations>0))), # Notes: l.nStations is numberOfMatchedStation, not numberOfMatches; also, muonBestTrackType!=2 is not available in nanoAODs
    muFullId = (lambda l, era : cuts["muRelaxedId"](l) and (l.isPFcand or (l.highPtId>0 and l.pt>200.))), # Note: full ID does not include isolation; FSR-corrected iso has to be applied on top.

    eleRelaxedId = (lambda l : l.pt > cuts["elePt"] and abs(l.eta) < 2.5 and abs(l.dxy) < cuts["dxy"] and abs(l.dz) < cuts["dz"] and abs(l.sip3d) < cuts["sip3d"]),
    eleFullId =  (lambda l, era : cuts["eleRelaxedId"](l) and abs(l.sip3d) < cuts["sip3d"] and passEleBDT(l, era)), #FIXME: BDT definition is available only for 2017!

    )


### Get processing customizations if defined by including .py; use defaults otherwise 
SAMPLENAME = getConf("SAMPLENAME", "test")
LEPTON_SETUP = getConf("LEPTON_SETUP", 2018)
if not (LEPTON_SETUP == 2016 or LEPTON_SETUP == 2017 or LEPTON_SETUP == 2018) :
    print("Invalid LEPTON_SETUP", LEPTON_SETUP)
    exit(1)
IsMC = getConf("IsMC", True)
PD = getConf("PD", "")
XSEC = getConf("XSEC", 1.)
SYNCMODE = getConf("SYNCMODE", False)
runMELA = getConf("runMELA", True)
bestCandByMELA = getConf("bestCandByMELA", True) # requires also runMELA=True

# ggH NNLOPS weight
APPLY_QCD_GGF_UNCERT = getConf("APPLY_QCD_GGF_UNCERT", False) 

# K factors for ggZZ (and old NLO ggH samples) 0:None; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
APPLY_K_NNLOQCD_ZZGG = getConf("APPLY_K_NNLOQCD_ZZGG", 0) 

# K factors for qqZZ
APPLY_K_NNLOQCD_ZZQQB = getConf("APPLY_K_NNLOQCD_ZZQQB", False) 
APPLY_K_NNLOEW_ZZQQB  = getConf("APPLY_K_NNLOEW_ZZQQB", False) 


# Preselection to speed up processing.
preselection = getConf("preselection", "nMuon + nElectron >= 4 &&" +
                       "Sum$(Muon_pt > {muPt}-2.) +" + # Allow for variations due to scale calib
                       "Sum$(Electron_pt > {elePt})" +
                       ">= 4").format(**cuts)

store = getConf("store","") #/eos/cms/ for files available on eos;
fileNames = getConf("fileNames", ["/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root",]) # sample ggH125 file
for i, file in enumerate(fileNames):
    fileNames[i] = store+file

localPath = os.environ['CMSSW_BASE']+"/src/ZZAnalysis/NanoAnalysis/"

jsonFile = None
if not IsMC :
    if LEPTON_SETUP == 2018 :
        jsonFile = localPath+"test/prod/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt"
    else:        
        exit(1) #2016-17 to be implemented

from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import muonScaleRes2016, muonScaleRes2017, muonScaleRes2018
muonScaleRes = {2016:muonScaleRes2016, 2017:muonScaleRes2017, 2018:muonScaleRes2018}


ZZSequence = [triggerAndSkim(isMC=IsMC, PD=PD, era=LEPTON_SETUP),
              muonScaleRes[LEPTON_SETUP](overwritePt=True, syncMode=SYNCMODE), # FIXME requires custom muonScaleResProducer.py
              lepFiller(cuts, LEPTON_SETUP), 
              ZZFiller(runMELA, bestCandByMELA, IsMC, LEPTON_SETUP),
              ]

if IsMC :
    from ZZAnalysis.NanoAnalysis.mcTruthAnalyzer import *
    ZZSequence.insert(0, mcTruthAnalyzer(dump=False))

    from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight_2016, puWeight_2017, puWeight_2018
    puWeight = {2016:puWeight_2016, 2017:puWeight_2017, 2018:puWeight_2018} # FIXME official weights are slightly different than the one we use (checked for 2018)
    ZZSequence.append(puWeight[LEPTON_SETUP]())

    from ZZAnalysis.NanoAnalysis.weightFiller import weightFiller
    ZZSequence.append(weightFiller(XSEC, APPLY_K_NNLOQCD_ZZGG, APPLY_K_NNLOQCD_ZZQQB, APPLY_K_NNLOEW_ZZQQB, APPLY_QCD_GGF_UNCERT))

branchsel_in = ""
if IsMC:
    branchsel_in  = localPath+"python/branchsel_in_MC.txt"
    branchsel_out = localPath+"python/branchsel_out_MC.txt"
else:
    branchsel_in  = localPath+"python/branchsel_in_Data.txt"
    branchsel_out = localPath+"python/branchsel_out_Data.txt"



from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
p = PostProcessor(".", fileNames,
                  prefetch=True, longTermCache=False,
                  cut=preselection, # pre-selection cuts (to speed up processing)
                  branchsel=branchsel_in, # select branches to be read
                  outputbranchsel=branchsel_out, # select branches to be written out
                  jsonInput=jsonFile, # path of json file for data
                  modules=ZZSequence,
                  noOut=False, # True = do not write out skimmed nanoAOD file
                  haddFileName="ZZ4lAnalysis.root", # name of output nanoAOD file
                  histFileName="histos.root", histDirName="plots", # file containing histograms
                  maxEntries=0, # Number of events to be read
                  firstEntry=0, # First event to be read
                  ) 

# p.run()
