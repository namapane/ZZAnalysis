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
from ZZAnalysis.NanoAnalysis.dumpEvents import dumpEvents

#import os
#import sys

### Definition of analysis cuts
cuts = dict(
    muPt = 5.,
    elePt = 7.,
    relIso = 0.35,
    sip3d = 4.,
    dxy =  0.5,
    dz = 1.,
    fsr_dRET2 = 0.012,
    fsr_Iso = 1.8,
    # Loose IDs includes SIP, but a version without SIP will be required for relaxed-SIP CRs
    muLooseId = (lambda l : l.pt > cuts["muPt"] and abs(l.eta) < 2.4 and abs(l.dxy) < cuts["dxy"] and abs(l.dz) < cuts["dz"] and abs(l.sip3d) < cuts["sip3d"] and (l.isGlobal or l.isTracker)),  #FIXME: l.nStations>0) is not equivalent to numberOfMatches>0); also, muonBestTrackType!=2 is not available
    muTightId = (lambda l : cuts["muLooseId"](l) and (l.isPFcand or (l.highPtId>0 and l.pt>200.))),  # Note: FSR-corrected isolation has also to be applied in addition to tight ID  #FIXME: highPtId does not match exactly our definition.     


    eleLooseId = (lambda l : l.pt > cuts["elePt"] and abs(l.eta) < 2.5 and abs(l.dxy) < cuts["dxy"] and abs(l.dz) < cuts["dz"] and abs(l.sip3d) < cuts["sip3d"]),
    eleTightId =  (lambda l : cuts["eleLooseId"](l) and abs(l.sip3d) < cuts["sip3d"] and passEleBDT(l.pt, l.eta+l.deltaEtaSC, l.mvaFall17V2Iso)), #FIXME: We use different BDT definitions for 2016 and 2018

    test = (lambda x : cuts["elePt"]),
    )


### Get processing customizations if defined by includer macro; use defaults otherwise 
SAMPLENAME = getConf("SAMPLENAME", "test")
LEPTON_SETUP = getConf("LEPTON_SETUP", 2018)
IsMC = getConf("IsMC", True)
PD = getConf("PD", "")
XSEC = getConf("XSEC", 1.)
runMELA = getConf("runMELA", False)
bestCandByMELA = getConf("bestCandByMELA", False) # requires also runMELA=True

# Preselection to speed up processing.
# Note that this is done before muon momentum scale calibration
preselection = getConf("preselection", "nMuon + nElectron >= 3 &&" +
                       "Sum$(Muon_pt > {muPt}) +" +
                       "Sum$(Electron_pt > {elePt})" +
                       ">= 3").format(**cuts)

store = getConf("store","") #/eos/cms/ for files available on eos;
fileNames = getConf("fileNames", ["/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root",]) # sample ggH125 file
for i, file in enumerate(fileNames):
    fileNames[i] = store+file
jsonFile = getConf("jsonFile", None)




ZZSequence = [triggerAndSkim(isMC=IsMC, PD=PD, era=LEPTON_SETUP),
              lepFiller(cuts),
              ZZProducer(runMELA, bestCandByMELA, cuts, IsMC, XSEC),
#              dumpEvents(dump=False)
              ]


if IsMC :
    from ZZAnalysis.NanoAnalysis.mcTruthAnalyzer import *
    ZZSequence.insert(0, mcTruthAnalyzer(dump=False))
#    ZZSequence.append(weightProducer())


branchsel_in = ""
if IsMC:
    branchsel_in = os.environ['CMSSW_BASE']+"/src/ZZAnalysis/NanoAnalysis/python/branchsel_in_MC.txt"
else:
    branchsel_in = os.environ['CMSSW_BASE']+"/src/ZZAnalysis/NanoAnalysis/python/branchsel_in_Data.txt"

branchsel_out = os.environ['CMSSW_BASE']+"/src/ZZAnalysis/NanoAnalysis/python/branchsel_out.txt"


from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
p = PostProcessor(".", fileNames,
                  prefetch=False, longTermCache=False, #FIXME prefetch should be used for remote samples
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
