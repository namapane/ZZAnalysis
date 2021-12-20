#!/usr/bin/env python
###
# Example for running the analysis locally, after customizing variables.
###
from __future__ import print_function
from ZZAnalysis.NanoAnalysis.tools import setConf, getConf


### Customize processing variables

### ggH125
setConf("SAMPLENAME", "ggH125")
setConf("XSEC", 48.58*0.0002745)
setConf("store","/eos/cms")
setConf("fileNames",[
         "/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root",
#         "/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/445F4788-F322-C443-AB54-699C7716976A.root",
#         "/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/1DA1D554-E919-0146-9873-77E0B1B08DB5.root",
         ])

### ZH125
# setConf("SAMPLENAME", "ZH125")
# setConf("store","root://cms-xrd-global.cern.ch/")
# setConf("fileNames",
#         ["/store/mc/RunIIAutumn18NanoAODv7/ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/C670B342-C02E-1848-AE6A-0B5550E3DFE3.root"])

### qqZZ: /ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/NANOAODSIM
# setConf("SAMPLENAME", "qqZZ")
# setConf("XSEC", 1.256)
# setConf("fileNames",
#         ["/store/mc/RunIIAutumn18NanoAODv7/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/100000/3BD1674D-804E-0E45-B106-EAB3C0D59E32.root"])

### High mass sample
# setConf("SAMPLENAME", "VBFH3000")
# setConf("fileNames",
#         ["/store/mc/RunIIAutumn18NanoAODv7/VBF_HToZZTo4L_M3000_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/A82C31DF-0A6F-B44F-857B-89BFD72CEFEA.root"])



# Select specific events to debug
#setConf("preselection","event.eventId == 840922")


### this should be done after all customizations
from ZZAnalysis.NanoAnalysis.nanoZZ4lAnalysis import *
# Edit postprocessor parameters as necessary
p.haddFileName=getConf("SAMPLENAME", "ZZAnalysis.root")+".root"
#p.maxEntries = 1000

p.run()
