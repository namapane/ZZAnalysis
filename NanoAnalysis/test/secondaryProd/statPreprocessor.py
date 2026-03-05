# Preprocessor for the statistical analysis.
# For use in batch processing, use:
# batch_Condor.py samples_secondary.csv -i statPreprocessor.py


from ZZAnalysis.NanoAnalysis.tools import setConf, getConf

NANOVERSION = getConf("NANOVERSION", 12)
LEPTON_SETUP = getConf("LEPTON_SETUP", 2022)
fileNames = getConf("fileNames", ["/eos/home-a/atarabin/STXS_samples/PROD_samplesNano_2022_MC_8d4c03f7/ggH125/ZZ4lAnalysis.root",])


from ZZAnalysis.NanoAnalysis.ZZCategorization import *
from ZZAnalysis.NanoAnalysis.flattenCand import *
sequence = [ZZCategorization(),
            flattenCand(),
]

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
p = PostProcessor(".", fileNames,
                  prefetch=False, longTermCache=False,
                  cut=None,
                  branchsel=None,
                  outputbranchsel= ['drop *',
                                    'keep BestCand_*',
                                    ],
                  jsonInput=None,
                  modules=sequence,
                  noOut=False,
                  #haddFileName="ZZ4lAnalysis_ext.root",
                  maxEntries=0,
                  firstEntry=0,
                  friend=False,
#                  postfix="_ext",
                  provenance = False
                  )

p.run()

