### 
# Example of copying selected info for the best candidate to produce flat compact trees for the stat analysis
###

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

class flattenCand(Module):
    def __init__(self, processCR=False):
        self.processCR = processCR 
        print("***flattenCand", flush=True)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("BestCand_FS", "S", title="Selected cand final state")
        self.out.branch("BestCand_mass", "F", title="Selected cand mass")
        self.out.branch("BestCand_weight", "F", title="Selected cand category")
        self.out.branch("BestCand_category", "S", title="Selected cand category")

    def analyze(self, event) :
        cands = Collection(event, "ZZCand")
        theCand = cands[event.bestCandIdx]
        FS = 0
        if theCand.Z1flav == -121 :
            if abs(theCand.Z2flav) == 121 : # also cover ZLL SS CRs
                FS = 1 #4e
            else :
                FS = 2 #2e2mu
        else :
            if abs(theCand.Z2flav) == 121 :
                FS = 3 #2mu2e
            else :
                FS = 4 #4mu
            
            
        self.out.fillBranch("BestCand_FS", FS)
        self.out.fillBranch("BestCand_mass", theCand.mass)
        self.out.fillBranch("BestCand_weight", event.overallEventWeight*theCand.dataMCWeight) # includes sigma*BR*k_factors*genweight*PU_W*DataMC_W; must be normalized by sum of gen weights
        self.out.fillBranch("BestCand_category", theCand.category)

        return True
