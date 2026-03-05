### 
# Add categories to ZZCand and ZLLCand collections
###

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

class ZZCategorization(Module):
    def __init__(self, processCR=False):
        self.processCR = processCR
        print("***ZZCategorization", flush=True)


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.book("ZZCand")
        if self.processCR :
            self.book("ZLLCand")

    def analyze(self, event) :
        self.fill('ZZCand', event)
        if self.processCR :
            self.fill('ZLLCand', event)

        return True

    def book(self, collName) :
        theLenVar="n"+collName
        self.out.branch(collName+"_category", "S", lenVar=theLenVar, title="Categorization")

    def fill(self, collName, event) :
        cands = Collection(event, collName)
        cats = [-1]*len(cands)
        for iCand, aCand in enumerate(cands):
            cats[iCand] = 1 # FIXME: actual function to be called here, passing the candidate and all relevant variables
        self.out.fillBranch(collName+"_category", cats)
