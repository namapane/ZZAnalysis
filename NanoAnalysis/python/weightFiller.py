from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import *


class weightFiller(Module):
    def __init__(self, XS):
        self.writeHistFile = False
        self.XS = XS

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("totalWeight", "F")

    def analyze(self, event):        
        w_total = self.XS * event.Generator_weight * event.puWeight * event.ZZ_dataMCWeight
        self.out.fillBranch("totalWeight", w_total)

        return True


