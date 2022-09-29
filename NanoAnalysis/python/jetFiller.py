### TO BE IMPLEMENTED:
### -Jet-Lepton cross-cleaning
### -Additional JES, JEC

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

class jetFiller(Module):
    def __init__(self, era):

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Jet_lepVeto", "I", lenVar="nJet")
# ...

    def analyze(self, event):
        jets = Collection(event, 'Jet')
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        leps = list(muons)+ list(electrons)
        nlep=len(leps)
# ...        

        return True
