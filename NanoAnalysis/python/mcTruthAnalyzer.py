from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

### Analyze MC Truth
# Optionally dump MC history;
# Add ZZ gen valriables:
# -GenZZFinalState: product of IDs of the four gen ZZ leptons
# -FsrPhoton_genFsrIdx: index of the closest gen FSR from Z->l (e, mu)
# -GenHLeps (TBD)
#

# Find particle's real mother, skipping rows with the same pdgId
def Mother (part, gen) :
    idxMother= part.genPartIdxMother
    while idxMother>=0 and gen[idxMother].pdgId == part.pdgId:
        idxMother = gen[idxMother].genPartIdxMother
    idMother=0
    if idxMother >=0 : idMother = gen[idxMother].pdgId
    return idxMother, idMother


class mcTruthAnalyzer(Module):
    def __init__(self, dump=False):
        self.writeHistFile = False
        self.printGenHist  = dump # print MC history


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("GenZZFinalState", "I")
        self.out.branch("FsrPhoton_genFsrIdx", "I", lenVar="nFsrPhoton")


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        genpart=Collection(event,"GenPart")

        ### Print Gen history and LHE particles.
        ### See also: https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/common/hepmcDump.py
        if self.printGenHist : 
            print ("---Gen:")
            for i, gp in enumerate(genpart) :
                motherId=-1
                gmotherId=-1
                if gp.genPartIdxMother >= 0 : 
                    motherId = genpart[gp.genPartIdxMother].pdgId
                    if genpart[gp.genPartIdxMother].genPartIdxMother >= 0 :
                        gmotherId = genpart[genpart[gp.genPartIdxMother].genPartIdxMother].pdgId
                print (i, gp.pdgId, gp.genPartIdxMother, gp.pt, gp.eta, gp.status)
        
            print("---------LHEPart---------")
            LHEPart = Collection (event, "LHEPart")
            for i, Lp in enumerate(LHEPart):
                print(i, Lp.pdgId, Lp.pt, Lp.eta, Lp.status, Lp.incomingpz)


        ## Search for gen FSR from Z->ll (e, mu)
        genFSRIdxs = []
        for i, gp in enumerate(genpart) :
            midx = gp.genPartIdxMother
            if gp.pdgId==22 and gp.pt > 2. and midx >= 0 :
                mid = genpart[midx].pdgId
                if abs(mid) == 11 or abs(mid) == 13 :
                    mmidx, mmid = Mother(genpart[midx], genpart) # possibly skip intermediate rows
                    if mmid == 23 :
                        genFSRIdxs.append(i)

        # Store FsrPhoton_genFsrIdx
        fsrPhotons = Collection(event, "FsrPhoton")
        fsrPhotons_genFsrIdx = [-1]*len(fsrPhotons)
        for ifsr, fsr in enumerate(fsrPhotons) :
            dRmin = 0.3
            for igen in genFSRIdxs :
                genfsr = genpart[igen]
                dR = deltaR (fsr.eta, fsr.phi, genfsr.eta, genfsr.phi)
                if dR < dRmin :
                    dRmin = dR
                    fsrPhotons_genFsrIdx[ifsr] = igen

        self.out.fillBranch("FsrPhoton_genFsrIdx", fsrPhotons_genFsrIdx)


        # Search for leptons from H
        GenHLeps = filter(lambda f : 
                          (abs(f.pdgId)==11 or abs(f.pdgId)==13 or abs(f.pdgId)==15) and
                          f.genPartIdxMother >=0 and genpart[f.genPartIdxMother].pdgId == 23 and
                          genpart[f.genPartIdxMother].genPartIdxMother >= 0 and
                          genpart[ genpart[f.genPartIdxMother].genPartIdxMother].pdgId == 25, genpart) #leptons generated from H->ZZ

#        GenEl_acc = filter(lambda f: abs(f.pdgId)== 11 and abs(f.eta)<2.5 and f.pt > conf["elePt"], GenHLeps)
#        GenMu_acc = filter(lambda f: abs(f.pdgId)== 13 and abs(f.eta)<2.4 and f.pt > conf["muPt"], GenHLeps)
#        gen_p4 = ROOT.TLorentzVector()
        GenZZFinalState=1
        for lep in GenHLeps :
#            gen_p4 += lep.p4()
            GenZZFinalState *= lep.pdgId
#        print("Gen: {:} leps M={:.4g} In Acc: {:} e, {:} mu,  {:} FSR".format(len(GenHLeps),gen_p4.M(), len(GenEl_acc), len(GenMu_acc), len (GenFSR)), "\n")

        self.out.fillBranch("GenZZFinalState", GenZZFinalState)

        return True
