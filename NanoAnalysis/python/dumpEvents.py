from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
#from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

### Dump reconstructed quantities, for debug purposes

charge={1:"+", -1:"-"}

class dumpEvents(Module):
    def __init__(self, dump=False):
        self.writeHistFile = False
        self.dump=False

    def analyze(self, event):
        if not self.dump: return True

        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")

        for i,lep in enumerate(muons):
#            if not muLooseId(lep) and not DEBUG : continue # Print only leptons passing loose ID
            print(eventId+" ", end="")
            print('#{} mu{} pt={:.3g} eta={:.3g} phi={:.3g} dxy={:.3g} dz={:.3g} SIP={:.3g} cmask={}'.format(i, charge[lep.charge],
                                                                                                             lep.pt,
                                                                                                             lep.eta,
                                                                                                             lep.phi,
                                                                                                             abs(lep.dxy),
                                                                                                             abs(lep.dz),
                                                                                                             lep.sip3d,
                                                                                                             lep.cleanmask),
                        end="")
            print(' GLB={} TK={} PF={}'.format(int(lep.isGlobal),
                                               int(lep.isTracker),
                                               int(lep.isPFcand),
                                               #int(lep.highPtId>0),# 1=trk,2=global ;  isHighPtId={}
                                               ),
                  end="")
            print(' combRelIsoPF={:.3g} combRelIsoPFFSRCorr={:.3g}'.format(lep.pfRelIso03_all, lep.pfRelIso03FsrCorr),
                  end="")
            print(' isLoose={} isTight={}'.format(int(muLooseId(lep)),int(muTightId(lep))),
                  end="")
            if (lep.fsrPhotonIdx>=0):
                fsr=fsrPhotons[lep.fsrPhotonIdx]
                print(' FSR: pt={:.3g} eta={:.3g}, phi={:.3g}, dREt2={:.3g}, Iso={:.3g}, true={}'.format(fsr.pt, fsr.eta, fsr.phi, fsr.dROverEt2, fsr.relIso03, (fsr.genFsrIdx>=0)))
            else:
                print()

        for i,lep in enumerate(electrons):
#            if not eleLooseId(lep) and not DEBUG : continue # Print only leptons passing loose ID
            print(eventId+" ", end="")
            print('#{} e{}  pt={:.3g} ptcorr={:.3g} eta={:.3g} phi={:.3g} dxy={:.3g} dz={:.3g} SIP={:.3g}'.format(i, charge[lep.charge],
                                                                                                    lep.pt,lep.eCorr,
                                                                                                    lep.eta,
                                                                                                    lep.phi,
                                                                                                    abs(lep.dxy),
                                                                                                    abs(lep.dz),
                                                                                                    lep.sip3d),
                  end="")
            print(' scEta={:.3g} BDT={:.3g} passBDT={}'.format(lep.eta+lep.deltaEtaSC, lep.mvaFall17V2Iso,
                                                               passEleBDT(lep.pt, lep.eta+lep.deltaEtaSC, lep.mvaFall17V2Iso)))
