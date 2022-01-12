from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
#from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR

### Dump reconstructed quantities, for debug purposes

charge={1:"+", -1:"-"}

class dumpEvents(Module):
    def __init__(self, level=1):
        self.writeHistFile = False
        self.level=level
        
    def analyze(self, event):

        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")

        eventId='{}:{}:{}'.format(event.run,event.luminosityBlock,event.event)

        if self.level >= 3 : #full detail
            for i,lep in enumerate(muons):
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
               print(' isLoose={} isTight={}'.format(lep.isLoose,lep.isTightIso),
                     end="")
               if (lep.myFsrPhotonIdx>=0):
                   fsr=fsrPhotons[lep.myFsrPhotonIdx]
                   print(' FSR: pt={:.3g} eta={:.3g}, phi={:.3g}, dREt2={:.3g}, Iso={:.3g}, true={}'.format(fsr.pt, fsr.eta, fsr.phi, fsr.dROverEt2, fsr.relIso03, (fsr.genFsrIdx>=0)))
               else:
                   print()
   
            for i,lep in enumerate(electrons):
#  #            if not eleLooseId(lep) and not DEBUG : continue # Print only leptons passing loose ID
               print(eventId+" ", end="")
               print('#{} e{}  pt={:.3g} ptcorr={:.3g} eta={:.3g} phi={:.3g} dxy={:.3g} dz={:.3g} SIP={:.3g}'.format(i, charge[lep.charge],
                                                                                                       lep.pt,lep.eCorr,
                                                                                                       lep.eta,
                                                                                                       lep.phi,
                                                                                                       abs(lep.dxy),
                                                                                                       abs(lep.dz),
                                                                                                       lep.sip3d),
                     end="")
               # FIXME must refactor this passBDT = passEleBDT(lep.pt, lep.eta+lep.deltaEtaSC, lep.mvaFall17V2Iso) 
               passBDT=True
               print(' scEta={:.3g} BDT={:.3g} passBDT={}'.format(lep.eta+lep.deltaEtaSC, lep.mvaFall17V2Iso, passBDT))
   

        if self.level >=1 : # final candidates
            # Sync-style printout
            leps = list(muons)+ list(electrons)
            ZZFlav = (leps[event.Z1_l1Idx].pdgId*leps[event.Z2_l1Idx].pdgId)**2
            print ('{}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{}:{:.6g}:{:.3f}:{:.3f}:{:.3f}:{:.3f}'.format(
                eventId,
                event.ZZ_mass,
                event.Z1_mass,
                event.Z2_mass,
                event.ZZ_massPreFSR,
                ZZFlav,
                event.ZZ_Dbkgkin,
                event.Generator_weight,
                event.puWeight,
                event.ZZ_dataMCWeight,
                event.totalWeight,
                ))

        return True

