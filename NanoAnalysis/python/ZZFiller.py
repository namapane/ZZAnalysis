from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import *

from functools import cmp_to_key

#MELA
from importlib import import_module #FIXME?
from ctypes import *
from JHUGenMELA.MELA.mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar


ZmassValue = 91.1876;
DEBUG = False
DUMP = True

### Comparators to select the best candidate in the event 
def bestCandByZ1Z2(a,b): # Choose by abs(MZ1-MZ), or sum(PT) if same Z1
    if abs(a.Z1.M-b.Z1.M) < 1e-4 : # same Z1: choose the candidate with highest-pT Z2 leptons
        if a.Z2.sumpt() > b.Z2.sumpt() : return -1        
        else :
            return 1
    else :
        if abs(a.Z1.M-ZmassValue) < abs(b.Z1.M-ZmassValue): return -1 # else choose based on Z1 masses
        else :
            return 1

def bestCandByDbkgKin(a,b): # Choose by DbkgKin
    if abs((a.p4).M() - (b.p4).M())<1e-4 and a.Z1.finalState()*a.Z2.finalState()==b.Z1.finalState()*b.Z2.finalState(): # Equivalent: same masss (tolerance 100 keV) and same FS -> same leptons and FSR
        return bestCandByZ1Z2(a,b)
    if a.KD > b.KD : return 1 # choose by best dbkgkin
    else: return -1    


# Class to cache computed values for Z candidates.
# We may need to move Zs into the Event with a separate ZFiller module at some point, for CR studies
class ZCand: 
    def __init__(self, l1Idx, l2Idx, leps, fsrPhotons):
        # FIXME: we may want to set a default order for i, j (eg 1=-, 2=+)
        self.l1Idx = l1Idx
        self.l2Idx = l2Idx
        self.l1 = leps[l1Idx]
        self.l2 = leps[l2Idx]
        self.fsr1Idx = self.l1.myFsrPhotonIdx
        self.fsr2Idx = self.l2.myFsrPhotonIdx

        self.l1DressedP4 = self.l1.p4()
        self.l2DressedP4 = self.l2.p4()
        if self.fsr1Idx>=0 : self.l1DressedP4 += fsrPhotons[self.fsr1Idx].p4()
        if self.fsr2Idx>=0 : self.l2DressedP4 += fsrPhotons[self.fsr2Idx].p4()

        self.p4 = self.l1DressedP4 + self.l2DressedP4

        self.M = self.p4.M() # cache the mass as it is used often

    def sumpt() : # sum of lepton pTs, used to sort candidates
        return self.l1.pt+self.l2.pt

    def finalState(self) :
        return self.l1.pdgId*self.l2.pdgId


# Struct to cache computed values for ZZ candidates.
class ZZCand:
    def __init__(self, Z1, Z2, p_GG_SIG_ghg2_1_ghz1_1_JHUGen=0., p_QQB_BKG_MCFM=1.):
        self.Z1 = Z1
        self.Z2 = Z2
        self.p4 = Z1.p4+Z2.p4
        self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen = p_GG_SIG_ghg2_1_ghz1_1_JHUGen
        self.p_QQB_BKG_MCFM = p_QQB_BKG_MCFM
        self.KD = p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(p_GG_SIG_ghg2_1_ghz1_1_JHUGen+p_QQB_BKG_MCFM) # without c-constants, for candidate sorting

    def finalState(self) :
        return self.Z1.finalState()*self.Z2.finalState()


class ZZFiller(Module):
    def __init__(self, runMELA, bestCandByMELA, isMC, XS):
        self.writeHistFile = True
        self.runMELA = runMELA or bestCandByMELA
        self.bestCandByMELA = bestCandByMELA        
        self.isMC = isMC
        self.XS = XS

        if self.runMELA :
            self.mela = Mela(13, 125, TVar.ERROR)
            self.mela.setCandidateDecayMode(TVar.CandidateDecay_ZZ)      
            self.lib = CDLL('libZZAnalysisAnalysisStep.so')
            self.lib.D_bkg_kin.restype = c_float

    def beginJob(self,histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        # Histograms
        self.h_ZZMass = ROOT.TH1F('ZZMass','ZZMass',130,70,200)
        self.addObject(self.h_ZZMass)


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("ZZ_mass", "F")
        self.out.branch("ZZ_massPreFSR", "F")
        self.out.branch("Z1_mass", "F")
        self.out.branch("Z2_mass", "F")
        if self.isMC :
            # gen weight = Generator_weight
            self.out.branch("W_pu", "F")
            self.out.branch("W_dataMC", "F")
        self.out.branch("W_total", "F")


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        eventId='{}:{}:{}'.format(event.run,event.luminosityBlock,event.event)
        if DEBUG : print ('Event '+eventId)

        # Collections
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        leps = list(muons)+ list(electrons)
        nlep=len(leps)
        

        ### Z combinatorial over selected leps (after FSR-corrected ISO cut for muons)
        Zs = []
        for i,l1 in enumerate(leps):
            if l1.isTightIso :
                for j in range(i+1,nlep):
                    l2 = leps[j]
                    if l2.isTightIso and l1.pdgId == -l2.pdgId : #OS,SF
                        aZ = ZCand(i, j, leps, fsrPhotons)
                        zmass = aZ.M
                        if DEBUG: print('Z={:.4g} pt1={:.3g} pt2={:.3g}'.format(zmass, l1.pt, l2.pt), l1.myFsrPhotonIdx,  l2.myFsrPhotonIdx)
                        if (zmass>12. and zmass<120.):
                            Zs.append(aZ)

# At some point, we may want to store Z candidates in the event, in a dedicated ZFiller module; like:
#                         Z_d1[iZ]=aZ.l1Idx
#                         Z_d2[iZ]=aZ.l2Idx
#                         Z_fsr1Idx[iZ]=aZ.fsr1Idx
#                         Z_fsr2Idx[iZ]=aZ.fsr2Idx
#                         Z_pt[iZ]=aZ.p4.pt()
#                         Z_eta[iZ]=aZ.p4.eta()
#                         Z_phi[iZ]=aZ.p4.phi()
#                         Z_mass[iZ]=aZ.p4.M()
#                         Z_finalState[iZ]=aZ.finalState()


        ### Build all ZZ combinations passing the ZZ selection
        ZZs = []
        if len(Zs) >= 2:
            for iZ,aZ in enumerate(Zs):
                for jZ in range(iZ+1, len(Zs)):
                    # check that these Zs are mutually exclusive (not sharing the same lepton) 
                    if Zs[iZ].l1Idx==Zs[jZ].l1Idx or Zs[iZ].l2Idx==Zs[jZ].l2Idx or Zs[iZ].l2Idx==Zs[jZ].l1Idx or Zs[iZ].l2Idx==Zs[jZ].l2Idx: continue
                    #FIXME:  add DeltaR>0.02 cut among all leptons (still necessary)?
                    
                    # set Z1 and Z2
                    Z1Idx, Z2Idx = jZ, iZ
                    if abs(Zs[iZ].M-ZmassValue) < abs(Zs[jZ].M-ZmassValue):
                        Z1Idx, Z2Idx = Z2Idx, Z1Idx

                    Z1=Zs[Z1Idx]
                    Z2=Zs[Z2Idx]

                    # Z1 mass cut
                    if Z1.M <= 40. : continue

                    zzleps = [Z1.l1, Z1.l2, Z2.l1, Z2.l2]
                    lepPts = []
                    # QCD suppression on all OS pairs, regardelss of flavour
                    passQCD = True
                    for k in range(0,4):
                        lepPts.append(zzleps[k].pt)
                        for l in range (k+1,4):
                            if zzleps[k].pdgId*leps[l].pdgId < 0 and (leps[k].p4()+leps[l].p4()).M()<=4.:
                                passQCD = False
                                break
                    if not passQCD: continue

                    # trigger acceptance cuts (20,10 GeV)
                    lepPts.sort()
                    if not (lepPts[3]>20. and lepPts[2]>10.) : continue

                    #"Smart cut" on alternate pairings for same-sign candidates
                    if abs(Z1.l1.pdgId) == abs(Z2.l1.pdgId):
                        mZa, mZb = 0., 0.
                        if Z1.l1.pdgId == -Z2.l1.pdgId:
                            mZa=(Z1.l1DressedP4+Z2.l1DressedP4).M()
                            mZb=(Z1.l2DressedP4+Z2.l2DressedP4).M()
                        elif Z1.l1.pdgId == -Z2.l2.pdgId:
                            mZa=(Z1.l1DressedP4+Z2.l2DressedP4).M()
                            mZb=(Z1.l2DressedP4+Z2.l1DressedP4).M()
                        if (abs(mZa-ZmassValue)>abs(mZb-ZmassValue)) : mZa, mZb = mZb, mZa
                        if (abs(mZa-ZmassValue)<abs(Z1.M-ZmassValue)) and mZb < 12.: continue

                    #Compute D_bkg^kin
                    p_GG_SIG_ghg2_1_ghz1_1_JHUGen = 0.
                    p_QQB_BKG_MCFM = 1.
                    if self.runMELA:
                        daughters = SimpleParticleCollection_t()
                        daughters.push_back(SimpleParticle_t(Z1.l1.pdgId, Z1.l1DressedP4))
                        daughters.push_back(SimpleParticle_t(Z1.l2.pdgId, Z1.l2DressedP4))
                        daughters.push_back(SimpleParticle_t(Z2.l1.pdgId, Z2.l1DressedP4))
                        daughters.push_back(SimpleParticle_t(Z2.l2.pdgId, Z2.l2DressedP4))
                        #print("MELA ", daughters[ilep].first, daughters[ilep].second.Pt(), daughters[ilep].second.Eta(), daughters[ilep].second.Phi(),daughters[ilep].second.M())
                        self.mela.setInputEvent(daughters, 0, 0, 0)
                        self.mela.setProcess(TVar.HSMHiggs, TVar.JHUGen, TVar.ZZGG)
                        p_GG_SIG_ghg2_1_ghz1_1_JHUGen = self.mela.computeP(True)
                        self.mela.setProcess(TVar.bkgZZ, TVar.MCFM, TVar.ZZQQB)
                        p_QQB_BKG_MCFM = self.mela.computeP(True)
                        self.mela.resetInputEvent()
                    ZZ=ZZCand(Z1, Z2, p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM)
                    if DEBUG: print("ZZ=", (Z1.p4+Z2.p4).M(), p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM, ZZ.KD)
                    ZZs.append(ZZ)

            # Choose best ZZ candidate
            if len(ZZs) > 0:
                if self.bestCandByMELA : bestZZ = min(ZZs, key = cmp_to_key(bestCandByDbkgKin))
                else: bestZZ = min(ZZs, key = cmp_to_key(bestCandByZ1Z2))

                Z1 = bestZZ.Z1
                Z2 = bestZZ.Z2

                # the four leptons
                ZZ_mass = bestZZ.p4.M()

                # Compute ZZ mass without adding FSR
                zzleps = [Z1.l1, Z1.l2, Z2.l1, Z2.l2] 
                m4l_p4=ROOT.TLorentzVector()
                for i in range (0,4): m4l_p4 += (zzleps[i].p4())
                ZZ_massPreFSR=m4l_p4.M()

                ZZFlav = bestZZ.finalState() # product of the IDs of the 4 leps (Assuming SF)

                Dbkgkin = 0.
                if self.runMELA:
                    Dbkgkin = self.lib.D_bkg_kin(c_float(bestZZ.p_GG_SIG_ghg2_1_ghz1_1_JHUGen), c_float(bestZZ.p_QQB_BKG_MCFM), c_int(int(ZZFlav)),c_float(ZZ_mass)) #including c-constants

                # Sync-style printout
                if DUMP: print ('{}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{}:{:.6g}'.format(eventId, ZZ_mass, Z1.M, Z2.M, ZZ_massPreFSR, ZZFlav, Dbkgkin)) #event, ZZMass, Z1Mass, Z2Mass, ZZMassPreFSR, ZZFlav, Dbkgkin


                ### Fill output tree
                self.out.fillBranch("ZZ_mass", ZZ_mass)
                self.out.fillBranch("ZZ_massPreFSR", ZZ_massPreFSR)
                self.out.fillBranch("Z1_mass", Z1.M)
                self.out.fillBranch("Z2_mass", Z2.M)

                # Compute weights #FIXME to be moved to a separate module
                w_total = 1.

                if self.isMC: 
                    w_dataMC = 1. 
                    w_pu = 1. #FIXME This will have to me moved in a further module.
                    w_total = event.Generator_weight*w_dataMC*w_pu*self.XS
                    self.out.fillBranch("W_pu", w_pu)
                    self.out.fillBranch("W_dataMC", w_dataMC)
                # Write this for data as well for simplicity
                self.out.fillBranch("W_total", w_total)
                                
                ### Fill plots for best candidate
                self.h_ZZMass.Fill(ZZ_mass) # Mass

                return True

        return False # No candidate found, do not write event

