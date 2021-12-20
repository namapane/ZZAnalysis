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

### Comparators to search for the best candidate in the event 
def bestCandByZ1Z2(a,b): # a=abs(MZ1-MZ); b-sum(PT) 
    if abs(a.Z1.M-b.Z1.M) < 1e-4 : # same Z1: choose the candidate with highest-pT Z2 leptons
        if a.Z2.sumpt > b.Z2.sumpt : return -1        
        else :
            return 1
    else :
        if abs(a.Z1.M-ZmassValue) < abs(b.Z1.M-ZmassValue): return -1 # else choose based on Z1 masses
        else :
            return 1

#FIXME
# def bestCandByDbkgKin(a,b):
#     if abs((a[0][4]+a[1][4]).M() - (b[0][4]+b[1][4]).M())<1e-4 and a[0][7]*a[1][7]==b[0][7]*b[1][7]: # Equivalent: same masss (tolerance 100 keV) and same FS -> same leptons and FSR
#         return bestCandByZ1Z2(a,b)
#     if a[4] > b [4] : return 1 # choose by best dbkgkin
#     else: return -1    

def bestCandByDbkgKin(a,b):
    if abs((a.p4).M() - (b.p4).M())<1e-4 and a.Z1.FS*a.Z2.FS==b.Z1.FS*b.Z2.FS: # Equivalent: same masss (tolerance 100 keV) and same FS -> same leptons and FSR
        return bestCandByZ1Z2(a,b)
    if a.KD > b.KD : return 1 # choose by best dbkgkin
    else: return -1    


# Struct to cache computed values for Z candidates.
# We may need to move this into the Event with a separate ZFiller module at some point.
class ZCand: 
    def __init__(self, p4, l1Idx, l2Idx, fsr1Idx, fsr2Idx, sumpt, FS):
        self.p4 = p4
        self.l1Idx = l1Idx
        self.l2Idx = l2Idx
        self.fsr1Idx = fsr1Idx
        self.fsr2Idx = fsr2Idx
        self.sumpt = sumpt
        self.FS=FS
        self.M = p4.M()

# Struct to cache computed values for ZZ candidates.
class ZZCand:
    def __init__(self, Z1, Z2, p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM):
        self.Z1 = Z1
        self.Z2 = Z2
        self.p4 = Z1.p4+Z2.p4
        self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen = p_GG_SIG_ghg2_1_ghz1_1_JHUGen
        self.p_QQB_BKG_MCFM = p_QQB_BKG_MCFM
        self.KD = p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(p_GG_SIG_ghg2_1_ghz1_1_JHUGen+p_QQB_BKG_MCFM) # without c-constants, for candidate sorting

class ZZProducer(Module):
    def __init__(self, runMELA, bestCandByMELA, XS, cuts):
        self.writeHistFile=True
        self.runMELA=runMELA
        self.bestCandByMELA = bestCandByMELA
        self.XS = XS
        self.muLooseId=cuts["muLooseId"]
        self.eleLooseId=cuts["eleLooseId"]

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


#    def endJob(self):
#        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("ZZ_mass", "F")
        self.out.branch("ZZ_massPreFSR", "F")
        self.out.branch("Z1_mass", "F")
        self.out.branch("Z2_mass", "F")
        self.out.branch("W_pu", "F")
        self.out.branch("W_dataMC", "F")
        self.out.branch("W_total", "F")


        
        # gen weight = Generator_weight

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        eventId='{}:{}:{}'.format(event.run,event.luminosityBlock,event.event)
        if DEBUG : print ('Event '+eventId)

        # Collections
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        leps = filter(self.muLooseId, muons)
        leps += filter(self.eleLooseId, electrons) # muons+electrons passing loose ID (without SIP)
        nlep=len(leps)
        

        ### Z combinatorial over selected leps (after FSR-corrected ISO cut for muons)
        Zs = []
        for i,l1 in enumerate(leps):
            if not l1.isTightIso : continue
            fsr1_p4 = fsrPhotons[l1.myFsrPhotonIdx].p4() if l1.myFsrPhotonIdx>=0 else ROOT.TLorentzVector()
                
            for j in range(i+1,nlep):
                l2 = leps[j]
                if not l2.isTightIso : continue
                if l1.pdgId == -l2.pdgId: #OS,SF
                    fsr2_p4 = fsrPhotons[l2.myFsrPhotonIdx].p4() if l2.myFsrPhotonIdx>=0 else ROOT.TLorentzVector()
                    Z_p4 = l1.p4() + l2.p4() + fsr1_p4 + fsr2_p4
                    zmass = Z_p4.M()
                    if DEBUG: print('Z={:.4g} pt1={:.3g} pt2={:.3g}'.format(zmass, l1.pt, l2.pt))
                    if (zmass>12. and zmass<120.):
                        # FIXME: we may want to set a default order for i, j.
                        Zs.append(ZCand(Z_p4, i, j, l1.myFsrPhotonIdx,  l2.myFsrPhotonIdx, l1.pt+l2.pt, l1.pdgId*l2.pdgId))
# AT some point, we may want to store Z candidates in the event, in a dedicated ZFiller module; like:
#                         Z_d1[iZ]=i
#                         Z_d2[iZ]=j
#                         Z_fsr1Idx[iZ]=l1.myFsrPhotonIdx
#                         Z_fsr2Idx[iZ]=l2.myFsrPhotonIdx
#                         Z_pt[iZ]=Z_p4.pt()
#                         Z_eta[iZ]=Z_p4.eta()
#                         Z_phi[iZ]=Z_p4.phi()
#                         Z_mass[iZ]=Z_p4.M()
#                         Z_finalState[iZ]=l1.pdgId*l2.pdgId


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

                    # indexes of four leptons (Z1_1, Z1_2, Z2_1, Z2_2)
                    lIdxs  =[Z1.l1Idx,Z1.l2Idx,Z2.l1Idx,Z2.l2Idx]
                    # indexes of corresponding FRSs
                    fsrIdxs =[Z1.fsr1Idx,Z1.fsr2Idx,Z2.fsr1Idx,Z2.fsr2Idx] 
                    
                    lepPts =[]
                    # QCD suppression on all OS pairs, regardelss of flavour
                    passQCD = True
                    for k in range(0,4):
                        lepPts.append(leps[lIdxs[k]].pt)
                        for l in range (k+1,4):
                            if leps[lIdxs[k]].pdgId*leps[lIdxs[l]].pdgId < 0 and (leps[lIdxs[k]].p4()+leps[lIdxs[l]].p4()).M()<=4.: passQCD = False
                    if not passQCD: continue

                    Z2ptsum=lepPts[2]+lepPts[3]

                    # trigger acceptance cuts (20,10 GeV)
                    lepPts.sort()
                    if not (lepPts[3]>20. and lepPts[2]>10.) : continue

                    #"Smart cut" on alternate pairings for same-sign candidates
                    if abs(leps[lIdxs[0]].pdgId) == abs(leps[lIdxs[2]].pdgId):
                        mZa, mZb = 0., 0.
                        if leps[lIdxs[0]].pdgId == -leps[lIdxs[2]].pdgId:
                            mZa=(leps[lIdxs[0]].p4()+leps[lIdxs[2]].p4()).M()
                            mZb=(leps[lIdxs[1]].p4()+leps[lIdxs[3]].p4()).M()
                        elif leps[lIdxs[0]].pdgId == -leps[lIdxs[3]].pdgId:
                            mZa=(leps[lIdxs[0]].p4()+leps[lIdxs[3]].p4()).M()
                            mZb=(leps[lIdxs[1]].p4()+leps[lIdxs[2]].p4()).M()
                        if (abs(mZa-ZmassValue)>abs(mZb-ZmassValue)) : mZa, mZb = mZb, mZa
                        if (abs(mZa-ZmassValue)<abs(Z1.M-ZmassValue)) and mZb < 12.: continue

                    #Compute D_bkg^kin
                    p_GG_SIG_ghg2_1_ghz1_1_JHUGen = 0.
                    p_QQB_BKG_MCFM = 1.
                    if self.runMELA:
                        daughters = SimpleParticleCollection_t()
                        for ilep in range (0,4):
                            lep = leps[lIdxs[ilep]]
                            lep_p4 = lep.p4()
                            if fsrIdxs[ilep]>=0 : lep_p4 += fsrPhotons[fsrIdxs[ilep]].p4
                            daughters.push_back(SimpleParticle_t(lep.pdgId, lep.p4()))                         
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
#                    ZZ.append([Z1, Z2, abs(Z1.M-ZmassValue), Z2ptsum, KD, p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM]) 

            # Choose best ZZ candidate
            bestZZ=[]
            if len(ZZs) > 0:
                if self.bestCandByMELA : bestZZ = min(ZZs, key = cmp_to_key(bestCandByDbkgKin))
                else: bestZZ = min(ZZs, key = cmp_to_key(bestCandByZ1Z2))

                Z1 = bestZZ.Z1
                Z2 = bestZZ.Z2

                # indexes of four leptons (Z1_1, Z1_2, Z2_1, Z2_2)
                lIdxs  =[Z1.l1Idx,Z1.l2Idx,Z2.l1Idx,Z2.l2Idx]
                # indixes of corresponding FRSs
                fsrIdxs =[Z1.fsr1Idx,Z1.fsr2Idx,Z2.fsr1Idx,Z2.fsr2Idx] 
 
                bestZZ_p4=Z1.p4+Z2.p4
                ZZ_mass = bestZZ_p4.M()

                # Compute ZZ mass without adding FSR
                m4l_p4=ROOT.TLorentzVector()
                for i in range (0,4): m4l_p4 += (leps[lIdxs[i]].p4())
                ZZ_massPreFSR=m4l_p4.M()


                ZZFlav = (leps[Z1.l1Idx].pdgId * leps[Z2.l2Idx].pdgId)**2 # product of the IDs of the 4 leps (Assuming SF)

                Dbkgkin = 0.
                if self.runMELA:
                    Dbkgkin = self.lib.D_bkg_kin(c_float(bestZZ.p_GG_SIG_ghg2_1_ghz1_1_JHUGen), c_float(bestZZ.p_QQB_BKG_MCFM), c_int(int(ZZFlav)),c_float(ZZ_mass)) #including c-constants

                # Sync-style printout
                if DUMP: print ('{}:{:.3f}:{:.3f}:{:.3f}:{:.3f}:{}:{:.6g}'.format(eventId, ZZ_mass, Z1.M, Z2.M, ZZ_massPreFSR, ZZFlav, Dbkgkin)) #event, ZZMass, Z1Mass, Z2Mass, ZZMassPreFSR, ZZFlav, Dbkgkin


                # Compute weights
                w_dataMC = 1. #FIXME
                w_pu = 1. #FIXME This will have to me moved in a further module.
                w_total = event.Generator_weight*w_dataMC*w_pu*self.XS

                ### Fill output tree
                self.out.fillBranch("ZZ_mass", ZZ_mass)
                self.out.fillBranch("ZZ_massPreFSR", ZZ_massPreFSR)
                self.out.fillBranch("Z1_mass", Z1.M)
                self.out.fillBranch("Z2_mass", Z2.M)
                self.out.fillBranch("W_pu", w_pu)
                self.out.fillBranch("W_dataMC", w_dataMC)
                self.out.fillBranch("W_total", w_total)
                                
                ### Fill plots for best candidate
                self.h_ZZMass.Fill(ZZ_mass) # Mass

                return True

        return False # No candidate found, do not write event

