#!/usr/bin/env python
#
#
#

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.tools import *
from importlib import import_module
from functools import cmp_to_key

#MELA
from JHUGenMELA.MELA.mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar

import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

charge={1:"+", -1:"-"}
ZmassValue = 91.1876;

DEBUG = False
DUMP = True          # dump candidates
isMC = False         # analyze for gen info
printGenHist = False # print MC history
runMELA = False
bestCandByMELA = False # requires also runMELA=True

conf = dict(
    muPt = 5.,
    elePt = 7.,
    relIso = 0.35,
    sip3d = 4.,
    dxy =  0.5,
    dz = 1.,
    fsr_dRET2 = 0.012,
    fsr_Iso = 1.8
)

def passEleBDT(pt, SCeta, BDT) :
    fSCeta = abs(SCeta)
    # FIXME: Using 2017 WP for ElectronMVAEstimatorRun2Fall17IsoV2Values. Other BDTs are not yet available
    return (pt<=10. and ((fSCeta<0.8 and BDT > 0.85216885148) or (fSCeta>=0.8 and fSCeta<1.479 and BDT >  0.82684550976) or (fSCeta>=1.479 and BDT >  0.86937630022))) or (pt>10. and  ((fSCeta<0.8 and BDT >  0.98248928759) or (fSCeta>=0.8 and fSCeta<1.479 and BDT >  0.96919224579) or (fSCeta>=1.479 and BDT >  0.79349796445)))



### Comparators to search for the best candidate in the event 
def bestCandByZ1Z2(a,b): # a=abs(MZ1-MZ); b-sum(PT) 
    if abs(a[2]-b[2]) < 1e-4 : # same Z1: choose the candidate with highest-pT Z2 leptons
        if a[3] > b[3] : return -1        
        else :
            return 1
    else :
        if a[2] < b[2]: return -1 # else choose based on Z1 masses
        else :
            return 1

def bestCandByDbkgKin(a,b):
    if abs((a[0][4]+a[1][4]).M() - (b[0][4]+b[1][4]).M())<1e-4 and a[0][7]*a[1][7]==b[0][7]*b[1][7]: # Equivalent: same masss (tolerance 100 keV) and same FS -> same leptons and FSR
        return bestCandByZ1Z2(a,b)
    if a[4] > b [4] : return 1 # choose by best dbkgkin
    else: return -1    

# Match Rec FSR with Gen FSR 
def fsrMatch (rec, gen) :
    FSRTrue = False
    for i,ph in enumerate (gen):
        if deltaR (ph.eta, ph.phi, rec.eta, rec.phi) <0.3 and ph.pt>2.:
            FSRTrue= True
        else:
            continue
    return FSRTrue

#Find particle's mother
def Mother (part, gen) :
    idxMother= part.genPartIdxMother
    while idxMother>=0 and gen[idxMother].pdgId == part.pdgId:
        idxMother = gen[idxMother].genPartIdxMother
    idMother=0
    if idxMother >=0 : idMother = gen[idxMother].pdgId
    return idxMother, idMother


class ZZProducer(Module):
    def __init__(self):
        self.writeHistFile=True
        if runMELA :
            self.mela = Mela(13, 125, TVar.ERROR)
            self.mela.setCandidateDecayMode(TVar.CandidateDecay_ZZ)      

    def beginJob(self,histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        # Histograms
        self.h_ZZMass = ROOT.TH1F('ZZMass','ZZMass',130,70,200)
        self.addObject(self.h_ZZMass)
        
        self.h_ptmFSR = ROOT.TH1F('muon_pt_FSR', 'muon_pt_FSR',70,0.,90.)
        self.addObject(self.h_ptmFSR)
        
        # Candidate mass, with FSR
        h4e=ROOT.TH1F('4e','4e',140,60.,200.)
        h4mu=ROOT.TH1F('4mu','4mu',140,60.,200.)
        h2e2mu=ROOT.TH1F('2e2mu','2e2mu',140,60.,200.)
        self.hZZ=[h4e,h4mu,h2e2mu]
        
        # Candidate mass with FSR, FSR events only
        h4eF=ROOT.TH1F('4e FSR','4e FSR',140,60.,200.)
        h4muF=ROOT.TH1F('4mu FSR','4mu FSR',140,60.,200.)
        h2e2muF=ROOT.TH1F('2e2mu FSR','2e2mu FSR',140,60.,200.)
        self.hZZFSR=[h4eF,h4muF,h2e2muF]

        #Candidate mass without FSR, FSR events only
        title= ['no correction for FSR']
        h4e_ncF=ROOT.TH1F('4e no correction FSR','4e no correction FSR',140,60.,200.)        
        h4mu_ncF=ROOT.TH1F('4mu no correction FSR','4mu no correction FSR',140,60.,200.)        
        h2e2mu_ncF=ROOT.TH1F('2e2mu no correction FSR','2e2mu no correction FSR',140,60.,200.) 
        self.hZZ_ncF=[h4e_ncF,h4mu_ncF,h2e2mu_ncF]
        
        for i in range (0,3):
            self.addObject(self.hZZ[i])
            self.addObject(self.hZZFSR[i])
            self.addObject(self.hZZ_ncF[i])

#    def endJob(self):
#        pass

#    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#        self.out = wrappedOutputTree
#        self.out.branch("EventMass", "F")

#    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#        pass


    # Recompute lepton isolation subtracting energy of the given FSR photons
    def isoFsrCorr(self, l, selectedFSR) : 
        combRelIsoPFFSRCorr = l.pfRelIso03_all
        for f in selectedFSR :
            dR = deltaR(l.eta, l.phi, f.eta, f.phi)
            if dR >0.01 and dR < 0.3 :
                combRelIsoPFFSRCorr = max(0., l.pfRelIso03_all-f.pt/l.pt)
#                if DEBUG : print ("FSR_iso_corr,", l.pt, f.pt)
        return combRelIsoPFFSRCorr


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        eventId='{}:{}:{}'.format(event.run,event.luminosityBlock,event.event)
        if DEBUG or DUMP: print ('Event '+eventId)
        

        ### MC Truth
        genZZMass=-1.
        GenFSR = []
        if(isMC):
            genZZMass=-1.
            genpart=Collection(event,"GenPart")
            if printGenHist : 
                print ("---Gen:")
                for i, gp in enumerate(genpart) :
                    motherId=-1
                    gmotherId=-1
                    if gp.genPartIdxMother >= 0 : 
                        motherId = genpart[gp.genPartIdxMother].pdgId
                        if genpart[gp.genPartIdxMother].genPartIdxMother >= 0 :
                            gmotherId = genpart[genpart[gp.genPartIdxMother].genPartIdxMother].pdgId
                    print (i, gp.pdgId, gp.genPartIdxMother, gp.pt, gp.eta, gp.status)
            
                ##LHEPart
                LHEPart = Collection (event, "LHEPart")
                print("---------LHEPart---------")
                for i, Lp in enumerate(LHEPart):
                    print(i, Lp.pdgId, Lp.pt, Lp.eta, Lp.status, Lp.incomingpz)
              
            for i, gp in enumerate(genpart) :
                midx = gp.genPartIdxMother
                if gp.pdgId==22 and gp.pt > 1. and midx >= 0 :
                    mid = genpart[midx].pdgId
                    if abs(mid) == 11 or abs(mid) == 13 :
                        mmidx, mmid = Mother(genpart[midx], genpart) # possibly skip intermediate rows
                        if mmid == 23 :
                            GenFSR.append(gp)
                   
            
            GenHLeps = filter(lambda f : 
                              (abs(f.pdgId)==11 or abs(f.pdgId)==13 or abs(f.pdgId)==15) and
                              f.genPartIdxMother >=0 and genpart[f.genPartIdxMother].pdgId == 23 and
                              genpart[f.genPartIdxMother].genPartIdxMother >= 0 and
                              genpart[ genpart[f.genPartIdxMother].genPartIdxMother].pdgId == 25, genpart) #leptons generated from H->ZZ
            GenEl_acc = filter(lambda f: abs(f.pdgId)== 11 and abs(f.eta)<2.5 and f.pt > conf["elePt"], GenHLeps)
            GenMu_acc = filter(lambda f: abs(f.pdgId)== 13 and abs(f.eta)<2.4 and f.pt > conf["muPt"], GenHLeps)
            gen_p4 = ROOT.TLorentzVector()
            for lep in GenHLeps :
                gen_p4 += lep.p4()
            print("Gen: {:} leps M={:.4g} In Acc: {:} e, {:} mu,  {:} FSR".format(len(GenHLeps),gen_p4.M(), 
                                                                                                     len(GenEl_acc), len(GenMu_acc), len (GenFSR)), "\n")
           
        ### good PV filter
        if event.PV_npvsGood == 0 : return False


        ### Trigger (FIXME: this is for the 2018 menu!)
        passSingleEle = event.HLT_Ele32_WPTight_Gsf
        passSingleMu = event.HLT_IsoMu24
        passDiEle = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_DoubleEle25_CaloIdL_MW
        passDiMu = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
        passMuEle = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ or event.HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ
        passTriEle = False
        passTriMu = event.HLT_TripleMu_10_5_5_DZ or event.HLT_TripleMu_12_10_5

        passTrigger = passDiEle or passDiMu or passMuEle or passTriEle or passTriMu or passSingleEle or passSingleMu # FIXME: for MC; for data must be done by PD

        if not passTrigger : return False


        # Select Tight leptons
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        leps = filter(muTightId, muons) 
        leps += filter(eleTightId, electrons)
        nlep=len(leps)

        # Select FSR photons associated to muons passing loosse ID+SIP, to correct isolation
        fsrPhotons = Collection(event, "FsrPhoton")
        selectedFSR = filter(lambda f : fsrSel(f) and muLooseId(muons[f.muonIdx]) and muons[f.muonIdx].sip3d<conf["sip3d"], fsrPhotons) 
        combRelIsoPFFSRCorr = [0.]*nlep
        for i,l in enumerate(leps) :
            combRelIsoPFFSRCorr[i] = self.isoFsrCorr(l, selectedFSR)            

        ### Lepton histograms
        # Fill muon pt histogram for muons w FSR
        for i,muon in enumerate(muons):
            if(muon.fsrPhotonIdx>=0):
                self.h_ptmFSR.Fill(muon.pt)

        ### Print lepton info
        if DUMP : 
            for i,lep in enumerate(muons):
#                if not muLooseId(lep) and not DEBUG : continue # Print only leptons passing loose ID
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
                print(' combRelIsoPF={:.3g} combRelIsoPFFSRCorr={:.3g}'.format(lep.pfRelIso03_all,self.isoFsrCorr(lep, selectedFSR)),
                      end="")
                print(' isLoose={} isTight={}'.format(int(muLooseId(lep)),int(muTightId(lep))),
                      end="")
                if (lep.fsrPhotonIdx>=0):
                    fsr=fsrPhotons[lep.fsrPhotonIdx]

                    print(' FSR: pt={:.3g} eta={:.3g}, dREt2={:.3g}, Iso={:.3g}, true={}'.format(fsr.pt, fsr.eta, fsr.dROverEt2, fsr.relIso03, int(fsrMatch(fsr,GenFSR))))
                        
                else:
                    print()

            for i,lep in enumerate(electrons):
#                if not eleLooseId(lep) and not DEBUG : continue # Print only leptons passing loose ID
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
        

        ### Z combinatorial over selected leps (after FSR-corrected ISO cut for muons)
        Zs = []
        for i,l1 in enumerate(leps):
            if abs(l1.pdgId) == 13 and combRelIsoPFFSRCorr[i] >= conf["relIso"] : continue
            for j in range(i+1,nlep):
                l2 = leps[j]
                if abs(l2.pdgId) == 13 and combRelIsoPFFSRCorr[j] >= conf["relIso"] : continue
                if l1.pdgId == -l2.pdgId: #OS,SF
                    myfsr = {i:ROOT.TLorentzVector(), j:ROOT.TLorentzVector()}
                    Z_p4 = (l1.p4() + l2.p4())
                    for k in [i, j]: 
                        if abs(leps[k].pdgId)==13 and leps[k].fsrPhotonIdx>=0: # add FSR if present
                            fsr = fsrPhotons[leps[k].fsrPhotonIdx]
                            if fsrSel(fsr) :
                                fsr_p4 = ROOT.TLorentzVector() # note: fsr.p4() does not work as it relies on fsr.M, which is not stored
                                fsr_p4.SetPtEtaPhiM(fsr.pt,fsr.eta,fsr.phi,0.) 
                                myfsr[k] = fsr_p4
                                Z_p4 += fsr_p4
                    zmass = Z_p4.M()
                    if DUMP: print('Z={:.4g} pt1={:.3g} pt2={:.3g}'.format(zmass, l1.pt, l2.pt))
                    if (zmass>12. and zmass<120.):
                        Zs.append([zmass, l1.pt+l2.pt, i, j, Z_p4, myfsr[i], myfsr[j], l1.pdgId*l2.pdgId]) # mass_with_FSR, sumpt, l_i, l_j, p4, i_FSRp4, j_FSRp4, FS

        # Build all ZZ combinations passing the ZZ selection
        ZZ = []
        if len(Zs) >= 2:
            for iZ,aZ in enumerate(Zs):
                for jZ in range(iZ+1, len(Zs)):
                    # check that these Zs are mutually exclusive (not sharing the same lepton) 
                    if Zs[iZ][2]==Zs[jZ][2] or Zs[iZ][2]==Zs[jZ][3] or Zs[iZ][3]==Zs[jZ][2] or Zs[iZ][3]==Zs[jZ][3]: continue
                    #FIXME:  add DeltaR>0.02 cut among all leptons (still necessary)?
                    
                    # set Z1 and Z2
                    Z1Idx, Z2Idx = jZ, iZ
                    if abs(Zs[iZ][0]-ZmassValue) < abs(Zs[jZ][0]-ZmassValue):
                        Z1Idx, Z2Idx = Z2Idx, Z1Idx

                    Z1=Zs[Z1Idx]
                    Z2=Zs[Z2Idx]

                    # Z1 mass cut
                    if Z1[0] <= 40. : continue

                    # indexes of four leptons (Z1_1, Z1_2, Z2_1, Z2_2)
                    lIdxs  =[Z1[2],Z1[3],Z2[2],Z2[3]] 
                    # p4s of corresponding FRSs
                    fsrP4s =[Z1[5],Z1[6],Z2[5],Z2[6]] 
                    
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
                        if (abs(mZa-ZmassValue)<abs(Z1[0]-ZmassValue)) and mZb < 12.: continue

                    #Compute D_bkg^kin
                    p_GG_SIG_ghg2_1_ghz1_1_JHUGen = 0.
                    p_QQB_BKG_MCFM = 1.
                    if runMELA:
                        daughters = SimpleParticleCollection_t()
                        for ilep in range (0,4):
                            lep = leps[lIdxs[ilep]]
                            lep_p4 = lep.p4()
                            daughters.push_back(SimpleParticle_t(lep.pdgId, lep.p4() + fsrP4s[ilep]))                            
                            #print("MELA ", daughters[ilep].first, daughters[ilep].second.Pt(), daughters[ilep].second.Eta(), daughters[ilep].second.Phi(),daughters[ilep].second.M())
                        self.mela.setInputEvent(daughters, 0, 0, 0)
                        self.mela.setProcess(TVar.HSMHiggs, TVar.JHUGen, TVar.ZZGG)
                        p_GG_SIG_ghg2_1_ghz1_1_JHUGen = self.mela.computeP(True)
                        self.mela.setProcess(TVar.bkgZZ, TVar.MCFM, TVar.ZZQQB)
                        p_QQB_BKG_MCFM = self.mela.computeP(True)
                        self.mela.resetInputEvent()
                    dBkgKin = p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(p_GG_SIG_ghg2_1_ghz1_1_JHUGen+p_QQB_BKG_MCFM)

                    if DEBUG: print("ZZ=", (Z1[4]+Z2[4]).M(), p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM, dBkgKin)

                    ZZ.append([Z1, Z2, abs(Z1[0]-ZmassValue), Z2ptsum, dBkgKin]) 

            # Choose best ZZ candidate
            bestZZ=[]
            if len(ZZ) > 0:
                if bestCandByMELA : bestZZ = min(ZZ, key = cmp_to_key(bestCandByDbkgKin))
                else: bestZZ = min(ZZ, key = cmp_to_key(bestCandByZ1Z2))
                bestZZ_p4=bestZZ[0][4]+bestZZ[1][4]
                ZZmass = bestZZ_p4.M()

                # Sync-style printout
                if DUMP: print ('{}:{:.3f}:{:.3f}:{:.3f}:{:.6g}'.format(eventId, ZZmass, bestZZ[0][0], bestZZ[1][0], bestZZ[4])) #event, ZZMass, Z1Mass, Z2Mass, Dbkgkin
                
                ### Fill plots for best candidate
                self.h_ZZMass.Fill(ZZmass) # Mass

                ZZflavour = abs(leps[lIdxs[0]].pdgId * leps[lIdxs[2]].pdgId)
                plotFlavours={11**2:0, 13**2:1, 11*13:2} # 4e, 4mu, 2e2mu
                pIdx=plotFlavours[ZZflavour]
                
                self.hZZ[pIdx].Fill(ZZmass) # Mass by final state

                # Mass for events with FSR photons
                hasFSR = reduce(lambda ret,p4 : ret and p4.Pt()>0, fsrP4s, False) #check if any of the leptons has FSR
                if(hasFSR):
                    # Compute ZZ mass without adding FSR
                    m4l_p4=ROOT.TLorentzVector()
                    for i in range (0,4):
                        m4l_p4 += (leps[lIdxs[i]].p4())
                    ZZmass_lep=m4l_p4.M()

                    self.hZZFSR[pIdx].Fill(ZZmass)
                    self.hZZ_ncF[pIdx].Fill(ZZmass_lep)                



        return True



muLooseId = lambda l : l.pt > conf["muPt"] and abs(l.eta) < 2.4 and abs(l.dxy) < conf["dxy"] and abs(l.dz) < conf["dz"] and (l.isGlobal or l.isTracker) #FIXME: l.nStations>0) is not equivalent to numberOfMatches>0); also, muonBestTrackType!=2 is not available
eleLooseId = lambda l : l.pt > conf["elePt"] and abs(l.eta) < 2.5 and abs(l.dxy) < conf["dxy"] and abs(l.dz) < conf["dz"]

muTightId = lambda l : muLooseId(l) and (l.isPFcand or (l.highPtId>0 and l.pt>200.)) and abs(l.sip3d) < conf["sip3d"] #FIXME: highPtId does not match exactly our definition.
    
eleTightId =  lambda l : eleLooseId(l) and abs(l.sip3d) < conf["sip3d"] and passEleBDT(l.pt, l.eta+l.deltaEtaSC, l.mvaFall17V2Iso) #FIXME: We use different BDTs for 2016 and 2018

fsrSel = lambda f : f.dROverEt2 < conf["fsr_dRET2"] and f.relIso03 < conf["fsr_Iso"] and f.pt > 2. and abs(f.eta) < 2.4

ZZSequence = [ZZProducer()]

preselection = ("nMuon + nElectron >= 2 &&" +
                "Sum$(Muon_pt > {muPt}) +" +
                "Sum$(Electron_pt > {elePt})" +
                ">= 2").format(**conf)

# Select specific events to debug
#preselection = ("event.eventId == 840922")

localstore = "/eos/cms/"
aaastore   = "root://cms-xrd-global.cern.ch/"
# ggH125
files = [localstore+"/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root"]

# ZH125
#files = [aaastore+"/store/mc/RunIIAutumn18NanoAODv7/ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/C670B342-C02E-1848-AE6A-0B5550E3DFE3.root"]

#High mass sample
#files = [aaastore+"/store/mc/RunIIAutumn18NanoAODv7/VBF_HToZZTo4L_M3000_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/A82C31DF-0A6F-B44F-857B-89BFD72CEFEA.root"]

p = PostProcessor(".", files, cut=preselection, branchsel=None, modules=ZZSequence, noOut=True, histFileName="histOut.root", histDirName="plots", maxEntries=0)
p.run()
