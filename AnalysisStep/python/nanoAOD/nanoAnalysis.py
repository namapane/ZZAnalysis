#!/usr/bin/env python
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.tools import *
from importlib import import_module
from functools import cmp_to_key
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

charge={1:"+", -1:"-"}
ZmassValue = 91.1876;

DEBUG = True

conf = dict(
    muPt = 5.,
    elePt = 7.,
    relIso = 0.35,
    sip3d = 4.,
    dxy =  0.5,
    dz = 1.,
    eleId = "mvaFall17V2noIso_WPL",
    muId = ""
)

# Comparator to search for the best candidate in the event 
def bestCandByZ1Z2(a,b):
    if abs(a[2]-b[2]) < 1e-4: # same Z1: choose the candidate with highest-pT Z2 leptons
        if a[3] > b[3] : return 1
        else: return -1
    else:
        if a[2] < b[2]: return 1 # else choose based on Z1 masses
        else: return -1



class ZZProducer(Module):
    def __init__(self, muSel=None, eleSel=None):
        self.muSel = muSel
        self.eleSel = eleSel
      
#        self.branches = ["Z"]
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

#    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#        self.out = wrappedOutputTree
#        self.out.branch("EventMass", "F")

#    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")


        leps = filter(self.muSel, muons)
        leps += filter(self.eleSel, electrons)
        nlep=len(leps)


        if DEBUG: print ('Event {}:{}:{}'.format(event.run,event.luminosityBlock,event.event))


        # Select FSR photons associated to muons passing loosse ID +SIP, to correct isolation
        selectedFSR = filter(lambda f : muLooseId(muons[f.muonIdx]) and muons[f.muonIdx].sip3d<conf["sip3d"], fsrPhotons) 
        combRelIsoPFFSRCorr = [0.]*nlep
        for i,l in enumerate(leps) :
            combRelIsoPFFSRCorr[i] = l.pfRelIso03_all
            for f in selectedFSR :
                dR = deltaR(l.eta, l.phi, f.eta, f.phi)
                if dR >0.01 and dR < 0.3 :
                    combRelIsoPFFSRCorr[i] = max(0., l.pfRelIso03_all-f.pt/l.pt)
                    if DEBUG : print ("FSR_iso_corr,", l.pt, f.pt) 
                    
                                                 
        if DEBUG: # Print lepton info
            for i,lep in enumerate(muons):

                print('mu{}  pt={:.3g} eta={:.3g} phi={:.3g} dxy={:.3g} dz={:.3g} SIP={:.3g}'.format(charge[lep.charge],
                                                                                                     lep.pt,
                                                                                                     lep.eta,
                                                                                                     lep.phi,
                                                                                                     abs(lep.dxy),
                                                                                                     abs(lep.dz),
                                                                                                     lep.sip3d),
                      end="")
                print(' GLB={} TK={} isPF={} isHighPtId={}'.format(int(lep.isGlobal),
                                                                   int(lep.isTracker),
                                                                   int(lep.isPFcand),
                                                                   int(lep.highPtId>0),# 1=trk,2=global
                                                                   ),
                      end="")
                print(' combRelIsoPF={:.3g} combRelIsoPFFSRCorr={:.3g}'.format(lep.pfRelIso03_all,0.), #FIXME combRelIsoPFFSRCorr[i]
                      end="")
                print(' isLoose={} isTight={}'.format(int(muLooseId(lep)),int(muTightId(lep))),
                      end="")
                if (lep.fsrPhotonIdx>=0):
                    fsr=fsrPhotons[lep.fsrPhotonIdx]
                    print(' FSRpt={:.3g}'.format(fsr.pt))
                else:
                    print()

        

        # Z combinatorial over selected leps + ISO for muons
        Zs = []
        for i,l1 in enumerate(leps):
            if abs(l1.pdgId) == 13 and combRelIsoPFFSRCorr[i] >= conf["relIso"] : continue
            for j in range(i+1,nlep):
                l2 = leps[j]
                if abs(l2.pdgId) == 13 and combRelIsoPFFSRCorr[j] >= conf["relIso"] : continue
                if l1.pdgId == -l2.pdgId:
                    Z_p4 = (l1.p4() + l2.p4())
                    for k in [i, j]: 
                        if abs(leps[k].pdgId)==13 and leps[k].fsrPhotonIdx>=0: # add FSR. FIXME more compact way to do the same?
                            fsr = fsrPhotons[leps[k].fsrPhotonIdx]
                            fsr_p4 = ROOT.TLorentzVector()
                            fsr_p4.SetPtEtaPhiM(fsr.pt,fsr.eta,fsr.phi,0.) 
                            Z_p4 += fsr_p4
                    zmass = Z_p4.M()

                    print('Z={:.4g} pt1={:.3g} pt2={:.3g}'.format(zmass, l1.pt, l2.pt))
                    if (zmass>12. and zmass<120.):
                        Zs.append([zmass, l1.pt+l2.pt, i, j, Z_p4]) # mass, sumpt, l_i, l_j, p4

        # Build all ZZ combinations passing the ZZ selection
        ZZ = []
        if len(Zs) >= 2:
            for i,aZ in enumerate(Zs):
                for j in range(i+1, len(Zs)):
                    # check that these Zs are mutually exclusive
                    if Zs[i][2]==Zs[j][2] or Zs[i][2]==Zs[j][3] or Zs[i][3]==Zs[j][2] or Zs[i][3]==Zs[j][3]: continue
                    #FIXME:  add DeltaR>0.02 cut among all leptons (still required)?

                    
                    # set Z1 and Z2
                    Z1Idx=j
                    Z2Idx=i
                    if abs(Zs[i][0]-ZmassValue) < abs(Zs[j][0]-ZmassValue):
                        Z1Idx, Z2Idx = Z2Idx, Z1Idx

                    # Z1 mass cut
                    if Zs[Z1Idx][0] <= 40. : continue

                    lIdxs=[Zs[Z1Idx][2],Zs[Z1Idx][3],Zs[Z2Idx][2],Zs[Z2Idx][3]]

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
                        if (abs(mZa-ZmassValue)<abs(Zs[Z1Idx][0]-ZmassValue)) and mZb < 12.: continue
                    
                    ZZ.append([Z1Idx, Z2Idx, abs(Zs[Z1Idx][0]-ZmassValue), Z2ptsum])
                    #DEBUG print ([Z1Idx, Z2Idx, abs(Zs[Z1Idx][0]-ZmassValue), Z2ptsum])
            # Choose best ZZ candidate. FIXME: implement selection by D_bkg^kin
            if len(ZZ) > 0:
                bestZZ = min(ZZ, key = cmp_to_key(bestCandByZ1Z2))
                bestZZ_p4=Zs[bestZZ[0]][4]+Zs[bestZZ[1]][4]
                print ('{}:{}:{}:{:.4g}:{:.3g}:{:.3g}:'.format(event.run,event.luminosityBlock,event.event,
                                                               bestZZ_p4.M(),Zs[ZZ[0][0]][0], Zs[ZZ[0][1]][0]))


        return True


class lepSkim(Module): #This is just an example; non actually used
    def __init__(self):
        self.writeHistFile = True

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        self.h_vpt = ROOT.TH1F('sumpt', 'sumpt', 100, 0, 1000)
        self.addObject(self.h_vpt)

    def analyze(self, event):
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = Collection(event, "Jet")
        fsrPhotons = Collection(event, "FsrPhoton")
        eventSum = ROOT.TLorentzVector()
        
        # select events with at least 2 muons
        if len(muons) >= 2:
            for lep in muons:  # loop on muons
                eventSum += lep.p4()
            for lep in electrons:  # loop on electrons
                eventSum += lep.p4()
#             for j in jets:  # loop on jets
#                 eventSum += j.p4()
            self.h_vpt.Fill(eventSum.M())  # fill histogram

        return True


muLooseId = lambda l : l.pt > conf["muPt"] and abs(l.eta) < 2.4 and abs(l.dxy) < conf["dxy"] and abs(l.dz) < conf["dz"] and (l.isGlobal or l.isTracker) #FIXME: isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) && muonBestTrackType!=2
eleLooseId = lambda l : l.pt > conf["elePt"] and abs(l.eta) < 2.5 and abs(l.dxy) < conf["dxy"] and abs(l.dz) < conf["dz"] #

muTightId = lambda l : muLooseId(l) and (l.isPFcand or (l.highPtId>0 and l.pt>200.)) and abs(l.sip3d) < conf["sip3d"] #FIXME: highPtId does not match exactly our definituion.
    
eleTightId =  lambda l : eleLooseId(l) and abs(l.sip3d) < conf["sip3d"] #FIXME add BDT

ZZSequence = [lepSkim(), ZZProducer(muSel = muTightId, eleSel = eleTightId)]

#preselection = "Jet_pt[0] > 250"
preselection = ("nMuon + nElectron >= 2 &&" + 
                "Sum$(Muon_pt > {muPt}) +" + #&& Muon_miniPFRelIso_all < {miniRelIso} && Muon_sip3d < {sip3d}"
                "Sum$(Electron_pt > {elePt})" + #&& Electron_miniPFRelIso_all < {miniRelIso} && Electron_sip3d < {sip3d} && Electron_{eleId}"
                ">= 2").format(**conf)

store = "/eos/cms/" #if not at CERN: "root://cms-xrd-global.cern.ch/"
files = [store+"/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root"]

p = PostProcessor(".", files, cut=preselection, branchsel=None, modules=ZZSequence, noOut=True, histFileName="histOut.root", histDirName="plots", maxEntries=100)
p.run()
