# Produce FSR-related plots from CJLST trees.
# For some of the plots, extended information is expected; trees must have been processed with addFSRDetails = true. 
#
from __future__ import print_function
import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

extendedTree = True # Trees include extended FSR info (addFSRDetails = true)
nEvents = 2e9

def loop():
#    inFileName = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MC_2018/ggH125/ZZ4lAnalysis.root" #official trees
    inFileName = "/eos/user/n/namapane/HZZ/FSR/ZZ4lAnalysis.root" # private trees with extended info
    outFileName = "fsrHisto.root"

    print ("Processing file: ",inFileName,"...")

    tree = ROOT.TChain("ZZTree/candTree")
    tree.Add(inFileName)
    of=ROOT.TFile(outFileName,"recreate")
    

    h_fsrMuPt=ROOT.TH1F('fsrMuPt', 'fsrMuPt', 200, 0., 200.)
    h_fsrElePt=ROOT.TH1F('fsrElePt', 'fsrElePt', 200, 0., 200.)
    h_fsrMuPtVsFakePt=ROOT.TH2F('fsrMuPtVsFakePt', 'fsrMuPtVsFakePt', 200, 0., 200., 100, 0., 100.)
    h_fsrElePtVsFakePt=ROOT.TH2F('fsrElePtVsFakePt', 'fsrElePtVsFakePt', 200, 0., 200., 100, 0., 100.)
    h_FSR_ET2vsDR_true=ROOT.TH2F('fsrET2vsDR_true', 'fsrET2vsDR_true', 100, 0., 400, 100, 0., 0.5)
    h_FSR_ET2vsDR_fake=ROOT.TH2F('fsrET2vsDR_fake', 'fsrET2vsDR_fake', 100, 0., 400, 100, 0., 0.5)


    minM=40
    maxM=200
    nbins=160
    h_ZZMassFSR=[ROOT.TH1F('ZZMassFSR_4e', 'ZZMassFSR_4e', nbins, minM, maxM),
                 ROOT.TH1F('ZZMassFSR_4mu', 'ZZMassFSR_4mu', nbins, minM, maxM),
                 ROOT.TH1F('ZZMassFSR_2e2mu', 'ZZMassFSR_2e2mu', nbins, minM, maxM)]
    h_ZZMassNoFSR=[ROOT.TH1F('ZZMassNoFSR4e', 'ZZMassNoFSR4e', nbins, minM, maxM),
                   ROOT.TH1F('ZZMassNoFSR4mu', 'ZZMassNoFSR4mu', nbins, minM, maxM),
                   ROOT.TH1F('ZZMassNoFSR2e2mu', 'ZZMassNoFSR2e2mu', nbins, minM, maxM)]
    h_ZZMassFSRFake=[ROOT.TH1F('ZZMassFSRFake4e', 'ZZMassFSRFake4e', nbins, minM, maxM),
                     ROOT.TH1F('ZZMassFSRFake4mu', 'ZZMassFSRFake4mu', nbins, minM, maxM),
                     ROOT.TH1F('ZZMassFSRFake2e2mu', 'ZZMassFSRFake2e2mu', nbins, minM, maxM)]
    h_ZZMassFSRTrue=[ROOT.TH1F('ZZMassFSRTrue4e', 'ZZMassFSRTrue4e', nbins, minM, maxM),
                     ROOT.TH1F('ZZMassFSRTrue4mu', 'ZZMassFSRTrue4mu', nbins, minM, maxM),
                     ROOT.TH1F('ZZMassFSRTrue2e2mu', 'ZZMassFSRTrue2e2mu', nbins, minM, maxM)]



    nCandSel=[0]*3 #4e, 4mu, 2e2mu
    nCandWFSR=[0]*3 #4e, 4mu, 2e2mu
    nIsoRecovered=[0]*3 #4e, 4mu, 2e2mu
    nIsoRecovered_wFake=[0]*3 #4e, 4mu, 2e2mu
    nhasFake=[0]*3 #4e, 4mu, 2e2mu
    nLepFSRRecovered_self=0
    nLepFSRRecovered_other=0

    iEntry=0
    while iEntry<nEvents and tree.GetEntry(iEntry):
        iEntry+=1
        if iEntry%1000 == 0 : print("Processing", iEntry)

        pass_selection = False
        if tree.ZZsel>=90 : pass_selection = True

        if not pass_selection:
            print (tree.EventNumber, "Fail sel")
            continue
        
        finalState=-1
        if tree.Z1Flav*tree.Z2Flav == 11*11*11*11 : finalState = 0 #4e
        elif tree.Z1Flav*tree.Z2Flav == 13*13*13*13 : finalState = 1 #4mu
        elif tree.Z1Flav*tree.Z2Flav == 11*11*13*13 : finalState = 2 #2e2mu
        else : print("ERROR!!! ZFlav")
        
        nCandSel[finalState]+=1

        # Select events affected by FSR
        hasFSR = False
        hasFake = False
        nRecoveredLeps = [0,0] #true, false 
        if len(tree.fsrPt) > 0 :
            hasFSR = True # FSR associated to one of the leptions             
            nCandWFSR[finalState]+=1
                
#            print (tree.EventNumber, tree.ZZMass, tree.ZZMassPreFSR, tree.Z1Flav, tree.Z2Flav, tree.LepPt[0], tree.passIsoPreFSR)
            h_ZZMassFSR[finalState].Fill(tree.ZZMass)
            h_ZZMassNoFSR[finalState].Fill(tree.ZZMassPreFSR)
            for i,f in enumerate(tree.fsrPt):                
                lepIdx=tree.fsrLept[i]-1
#                print (" ",tree.fsrPt[i], tree.fsrLept[i], tree.LepPt[lepIdx], tree.LepLepId[lepIdx])
                if abs(tree.LepLepId[lepIdx]) == 13 : h_fsrMuPt.Fill(tree.LepPt[lepIdx])
                elif abs(tree.LepLepId[lepIdx]) == 11 : h_fsrElePt.Fill(tree.LepPt[lepIdx])                
                if (extendedTree) :
                    ET2=tree.fsrPt[i]*tree.fsrPt[i]
                    isFake = tree.fsrGenPt[i]<0.
                    if isFake :
                        hasFake = True
                        h_FSR_ET2vsDR_fake.Fill(ET2,tree.fsrDR[i])
                        #Plot of lepton pT:photon pT for fakes
                        if abs(tree.LepLepId[lepIdx]) == 13 : h_fsrMuPtVsFakePt.Fill(tree.LepPt[lepIdx], f)
                        elif abs(tree.LepLepId[lepIdx]) == 11 : h_fsrElePtVsFakePt.Fill(tree.LepPt[lepIdx], f)                

                    else:
                        h_FSR_ET2vsDR_true.Fill(ET2,tree.fsrDR[i])

            #Check how often a lepton is recovered by its own or another photon
            for j in range(0,4):
                if not bool(tree.LepIsoPreFSR[j]) : # lepton was recovered by FSR
                    ownFSR = False;
                    for f in (tree.fsrLept):                        
                        if j == f-1 :
                            nLepFSRRecovered_self+=1
                            ownFSR = True
                    if not ownFSR : nLepFSRRecovered_other+=1

            if hasFake :
                h_ZZMassFSRFake[finalState].Fill(tree.ZZMass)
                nhasFake[finalState]+=1
            else :
                h_ZZMassFSRTrue[finalState].Fill(tree.ZZMass)



        if pass_selection and not tree.passIsoPreFSR :
            nIsoRecovered[finalState]+=1
            if not hasFSR : print (" FSR-ISO-RECOVER") # do we recover candidates due to iso correction of FSR from (eg loose) leptons that are not part of the candidate? Apparently, not.        
            if hasFake : nIsoRecovered_wFake[finalState]+=1

        # Consistency checks
        if hasFSR != (abs(tree.ZZMassPreFSR-tree.ZZMass) > 1e-4) : 
            print ("Inconsistent FSR: ", hasFSR, tree.ZZMassPreFSR, tree.ZZMass)

        if (extendedTree) : 
            chkPassIsoPreFSR = True
            for i in range(0,4):
                chkPassIsoPreFSR = chkPassIsoPreFSR and bool(tree.LepIsoPreFSR[i])
            if chkPassIsoPreFSR != tree.passIsoPreFSR :
                print ("Error passIsoPreFSR: ", chkPassIsoPreFSR, tree.passIsoPreFSR)

                
    of.Write()    
    print ("nCandSel =   ", sum(nCandSel), nCandSel) 
    print ("nCandWFSR =  {:} [4e, 4mu, 2e2mu]: {:} [{:.1f}% {:.1f}% {:.1f}%]".format( sum(nCandWFSR), nCandWFSR, 100.*nCandWFSR[0]/nCandSel[0], 100.*nCandWFSR[1]/nCandSel[1], 100.*nCandWFSR[2]/nCandSel[2]))
    print ("nhasFake  =  {:} [4e, 4mu, 2e2mu]: {:} [{:.1f}% {:.1f}% {:.1f}%]".format(sum(nhasFake), nhasFake, 100.*nhasFake[0]/nCandSel[0],100.*nhasFake[1]/nCandSel[1],100.*nhasFake[2]/nCandSel[2]))
    print ("nRecovered = {:} [4e, 4mu, 2e2mu]: {:} [{:.1f}% {:.1f}% {:.1f}%]".format(sum(nIsoRecovered), nIsoRecovered, 100.*nIsoRecovered[0]/nCandSel[0],100.*nIsoRecovered[1]/nCandSel[1],100.*nIsoRecovered[2]/nCandSel[2]))
    print ("nRecovered_fake = {:} [4e, 4mu, 2e2mu]: {:} [{:.1f}% {:.1f}% {:.1f}%]".format(sum(nIsoRecovered_wFake), nIsoRecovered_wFake, 100.*nIsoRecovered_wFake[0]/nCandSel[0],100.*nIsoRecovered_wFake[1]/nCandSel[1],100.*nIsoRecovered_wFake[2]/nCandSel[2]))
    print ("Lepton ISO recovered by own photon = ", nLepFSRRecovered_self) #FIXME true/fake
    print ("Lepton ISO recovered by other phot = ", nLepFSRRecovered_other)


if __name__ == "__main__" :
    loop()
