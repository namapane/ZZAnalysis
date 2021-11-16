# Produce FSR-related plots from CJLST trees.
# For some of the plots, extended information is expected; trees must have been processed with addFSRDetails = true. 
#
from __future__ import print_function
import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

extendedTree = True # Trees include extended FSR info (addFSRDetails = true)
nEvents = 1e9

def binError(num, den) :
    eff=float(num)/den
    return math.sqrt((eff*(1-eff))/den)

def printEff(num, den):
    return  "({:.3f}+{:.3f})%".format(100.*num/den,100.*binError(num,den))

def loop():

    ### official production
    #inFileName = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MC_2018/ggH125/ZZ4lAnalysis.root" 

    ### removing FSR for e and mu < 20 GeV as in nanoAOD
    #inFileName = "/eos/user/n/namapane/HZZ/FSR/ZZ4lAnalysis-emulateNanoFSR.root" 
    #outFileName = "fsrHisto_nano.root"

    ### private trees; std sel; with extended FSR info (full ggH125)
    inFileName = "/eos/user/n/namapane/HZZ/FSR/ZZ4lAnalysis.root" 
    outFileName= "fsrHisto.root"


    print ("Processing file: ",inFileName,"...")

    tree = ROOT.TChain("ZZTree/candTree")
    tree.Add(inFileName)
    of=ROOT.TFile(outFileName,"recreate")
    

    h_fsrMuPt=ROOT.TH1F('fsrMuPt', 'fsrMuPt', 40, 0., 200.)
    h_fsrElePt=ROOT.TH1F('fsrElePt', 'fsrElePt', 40, 0., 200.)
    h_fsrMuPtVsFakePt=ROOT.TH2F('fsrMuPtVsFakePt', 'fsrMuPtVsFakePt', 200, 0., 200., 100, 0., 100.)
    h_fsrElePtVsFakePt=ROOT.TH2F('fsrElePtVsFakePt', 'fsrElePtVsFakePt', 200, 0., 200., 100, 0., 100.)
    h_FSR_ET2vsDR_true=ROOT.TH2F('fsrET2vsDR_true', 'fsrET2vsDR_true', 100, 0., 400, 100, 0., 0.5)
    h_FSR_ET2vsDR_fake=ROOT.TH2F('fsrET2vsDR_fake', 'fsrET2vsDR_fake', 100, 0., 400, 100, 0., 0.5)

    #purity
    h_fsrMuPtFake=ROOT.TH1F("fsrMuPtFake","fsrMuPtFake",200,0.,200.)
    h_fsrMuPtTrue=ROOT.TH1F("fsrMuPtTrue","fsrMuPtTrue",40,0.,200.)
    h_fsrElePtFake=ROOT.TH1F("fsrElePtFake","fsrElePtFake",200,0.,200.)
    h_fsrElePtTrue=ROOT.TH1F("fsrElePtTrue","fsrElePtTrue",40,0.,200.)
    
    #purity vs eta
    h_MufsrEta= ROOT.TH1F("MufsrEta","MufsrEta",25,0.,2.5)
    h_MufsrEtaTrue= ROOT.TH1F("MufsrEtaTrue","MufsrEtaTrue",25,0.,2.5)


    minM=40
    maxM=200
    nbins=160
    h_ZZMassFSR=[ROOT.TH1F('ZZMassFSR_4e', 'ZZMassFSR_4e', nbins, minM, maxM),
                 ROOT.TH1F('ZZMassFSR_4mu', 'ZZMassFSR_4mu', nbins, minM, maxM),
                 ROOT.TH1F('ZZMassFSR_2e2mu', 'ZZMassFSR_2e2mu', nbins, minM, maxM)]
    h_ZZMassNoFSR=[ROOT.TH1F('ZZMassNoFSR4e', 'ZZMassNoFSR4e', nbins, minM, maxM),
                   ROOT.TH1F('ZZMassNoFSR4mu', 'ZZMassNoFSR4mu', nbins, minM, maxM),
                   ROOT.TH1F('ZZMassNoFSR2e2mu', 'ZZMassNoFSR2e2mu', nbins, minM, maxM)] #ZZ mass without correction for FSR
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
    nCandRecoverZ1Mass = [0]*3 #4e, 4mu, 2e2mu
    nCandRecoverZ2Mass = [0]*3 #4e, 4mu, 2e2mu
    nLepWithFSR=0
    nLepFSRRecovered_self=0
    nLepFSRRecovered_other=0
    nEleFSRRecovered_self= [0,0] 
    nMuFSRRecovered_self=[0,0]
    nEleFSRRecovered_other= [0,0] #[true,fakeq] 
    nMuFSRRecovered_other=[0,0] # [true, fake]
    
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
        
        #Window mass
        lep_p4 = ROOT.TLorentzVector()
        Z1_p4 = ROOT.TLorentzVector()
        Z2_p4 = ROOT.TLorentzVector()
        for jLep,l in enumerate (tree.LepPt):
            lep_p4.SetPtEtaPhiM(tree.LepPt[jLep],tree.LepEta[jLep],tree.LepPhi[jLep],0.)
            if jLep <2 : Z1_p4 += lep_p4
            else: Z2_p4 += lep_p4
        ZZ_mass= (Z1_p4+Z2_p4).M()
        #Check ZZMass 
        if(tree.ZZMassPreFSR-ZZ_mass > 0.01):
            print("ERROR Mass")
        if(Z1_p4.M()<40 or Z1_p4.M()>120):
            nCandRecoverZ1Mass[finalState]+=1
        if(Z2_p4.M()<12 or Z2_p4.M()>120):
            nCandRecoverZ2Mass[finalState]+=1


        if len(tree.fsrPt) > 0 :
            nLepWithFSR+=len(tree.fsrPt)
            hasFSR = True # FSR associated to one of the leptons             
            nCandWFSR[finalState]+=1 # number of candidates with FSR
                
#            print (tree.EventNumber, tree.ZZMass, tree.ZZMassPreFSR, tree.Z1Flav, tree.Z2Flav, tree.LepPt[0], tree.passIsoPreFSR)
            h_ZZMassFSR[finalState].Fill(tree.ZZMass)
            h_ZZMassNoFSR[finalState].Fill(tree.ZZMassPreFSR)
            for i,f in enumerate(tree.fsrPt):                
                lepIdx=tree.fsrLept[i]-1
#                print (" ",tree.fsrPt[i], tree.fsrLept[i], tree.LepPt[lepIdx], tree.LepLepId[lepIdx])
                if abs(tree.LepLepId[lepIdx]) == 13 : 
                    h_fsrMuPt.Fill(tree.LepPt[lepIdx])
                    h_MufsrEta.Fill(abs(tree.fsrEta[i])) #eta
                elif abs(tree.LepLepId[lepIdx]) == 11 : h_fsrElePt.Fill(tree.LepPt[lepIdx])                
                if (extendedTree) :
                    ET2=tree.fsrPt[i]*tree.fsrPt[i]
                    isFake = tree.fsrGenPt[i]<0.
#                    print ("FSR", tree.RunNumber, tree.LumiNumber, tree.EventNumber, tree.fsrPt[i], tree.fsrDR[i]/ET2, isFake)# tree.fsrLept[i], tree.LepPt[lepIdx], tree.LepLepId[lepIdx])
                    
                    if isFake :
                        hasFake = True
                        h_FSR_ET2vsDR_fake.Fill(ET2,tree.fsrDR[i])
                        #Plot of lepton pT:photon pT for fakes
                        if abs(tree.LepLepId[lepIdx]) == 13 : 
                            h_fsrMuPtVsFakePt.Fill(tree.LepPt[lepIdx], f)
                            h_fsrMuPtFake.Fill(tree.LepPt[lepIdx]) 
                        elif abs(tree.LepLepId[lepIdx]) == 11 : 
                            h_fsrElePtVsFakePt.Fill(tree.LepPt[lepIdx], f)     
                            h_fsrElePtFake.Fill(tree.LepPt[lepIdx])
     

                    else:
                        h_FSR_ET2vsDR_true.Fill(ET2,tree.fsrDR[i])
                        if abs(tree.LepLepId[lepIdx]) == 13 :
                            h_fsrMuPtTrue.Fill(tree.LepPt[lepIdx])
                            h_MufsrEtaTrue.Fill(abs(tree.fsrEta[i]))
                        else:
                            h_fsrElePtTrue.Fill(tree.LepPt[lepIdx])

            #Check how often a lepton is recovered by its own or another photon
            for j in range(0,4):
                if not bool(tree.LepIsoPreFSR[j]) : # lepton was recovered by FSR
                    ownFSR = False
                    for f in (tree.fsrLept):                        
                        if j == f-1 : # lepton has FSR attached (we then assume that FSR caused the lepton's recovery)
                            nLepFSRRecovered_self+=1 
                            if abs(tree.LepLepId[j]) == 11: 
                                if hasFake: # Event contains a fake FSR not necessarily this one!
                                    nEleFSRRecovered_self[1] += 1
                                else:
                                    nEleFSRRecovered_self[0] += 1

                            elif abs(tree.LepLepId[j]) == 13:
                                if hasFake:
                                    nMuFSRRecovered_self[1] += 1
                                else:
                                    nMuFSRRecovered_self[0] += 1
                            ownFSR = True

                    if not ownFSR : 
                        nLepFSRRecovered_other+=1
                        if abs(tree.LepLepId[j]) == 11: 
                                if hasFake:
                                    nEleFSRRecovered_other[1] += 1
                                else:
                                    nEleFSRRecovered_other[0] +=1
                        elif abs(tree.LepLepId[j]) == 13:
                                if hasFake:
                                    nMuFSRRecovered_other[1] += 1
                                else:
                                    nMuFSRRecovered_other[0] +=1 
                    

            
            if hasFake :
                h_ZZMassFSRFake[finalState].Fill(tree.ZZMass)
                nhasFake[finalState]+=1
            else :
                h_ZZMassFSRTrue[finalState].Fill(tree.ZZMass)




        if pass_selection and not tree.passIsoPreFSR : # not passIsoPreFSR = code will be executed when Isolation failed without FSR
            nIsoRecovered[finalState]+=1
            if not hasFSR : print (" FSR-ISO-RECOVER") # do we recover candidates due to iso correction of FSR from (eg loose) leptons that are not part of the candidate? Apparently, not.        
            if hasFake : nIsoRecovered_wFake[finalState]+=1

        # Consistency checks
        if hasFSR != (abs(tree.ZZMassPreFSR-tree.ZZMass) > 1e-4) : # strange if there's FSR and mass difference is small. Check this situation
            print ("Inconsistent FSR: ", hasFSR, tree.ZZMassPreFSR, tree.ZZMass)

        if (extendedTree) : 
            chkPassIsoPreFSR = True
            for i in range(0,4):
                chkPassIsoPreFSR = chkPassIsoPreFSR and bool(tree.LepIsoPreFSR[i]) 
            if chkPassIsoPreFSR != tree.passIsoPreFSR :
                print ("Error passIsoPreFSR: ", chkPassIsoPreFSR, tree.passIsoPreFSR)

                
    of.Write()

#     nCandSel  = [22420, 43183, 55968]
#     nCandWFSR = [730, 2332, 2391]
#     nhasFake  =  [155, 391, 458]
#     nIsoRecovered = [0, 987, 601]
#     nIsoRecovered_fake = [0, 79, 48]    
    print ("nCandSel =   ", sum(nCandSel), nCandSel) 
    print ("nCandWFSR =  {:} [4e, 4mu, 2e2mu]: {:} [{} {} {}]".format( sum(nCandWFSR), nCandWFSR, printEff(nCandWFSR[0],nCandSel[0]), printEff(nCandWFSR[1],nCandSel[1]), printEff(nCandWFSR[2],nCandSel[2])))
    print ("nhasFake  =  {:} [4e, 4mu, 2e2mu]: {:} [{:.1f}% {:.1f}% {:.1f}%]".format(sum(nhasFake), nhasFake, 100.*nhasFake[0]/nCandSel[0],100.*nhasFake[1]/nCandSel[1],100.*nhasFake[2]/nCandSel[2]))
    print ("Ev purity  =  [4e, 4mu, 2e2mu]: [{} {} {}]".format(printEff(nCandWFSR[0]-nhasFake[0],nCandWFSR[0]),printEff(nCandWFSR[1]-nhasFake[1],nCandWFSR[1]),printEff(nCandWFSR[2]-nhasFake[2],nCandWFSR[2])))
    print ("nRecovered = {:} [4e, 4mu, 2e2mu]: {:} [{:.1f}% {:.1f}% {:.1f}%]".format(sum(nIsoRecovered), nIsoRecovered, 100.*nIsoRecovered[0]/nCandSel[0],100.*nIsoRecovered[1]/nCandSel[1],100.*nIsoRecovered[2]/nCandSel[2]))
    print ("nRecovered_fake = {:} [4e, 4mu, 2e2mu]: {:} [{:.1f}% {:.1f}% {:.1f}%]".format(sum(nIsoRecovered_wFake), nIsoRecovered_wFake, 100.*nIsoRecovered_wFake[0]/nCandSel[0],100.*nIsoRecovered_wFake[1]/nCandSel[1],100.*nIsoRecovered_wFake[2]/nCandSel[2]))
    print ("Yield var (iso)  =  [4e, 4mu, 2e2mu]: [{} {} {}]".format(printEff(nIsoRecovered[0],nCandSel[0]-nIsoRecovered[0]),printEff(nIsoRecovered[1],nCandSel[1]-nIsoRecovered[1]),printEff(nIsoRecovered[2],nCandSel[2]-nIsoRecovered[2])))

    print ("Total number of leptons with FSR   = ", nLepWithFSR) 
    print ("Lepton ISO recovered by own photon = ", nLepFSRRecovered_self)
    print ("Lepton ISO recovered by other phot = ", nLepFSRRecovered_other)
    print ("splitted by events containing all true FSR or containing fakes:")
    print ("Muon ISO recovered by own photon_true =", nMuFSRRecovered_self[0])
    print ("Muon ISO recovered by own photon_fake = ",nMuFSRRecovered_self[1])
    print ("Muon ISO recovered by oth photon_true = ", nMuFSRRecovered_other[0])
    print ("Muon ISO recovered by oth photon_fake = ", nMuFSRRecovered_other[1])
    #Just to check
    print("Ele ISO recovered other photon_true = ", nEleFSRRecovered_other[0])
    print("Ele ISO recovered other photon_fake = ", nEleFSRRecovered_other[1])
    
    #Mass Window
    print(" nCandRecoverZ1Mass =  [4e, 4mu, 2e2mu]: {:} [{} {} {}] ".format(nCandRecoverZ1Mass,printEff(nCandRecoverZ1Mass[0],nCandSel[0]), printEff(nCandRecoverZ1Mass[1],nCandSel[1]), printEff(nCandRecoverZ1Mass[2],nCandSel[2])))
    print(" nCandRecoverZ2Mass =  [4e, 4mu, 2e2mu]: {:} [{} {} {}] ".format(nCandRecoverZ2Mass,printEff(nCandRecoverZ2Mass[0],nCandSel[0]), printEff(nCandRecoverZ2Mass[1],nCandSel[1]), printEff(nCandRecoverZ2Mass[2],nCandSel[2])))
    print(" nCandRecoverZZMass =  [4e, 4mu, 2e2mu]: [{} {} {}]".format(printEff(nCandRecoverZ1Mass[0]+nCandRecoverZ2Mass[0],nCandSel[0]), printEff(nCandRecoverZ1Mass[1]+nCandRecoverZ2Mass[1],nCandSel[1]), printEff(nCandRecoverZ1Mass[2]+nCandRecoverZ2Mass[2],nCandSel[2])))

if __name__ == "__main__" : #Check me 
    loop()
