# Plot files produced with FSRStandaloneAnalyzer
# Extract FSR candidates with:
# grep -e "^FSR:" out_FSR.txt -a | sed -e s/FSR:// -e s/:/" "/g > FSR_tree.txt
# zgrep -ah "^FSR:" AAAOK/ggH125_Chunk*/log.txt.gz | sed -e s/FSR:// -e s/:/" "/g > FSR_tree.txt
from __future__ import print_function
import math
from ROOT import TAttLine,TCanvas, TFile, TColor, gStyle, TH1F, TH1, gROOT, TNtuple, gDirectory
from ROOT import kBlack, kBlue, kRed, kOrange, kGreen, kPink,kCyan

def eff(num, den) :
    e = num/den
    return e, math.sqrt(e*(1-e)/den)

def printCanvases(type="png", path=".") :
    canvases = gROOT.GetListOfCanvases()
    for c in canvases :
        c.Print(path+"/"+c.GetTitle()+"."+type)

gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
#gStyle.SetHistLineWidth(4)
gStyle.SetPadGridX(True);
gStyle.SetPadGridY(True);
gStyle.SetGridColor(0);
gStyle.SetGridStyle(3);
gStyle.SetGridWidth(1);

TH1.SetDefaultSumw2()


# SCVetoTight : True = would be vetoed by tighter SC veto (as in nanoAOD)
# vetoPT :   pT of the electron causing a SC tight veto
# lID :ID of the reco l the photon is associated to
# lPT :TpT of the reco l the photon is associated to
# lTight  : lepton matched to gamma passes tight selection
# isLoose : gamma passes nanoAOD cuts (preselected->always true)
# isTight : gamma passes H4l cuts
# maskLoose : masked by another photon passing loose sel
# maskTight : masked by another photon passing tight sel
# isSelStd  : selected by standard sel: (isTight&&!maskTight)
# isSelNano : tight sel applied after matching loose photons to leptons (emulate tight on top of nanoAOD)
# muEleOverlap : photon can be associated (within the max DR cut) to both a mu and an ele
# stealer_ID : of loose lepton that would steal this photon (be careful of tight veto)
# stealer_pT
# stealer_DR
# pTGen : pT of the gen FSR photon matched to the reco photon. Only FSR photons from W,Z->l are considered.
# genFS: final state (4mu=0, 4e=1, 2e2mu=2)
# genMother: 25 = l from H->Z; 23 = Z, etc.


#n=TNtuple("test", "test", "run:ls:event:pT:eta:phi:gRelIso:SCVetoTight:lID:lTight:dR:dRET2:isLoose:isTight:maskLoose:maskTight:isSelStd:isSelNano:muEleOverlap:pTGen:etaGen:phiGen:genFS:genMother")
n=TNtuple("test", "test", "run:ls:event:pT:eta:phi:gRelIso:SCVetoTight:vetoPT:lID:lpT:lTight:dR:dRET2:isLoose:isTight:maskLoose:maskTight:isSelStd:isSelNano:muEleOverlap:stealer_ID:stealer_pT:stealer_DR:pTGen:etaGen:phiGen:genFS:genMother")
n.ReadFile("/afs/cern.ch/user/n/namapane/public/for_Valentina/FSR_tree5.txt");


selection = "isSelStd&&lTight&&genFS<3"
#selection = "isSelStd&&lTight&&genFS<3&&genMother==25" # only for matched lepton from H
#selection = "isSelStd&&lTight&&(!SCVetoTight)" # Add tighter SC veto

c01 =TCanvas ("SCveto","SC veto",800,800)

p_ele = n.Draw("pT",selection+"&&abs(lID)==11&&pTGen>0")/float(n.Draw("pT",selection+"&&abs(lID)==11"))
p_mu  = n.Draw("pT",selection+"&&abs(lID)==13&&pTGen>0")/float(n.Draw("pT",selection+"&&abs(lID)==13"))

print("purity e/mu: ", p_ele, p_mu)

# Effect of tight SC veto
n.Draw("pTGen>0:SCVetoTight>>h(2,0,2,2,0,2)",selection,"textcolz")
h = gROOT.FindObject('h')
h.SetXTitle("Tight SC veto")
h.SetYTitle("Fake/True")
fakesRej=h.GetBinContent(2,1)
truesRej=h.GetBinContent(2,2)
allFakes=h.GetBinContent(1,1)+fakesRej
allTrues=h.GetBinContent(1,2)+truesRej
deltaTrues=truesRej/allTrues
deltaPur=fakesRej/(allTrues+allFakes)
print("tight SC veto: Delta_eff = ", deltaTrues, "Delta_fakes", fakesRej/allFakes, "Delta_purity", deltaPur)

c02=TCanvas ("stealing","stealing",1600,800)
c02.Divide(2,1)
c02.cd(1)
n.Draw("pTGen>0:stealer_pT>0>>h2(2,0,2,2,0,2)",selection,"textcolz")
h2 = gROOT.FindObject('h2')
h2.SetXTitle("Stealed False/True")
h2.SetYTitle("Fake/True")
c02.cd(2)
n.Draw("pTGen>0:stealer_pT>0>>h3(2,0,2,2,0,2)",selection+"&&(!SCVetoTight)","textcolz")

c03 =TCanvas ("Masking","Masking",800,800)
n.Draw("pTGen>0:maskLoose>>h4(2,0,2,2,0,2)",selection,"textcolz")

# Discrepaancy in nanoAOD due to stealing
# >>> (121.+93)/(11857+3654+121+93)
#     0.013608903020667727
# For trues and fakes:
# >>> (121.)/(11857+121)
#     0.010101853397896142
# >>> (93.)/(3654+93)
#     0.024819855884707767


nbins_pt = 16
nMAX_pt=80.
nmin_pt = 0.

nbins_eta = 30
nMAX_eta=3.
nmin_eta = 0.

c1 = TCanvas("ctemp","ctemp",800,800)


h_pur_mu_pt=[TH1F("den_pt_mu","den_pt_mu",nbins_pt,nmin_pt,nMAX_pt),
             TH1F("num_pt_mu","num_pt_mu",nbins_pt,nmin_pt,nMAX_pt),
             TH1F("purity_pt_mu","purity_pt_mu",nbins_pt,nmin_pt,nMAX_pt)]
h_pur_el_pt =[TH1F("den_pt_el","den_pt_el",nbins_pt,nmin_pt,nMAX_pt),
              TH1F("num_pt_el","num_pt_el",nbins_pt,nmin_pt,nMAX_pt),
              TH1F("purity_pt_el","purity_pt_el",nbins_pt,nmin_pt,nMAX_pt)]

h_pur_mu_eta=[TH1F("den_eta_mu","den_eta_mu",nbins_eta,nmin_eta,nMAX_eta),
              TH1F("num_eta_mu","num_eta_mu",nbins_eta,nmin_eta,nMAX_eta),
              TH1F("purity_eta_mu","purity_eta_mu",nbins_eta,nmin_eta,nMAX_eta)]
h_pur_el_eta =[TH1F("den_eta_el","den_eta_el",nbins_eta,nmin_eta,nMAX_eta),
               TH1F("num_eta_el","num_eta_el",nbins_eta,nmin_eta,nMAX_eta),
               TH1F("purity_eta_el","purity_eta_el",nbins_eta,nmin_eta,nMAX_eta)]

p_ele = n.Draw("pT",selection+"&&abs(lID)==11&&pTGen>0")/float(n.Draw("pT",selection+"&&abs(lID)==11"))
p_mu  = n.Draw("pT",selection+"&&abs(lID)==13&&pTGen>0")/float(n.Draw("pT",selection+"&&abs(lID)==13"))

c1.Delete()


iEntry=0
nEntries = n.GetEntries()
while iEntry< nEntries and n.GetEntry(iEntry):
    if( n.isSelStd and n.lTight and n.genFS<3 and abs(n.lID)== 13):
        h_pur_mu_pt[0].Fill(n.pT)
        h_pur_mu_eta[0].Fill(abs(n.eta))

        if(n.pTGen>0):
           h_pur_mu_pt[1].Fill(n.pT)
           h_pur_mu_eta[1].Fill(abs(n.eta))
           
    if( n.isSelStd and n.lTight and n.genFS<3 and abs(n.lID)== 11):
        h_pur_el_pt[0].Fill(n.pT)
        h_pur_el_eta[0].Fill(abs(n.eta))
        if(n.pTGen>0):
            h_pur_el_pt[1].Fill(n.pT)
            h_pur_el_eta[1].Fill(abs(n.eta))  
    
    iEntry+=1

h_pur_mu_pt[2].Divide(h_pur_mu_pt[1],h_pur_mu_pt[0],1.,1.,"b")
h_pur_mu_pt[2].GetXaxis().SetTitle("pT [GeV]")
c1= TCanvas("p_mu_pt","p_mu_pt",800,800)
h_pur_mu_pt[2].Draw()

h_pur_mu_eta[2].Divide(h_pur_mu_eta[1],h_pur_mu_eta[0],1.,1.,"b")
h_pur_mu_eta[2].GetXaxis().SetTitle("#eta")
c2= TCanvas("p_mu_eta","p_mu_eta",800,800)
h_pur_mu_eta[2].Draw()

h_pur_el_pt[2].Divide(h_pur_el_pt[1],h_pur_el_pt[0],1.,1.,"b")
h_pur_el_pt[2].GetXaxis().SetTitle("pT [GeV]")
c3= TCanvas("p_el_pt","p_el_pt",800,800)
h_pur_el_pt[2].Draw()

h_pur_el_eta[2].Divide(h_pur_el_eta[1],h_pur_el_eta[0],1.,1.,"b")
h_pur_el_eta[2].GetXaxis().SetTitle("#eta")
c4= TCanvas("p_el_eta","p_el_eta",800,800)
h_pur_el_eta[2].Draw()




