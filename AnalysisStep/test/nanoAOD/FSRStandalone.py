###
# Plot files produced with FSRStandaloneAnalyzer.cc, which runs on miniAODs and prints 1 line for each FSR candidate (no event selection)
#
# To extract files from FSRStandaloneAnalyzer.cc's output:
# grep -e "^FSR:" out_FSR.txt -a | sed -e s/FSR:// -e s/:/" "/g > FSR_tree.txt
# zgrep -ah "^FSR:" AAAOK/ggH125_Chunk*/log.txt.gz | sed -e s/FSR:// -e s/:/" "/g > FSR_tree.txt
###
from __future__ import print_function
import math
from array import array
from ROOT import TAttLine,TCanvas, TFile, TColor, gStyle, TH1F, TH1, gROOT, TNtuple, gDirectory, TLegend, TF1, TH2F, TMarker, TLine, TGraph
from ROOT import kBlack, kBlue, kRed, kOrange, kGreen, kPink, kCyan, kMagenta, kGreen, kTeal
import numpy as np

#compute efficiency and its error
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
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridColor(0)
gStyle.SetGridStyle(3)
gStyle.SetGridWidth(1)
gStyle.SetPadLeftMargin(0.11)

TH1.SetDefaultSumw2()

### Tree format:
# pT, eta, phi, gRelIso: of the photon
# SCVetoTight : True = would be vetoed by tighter SC veto (as in nanoAOD)
# vetoPT :   pT of the electron causing a SC tight veto
# lID :ID of the reco l the photon is associated to
# lPT :TpT of the reco l the photon is associated to
# lTight  : the lepton matched to the photon passes tight ZZ selection
# isLoose : gamma passes nanoAOD cuts (would be written in nanoAODs)
# isTight : gamma passes H4l cuts (tighter than isLosse)
# maskLoose  :Would not be selected because there's another photon passing loose sel associated to the same muon 
# maskTight : Would not be selected because there's another photon passing tight sel associated to the same muon 
# isSelStd  : selected by standard ZZ selection: equivalent to (isTight&&!maskTight)
# isSelNano : emulate requiring tight selection on top of nanoAOD without re-matching photons to leptons
# muEleOverlap : photon can be associated (within the max DR cut) to both a mu and an ele
# stealer_ID, pT, DR : of loose lepton that would steal this photon (be careful of tight veto)
# pTGen : pT of the gen FSR photon matched to the reco photon. Only FSR photons from W,Z->l are considered.
# genFS: final state (4mu=0, 4e=1, 2e2mu=2)
# genMother: 25 = l from H->Z; 23 = Z, etc.



n=TNtuple("test", "test", "run:ls:event:pT:eta:phi:gRelIso:SCVetoTight:vetoPT:lID:lpT:lTight:dR:dRET2:isLoose:isTight:maskLoose:maskTight:isSelStd:isSelNano:muEleOverlap:stealer_ID:stealer_pT:stealer_DR:pTGen:etaGen:phiGen:genFS:genMother")

n.ReadFile("/afs/cern.ch/user/n/namapane/public/for_Valentina/FSR_tree_nocuts.txt") # This version includes all photons with g_pt>1; g_iso<2, no cut on dr/et2


### Select photon candidates to match as much as possible those we have in actual 4l events:
# lTight - we only look at photons attached to good leptons (passing full ID)
# genMother==25: only leptons that are matched to a H decay (avoid fake leptons and leptons from taus, which normally don't get selected in 4l events)
# genFS<3: only in 2e2mu/4e/4mu events at gen level (to avoid adding backround photons from tau decays)

selection_stdveto = "isSelStd&&lTight&&genFS<3&&genMother==25" # Consider only photons in 2e2mu/4e/4mu events, matched to a lepton from H decay passing tight cuts
selection = "isSelStd && (!SCVetoTight) && lTight && genFS<3 && genMother==25" # same, with tighter SC veto


### Print overall purity
c1 = TCanvas("ctemp","ctemp",800,800)
p_ele = n.Draw("pT",selection+"&&abs(lID)==11&&pTGen>0")/float(n.Draw("pT",selection+"&&abs(lID)==11"))
p_mu  = n.Draw("pT",selection+"&&abs(lID)==13&&pTGen>0")/float(n.Draw("pT",selection+"&&abs(lID)==13"))
print("purity e/mu: ", p_ele, p_mu)
c1.Delete()



### Effect of tight SC veto 
c01 =TCanvas ("SCveto","SC veto",800,800)
n.Draw("pTGen>0:SCVetoTight>>h(2,0,2,2,0,2)",selection_stdveto,"textcolz")
h = gROOT.FindObject('h')
h.SetXTitle("Tight SC veto")
h.SetYTitle("Fake/True")
fakesRej=h.GetBinContent(2,1)
truesRej=h.GetBinContent(2,2)
allFakes=h.GetBinContent(1,1)+fakesRej
allTrues=h.GetBinContent(1,2)+truesRej
deltaTrues=truesRej/allTrues
deltaPur=fakesRej/(allTrues+allFakes)
print("tight SC veto: Delta_eff = ", deltaTrues, "Delta_fakes = ", fakesRej/allFakes, "Delta_purity = ", deltaPur)


### Effects of masking and stealing (which can be corrected by re-associating photons to leptons)
if False: 
    c02=TCanvas ("stealing","stealing",1600,800)
    c02.Divide(2,1)
    c02.cd(1)
    n.Draw("pTGen>0:stealer_pT>0>>h2(2,0,2,2,0,2)",selection,"textcolz")
    h2 = gROOT.FindObject('h2')
    h2.SetXTitle("Stolen False/True")
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


### Plots of FSR purity vspT, eta
nbins_FSR_pt = 10
#nMAX_FSR_pt=250.
#nmin_FSR_pt = 0.
lowEdge= [2,10,20,30,40,50,60,80,100,150,250]

nbins_pt = 16
nMAX_pt=80.
nmin_pt = 0.

nbins_eta = 25
nMAX_eta=2.5
nmin_eta = 0.

h_FSRpT = TH1F("FSRpT","FSRpT",nbins_FSR_pt,array('d',lowEdge))
h_pur_mu_FSR_pt=[TH1F("den_FSR_pt_mu","den_FSR_pt_mu",nbins_FSR_pt,array('d',lowEdge)),
             TH1F("num_FSR_pt_mu","num_FSR_pt_mu",nbins_FSR_pt,array('d',lowEdge)),
             TH1F("purity_FSR_pt_mu","purity_FSR_pt_mu",nbins_FSR_pt,array('d',lowEdge))]
h_pur_el_FSR_pt =[TH1F("den_FSR_pt_el","den_FSR_pt_el",nbins_FSR_pt,array('d',lowEdge)),
              TH1F("num_FSR_pt_el","num_FSR_pt_el",nbins_FSR_pt,array('d',lowEdge)),
              TH1F("purity_FSR_pt_el","purity_FSR_pt_el",nbins_FSR_pt,array('d',lowEdge))]

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

iEntry=0
nEntries = n.GetEntries()
while iEntry< nEntries and n.GetEntry(iEntry):
    if( n.isSelStd and n.lTight and n.genFS<3 and n.genMother==25 and abs(n.lID)== 13):
        h_pur_mu_FSR_pt[0].Fill(n.pT)
        h_pur_mu_pt[0].Fill(n.lpT)
        h_pur_mu_eta[0].Fill(abs(n.eta))

        if(n.pTGen>0):
           h_pur_mu_FSR_pt[1].Fill(n.pT)
           h_pur_mu_pt[1].Fill(n.lpT)
           h_pur_mu_eta[1].Fill(abs(n.eta))
           
    if( n.isSelStd and n.lTight and n.genFS<3 and n.genMother==25 and abs(n.lID)== 11):
        h_pur_el_FSR_pt[0].Fill(n.pT)
        h_pur_el_pt[0].Fill(n.lpT)
        h_pur_el_eta[0].Fill(abs(n.eta))
        if(n.pTGen>0):
            h_pur_el_FSR_pt[1].Fill(n.pT)
            h_pur_el_pt[1].Fill(n.lpT)
            h_pur_el_eta[1].Fill(abs(n.eta))
    if( n.isSelStd and n.lTight and n.genFS<3):
        h_FSRpT.Fill(n.pT)
    
    iEntry+=1

for bin in range(1,11):
    h_FSRpT.SetBinContent(bin,h_FSRpT.GetBinContent(bin)/h_FSRpT.GetBinWidth(bin))
    h_FSRpT.SetBinError(bin,h_FSRpT.GetBinError(bin)/h_FSRpT.GetBinWidth(bin))

h_FSRpT.GetXaxis().SetTitle("pT_{#gamma} [GeV]")
h_FSRpT.GetYaxis().SetTitle("Events/GeV")
c00 =TCanvas ("FSRpT", "FSRpT", 800,800)
h_FSRpT.SetMarkerStyle(20)
h_FSRpT.Draw()

h_pur_mu_FSR_pt[2].Divide(h_pur_mu_FSR_pt[1],h_pur_mu_FSR_pt[0],1.,1.,"b")
h_pur_mu_FSR_pt[2].GetXaxis().SetTitle("pT_{#gamma} [GeV]")
c1= TCanvas("p_mu_FSR_pt","p_mu_FSR_pt",800,800)
h_pur_mu_FSR_pt[2].GetYaxis().SetTitle("Purity")
h_pur_mu_FSR_pt[2].SetMarkerStyle(20)
h_pur_mu_FSR_pt[2].GetYaxis().SetRangeUser(0,1)
h_pur_mu_FSR_pt[2].Draw()

h_pur_mu_eta[2].Divide(h_pur_mu_eta[1],h_pur_mu_eta[0],1.,1.,"b")
h_pur_mu_eta[2].GetXaxis().SetTitle("#eta_{#gamma}")
c2= TCanvas("p_mu_eta","p_mu_eta",800,800)
h_pur_mu_eta[2].GetYaxis().SetTitle("Purity")
h_pur_mu_eta[2].SetMarkerStyle(20)
h_pur_mu_eta[2].GetYaxis().SetRangeUser(0,1)
h_pur_mu_eta[2].Draw()

h_pur_el_FSR_pt[2].Divide(h_pur_el_FSR_pt[1],h_pur_el_FSR_pt[0],1.,1.,"b")
h_pur_el_FSR_pt[2].GetXaxis().SetTitle("pT_{#gamma} [GeV]")
c3= TCanvas("p_el_FSR_pt","p_el_FSR_pt",800,800)
h_pur_el_FSR_pt[2].GetYaxis().SetTitle("Purity")
h_pur_el_FSR_pt[2].SetMarkerStyle(20)
h_pur_el_FSR_pt[2].GetYaxis().SetRangeUser(0,1)
h_pur_el_FSR_pt[2].Draw()

h_pur_el_eta[2].Divide(h_pur_el_eta[1],h_pur_el_eta[0],1.,1.,"b")
h_pur_el_eta[2].GetXaxis().SetTitle("#eta_{#gamma}")
c4= TCanvas("p_el_eta","p_el_eta",800,800)
h_pur_el_eta[2].GetYaxis().SetTitle("Purity")
h_pur_el_eta[2].SetMarkerStyle(20)
h_pur_el_eta[2].GetYaxis().SetRangeUser(0,1)
h_pur_el_eta[2].Draw()

h_pur_mu_pt[2].Divide(h_pur_mu_pt[1],h_pur_mu_pt[0],1.,1.,"b")
#h_pur_mu_pt[2].GetXaxis().SetTitle("pT_{#mu} [GeV]")
h_pur_mu_pt[2].GetXaxis().SetTitle("pT_{lepton} [GeV]")
c5= TCanvas("p_mu_pt","p_mu_pt",800,800)
h_pur_mu_pt[2].GetYaxis().SetTitle("Purity")
h_pur_mu_pt[2].SetMarkerStyle(20)
h_pur_mu_pt[2].GetYaxis().SetRangeUser(0,1)
h_pur_mu_pt[2].Draw()

h_pur_el_pt[2].Divide(h_pur_el_pt[1],h_pur_el_pt[0],1.,1.,"b")
#h_pur_el_pt[2].GetXaxis().SetTitle("pT_{e} [GeV]")
#c6= TCanvas("p_el_pt","p_el_pt",800,800)
h_pur_el_pt[2].GetYaxis().SetTitle("Purity")
h_pur_el_pt[2].SetMarkerStyle(4)
h_pur_el_pt[2].GetYaxis().SetRangeUser(0,1)
h_pur_el_pt[2].SetLineColor(kRed)
h_pur_el_pt[2].Draw("same")

legend = TLegend()
legend.AddEntry(h_pur_mu_pt[2], "muon FSR purity", "lp")
legend.AddEntry(h_pur_el_pt[2], "electron FSR purity", "lp")
legend.Draw()

h_pur_mu_pt[1].GetXaxis().SetTitle("pT_{lepton} [GeV]")
c7= TCanvas("lepton pT true FSR","lepton pT true FSR",800,800)
c7.SetLeftMargin(1.4)
h_pur_mu_pt[1].GetYaxis().SetTitle("Entries")
h_pur_mu_pt[1].SetMarkerStyle(20)
h_pur_mu_pt[1].Draw()

h_pur_el_pt[1].SetMarkerStyle(4)
h_pur_el_pt[1].SetLineColor(kRed)
h_pur_el_pt[1].Draw("same")

legend1 = TLegend()
legend1.AddEntry(h_pur_mu_pt[1], "muon pT (true FSR)", "lp")
legend1.AddEntry(h_pur_el_pt[1], "electron pT (true FSR)", "lp")
legend1.Draw()


### Scatter plots: Et vs DRET2
fsel = "abs(eta)<2.5 && gRelIso<1.8 && (!SCVetoTight) && lTight&&genFS<3&&genMother==25" #Photon selection (without dret2 cut and check on masking from other photons)
fakeCol=kPink-9
trueCol=kBlue+1

c10 = TCanvas ("FSR_Et_DRET2", "FSR_Et_DRET2", 800,800)
c10.SetLogx()
n.Draw("dR:pT>>hFSRFake(20,1,100,20,0,0.5)","pTGen<0&&"+fsel)
n.Draw("dR:pT>>hFSRTrue(20,1,100,20,0,0.5)","pTGen>0&&"+fsel,"same")
c10.GetListOfPrimitives()[1].GetYaxis().SetTitle("#Delta R")
c10.GetListOfPrimitives()[1].GetXaxis().SetTitle("E_{T} [GeV]   ")
hFSRFake = c10.GetListOfPrimitives()[2]
hFSRTrue = c10.GetListOfPrimitives()[3]
hFSRFake.SetMarkerStyle(8)
hFSRFake.SetMarkerSize(0.18)
hFSRFake.SetMarkerColor(fakeCol)
hFSRTrue.SetMarkerStyle(8)
hFSRTrue.SetMarkerSize(0.18)
hFSRTrue.SetMarkerColor(trueCol)
hFSRFake.GetYaxis().SetTitle("#Delta R")
hFSRFake.GetXaxis().SetTitle("E_{t} [GeV]")
hFSRFake.Draw("P")
hFSRTrue.Draw("SAMEP")


f = TF1("f1","x*x*0.012",0,100)
f.SetLineColor(kBlack)
f.Draw("SAME")

l = TLine(2.,0.,2.,0.5)
l.SetLineColor(kBlack)
l.SetLineWidth(2)
l.SetLineStyle(9)
l.Draw("SAME")


mf=TMarker(0,0,8)
mf.SetMarkerColor(fakeCol)
mf.SetMarkerSize(2)
mt=TMarker(0,0,8)
mt.SetMarkerColor(trueCol)
mt.SetMarkerSize(2)

legend2 = TLegend()
legend2.SetMargin(0.25)
#legend2.SetTextFont(63)
legend2.SetTextSize(0.03)
legend2.AddEntry(mt, "True FSR", "p")
legend2.AddEntry(mf, "Fake FSR", "p")
legend2.AddEntry(f, "#DeltaR/E^{2}_{t}<0.012","l")
legend2.Draw()


c11 = TCanvas ("FSR_Et_DRET2_e", "FSR_Et_DRET2_e", 800,800)
c11.SetLogx()
n.Draw("dR:pT>>hFSRFake(20,1,100,20,0,0.5)","pTGen<0&&abs(lID)==11&&"+fsel)
n.Draw("dR:pT>>hFSRTrue(20,1,100,20,0,0.5)","pTGen>0&&abs(lID)==11&&"+fsel,"same")
c11.GetListOfPrimitives()[1].GetYaxis().SetTitle("#Delta R")
c11.GetListOfPrimitives()[1].GetXaxis().SetTitle("E_{T} [GeV]   ")
hFSRFake = c11.GetListOfPrimitives()[2]
hFSRTrue = c11.GetListOfPrimitives()[3]
hFSRFake.SetMarkerStyle(8)
hFSRFake.SetMarkerSize(0.18)
hFSRFake.SetMarkerColor(fakeCol)
hFSRTrue.SetMarkerStyle(8)
hFSRTrue.SetMarkerSize(0.18)
hFSRTrue.SetMarkerColor(trueCol)
hFSRFake.Draw("P")
hFSRTrue.Draw("SAMEP")

f.Draw("SAME")
l.Draw("SAME")
legend2.Draw()


c12 = TCanvas ("FSR_Et_DRET2_m", "FSR_Et_DRET2_m", 800,800)
c12.SetLogx()
n.Draw("dR:pT>>hFSRFake(20,1,100,20,0,0.5)","pTGen<0&&abs(lID)==13&&"+fsel)
n.Draw("dR:pT>>hFSRTrue(20,1,100,20,0,0.5)","pTGen>0&&abs(lID)==13&&"+fsel,"same")
c12.GetListOfPrimitives()[1].GetYaxis().SetTitle("#Delta R")
c12.GetListOfPrimitives()[1].GetXaxis().SetTitle("E_{T} [GeV]   ")
hFSRFake = c12.GetListOfPrimitives()[2]
hFSRTrue = c12.GetListOfPrimitives()[3]
hFSRFake.SetMarkerStyle(8)
hFSRFake.SetMarkerSize(0.18)
hFSRFake.SetMarkerColor(fakeCol)
hFSRTrue.SetMarkerStyle(8)
hFSRTrue.SetMarkerSize(0.18)
hFSRTrue.SetMarkerColor(trueCol)
hFSRFake.Draw("P")
hFSRTrue.Draw("SAMEP")

f.Draw("SAME")
l.Draw("SAME")
legend2.Draw()




### Pattern of fakes caused by FSR attached to leptons with genMother==0 - to be investigated better
c13 = TCanvas ("FSR_Et_DRET2_m_bkg", "Background pattern to be investigated", 800,800)
c13.SetLogx()
n.Draw("dR:pT>>hFSRFake(20,5,100,20,0,0.5)","pTGen<0&&abs(lID)==13&&lTight&&(!SCVetoTight)&&genFS<3 &&genMother==0")
c13.GetListOfPrimitives()[1].GetYaxis().SetTitle("#Delta R")
c13.GetListOfPrimitives()[1].GetXaxis().SetTitle("E_{T} [GeV]   ")
# Try to figure out the pattern we observe
f2 = TF1("f2","0.7/x",5,100)
f2.SetLineWidth(1)
f2.SetLineStyle(3)
f2.Draw("same")

### Purity vs Efficiency as a function of dR/ET2 cut, for different cuts on gRelIso
TotTrue_m =0
TotTrue_e =0
TotTrue =0

cutdRET2 = [0.03,0.025,0.02,0.018,0.016,0.015,0.014,0.013,0.012,0.011,0.010,0.009,0.007,0.006,0.005]
cutIso = [1.6, 1.8, 2.] # 2. is the preselection in our current file
nCuts = len(cutdRET2)
nIsoCuts = len(cutIso)
TrueCounter = [[0. for x in range(nCuts)] for y in range(nIsoCuts)]
FakeCounter = [[0. for x in range(nCuts)] for y in range(nIsoCuts)]

iEntry=0
nEntries = n.GetEntries()

while iEntry< nEntries and n.GetEntry(iEntry):
    iEntry+=1
    if(n.pT>2. and n.lTight and n.genFS<3 and n.genMother==25 and ((abs(n.lID)==13 and abs(n.eta)<2.4) or(abs(n.lID)==11 and abs(n.eta)<2.5))  ):
        if (n.pTGen>0):
            TotTrue+=1
        
        for i in range(nCuts):
            for j in range (nIsoCuts): 
                if (n.gRelIso<cutIso[j] and (not n.maskTight) and n.dRET2 < cutdRET2[i]) : # identical to standard sel
#                if (n.gRelIso<cutIso[j] and n.dRET2 < cutdRET2[i]) : #remove maskTight = treat photons as independent.
                    if n.pTGen>0 :
                        TrueCounter[j][i]+=1
                    else :
                        FakeCounter[j][i]+=1



eff = np.array([0.]*nCuts, dtype = 'double')
pur= np.array([0.]*nCuts, dtype = 'double')

c14 = TCanvas ("FSRPurVsEff","FSRPurVsEff",900,800)

g=[]
for j in range (nIsoCuts):
    print ("--- iso cut ", cutIso[j])
    for i in range (nCuts):
        print ("Eff cut =", cutdRET2[i], "" ,TrueCounter[j][i]/TotTrue)
        eff[i] = TrueCounter[j][i]/TotTrue
        print ("Pur cut =", cutdRET2[i], "",TrueCounter[j][i]/(TrueCounter[j][i]+FakeCounter[j][i]), "\n" )
        pur[i] = TrueCounter[j][i]/(TrueCounter[j][i]+FakeCounter[j][i])
    g.append(TGraph (nCuts,eff,pur))
    g[j].SetMarkerStyle(22+j)
    g[j].SetMarkerColor(j+1)
    if j==0:
        g[j].Draw("APC")
        g[j].GetXaxis().SetTitle("Efficiency")
        g[j].GetYaxis().SetTitle("Purity")
    else:
        g[j].Draw("samePC")

    if j==1 : 
        e1 = np.array([eff[8]], dtype = 'double')
        p1 = np.array([pur[8]], dtype = 'double')
        g1 = TGraph (1,e1,p1)

        g1.SetMarkerStyle(21)
        g1.SetMarkerSize(1.5)
        g1.SetMarkerColor(kMagenta)
        g1.Draw("sameP")

        legend = TLegend()
        legend.AddEntry(g1, "dR/E_{T}^{2}=0.012", "p")
        legend.Draw()
