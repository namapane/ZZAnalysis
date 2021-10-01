
# Plot files produced with FSRAnalysis.py
#
from __future__ import print_function
import math
from ROOT import TAttLine,TCanvas, TFile, TColor, gStyle, TH1F, TH1, gROOT
from ROOT import kBlack, kBlue, kRed, kOrange, kGreen, kPink,kCyan

def printCanvases(type="png", path=".") :
    canvases = gROOT.GetListOfCanvases()
    for c in canvases :
        c.Print(path+"/"+c.GetTitle()+"."+type)

ifile=TFile("fsrHisto.root")

gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetHistLineWidth(4)
gStyle.SetPadGridX(True);
gStyle.SetPadGridY(True);
gStyle.SetGridColor(0);
gStyle.SetGridStyle(3);
gStyle.SetGridWidth(1);


TH1.SetDefaultSumw2()

def plotLepPt(hMu,hEle,hMuFake,hEleFake):
    hMu.SetLineColor(kRed)
    hMu.GetXaxis().SetTitle("p_{t} [GeV/c^{2}]")
    hMu.GetYaxis().SetTitle("Arbitrary Units")
    hEle.SetLineColor(kBlue)
    hMuFake.SetLineColor(kPink-9)
    hEleFake.SetLineColor(kCyan-6)
    hMu.Draw()
    hEle.Draw("same")
    hMuFake.Draw("same")
    hEleFake.Draw("same")

c1 = TCanvas( 'FSRLepPts', 'FSRLepPts',800,800)
plotLepPt(ifile.fsrMuPt, ifile.fsrElePt, ifile.fsrMuPtFake, ifile.fsrElePtFake)


def plot(hFSR,hNoFSR,hTrue,hFalse) :
    hFSR.SetLineColor(kBlack)
    hFSR.GetXaxis().SetTitle("m_{4 #it{l} } [GeV/c^{2}]")
    hFSR.GetYaxis().SetTitle("Arbitrary Units")    
    hNoFSR.SetLineColor(kOrange)
    hTrue.SetLineColor(kRed)
    hFalse.SetLineColor(kGreen)
    hFSR.SetLineWidth(2)
    hNoFSR.SetLineWidth(2)
    hTrue.SetLineWidth(2)
    hFalse.SetLineWidth(2)
    hFSR.Draw()
    hNoFSR.Draw("same")
    hTrue.Draw("same")
    hFalse.Draw("same")
    hFSR.Draw("same") #redraw on top
    
c2 = TCanvas( 'CandMass_4e_FSR', 'CandMass_4e_FSR',800,800)
plot(ifile.ZZMassFSR_4e, ifile.ZZMassNoFSR4e, ifile.ZZMassFSRTrue4e, ifile.ZZMassFSRFake4e)

c3 = TCanvas( 'CandMass_4mu_FSR', 'CandMass_4mu_FSR',800,800)
plot(ifile.ZZMassFSR_4mu, ifile.ZZMassNoFSR4mu, ifile.ZZMassFSRTrue4mu, ifile.ZZMassFSRFake4mu)

c4 = TCanvas("CandMass_2e2mu_FSR","CandMass_2e2mu_FSR",800,800)
plot(ifile.ZZMassFSR_2e2mu,ifile.ZZMassNoFSR2e2mu,ifile.ZZMassFSRTrue2e2mu,ifile.ZZMassFSRFake2e2mu)



purity = 1.-ifile.fsrET2vsDR_fake.GetEntries()/(ifile.fsrET2vsDR_fake.GetEntries()+ifile.fsrET2vsDR_true.GetEntries())
print("FSR purity: = {:.3g}".format(purity))


h_purityMu = TH1F("FSR purity muon","FSR purity muon",40,0.,200.)
c5=TCanvas ("FSRPurity_muon","FSRPurity_muon",800,800)
h_purityMu.Divide(ifile.fsrMuPtTrue,ifile.fsrMuPt,1.,1.,"b")
h_purityMu.SetLineWidth(1)
h_purityMu.Draw("PE")

h_purityEle = TH1F("FSR purity ele","FSR purity ele",40,0.,200.)
c6=TCanvas ("FSRPurity_ele","FSRPurity_ele",800,800)
h_purityEle.Divide(ifile.fsrElePtTrue,ifile.fsrElePt,1.,1.,"b")
h_purityEle.SetLineWidth(1)
h_purityEle.Draw("PE")

h_purMufsrEta = TH1F("FSR Eta purity muon","FSR Eta purity muon",25,0.,2.5)
c7=TCanvas ("FSR Eta purity muon","FSR Eta purity muon",800,800)
h_purMufsrEta.Divide(ifile.MufsrEtaTrue,ifile.MufsrEta,1.,1.,"b")
h_purMufsrEta.SetLineWidth(1)
h_purMufsrEta.Draw("PE")
h_purMufsrEta.GetXaxis().SetTitle("#eta")
