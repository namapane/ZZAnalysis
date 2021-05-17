# Plot files produced with FSRtest.py
#
from __future__ import print_function
import math
from ROOT import TCanvas, TFile, TColor, gStyle
from ROOT import kBlack, kBlue, kRed, kOrange, kGreen

ifile=TFile("fsrHisto.root")

gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetHistLineWidth(4)

c1 = TCanvas( 'FSRLepPts', 'FSRLepPts',)
ifile.fsrMuPt.SetLineColor(kRed)
ifile.fsrElePt.SetLineColor(kBlue)
ifile.fsrMuPt.Draw()
ifile.fsrElePt.Draw("same")

c2 = TCanvas( 'FSRMass_4e', '4e Candidates with FSR',)
ifile.ZZMassFSR_4e.SetLineColor(kBlack)
ifile.ZZMassNoFSR4e.SetLineColor(kOrange)
ifile.ZZMassFSRTrue4e.SetLineColor(kRed)
ifile.ZZMassFSRFake4e.SetLineColor(kGreen)
ifile.ZZMassFSR_4e.SetLineWidth(2)
ifile.ZZMassNoFSR4e.SetLineWidth(2)
ifile.ZZMassFSRTrue4e.SetLineWidth(2)
ifile.ZZMassFSRFake4e.SetLineWidth(2)
ifile.ZZMassFSR_4e.Draw()
ifile.ZZMassNoFSR4e.Draw("same")
ifile.ZZMassFSRTrue4e.Draw("same")
ifile.ZZMassFSRFake4e.Draw("same")
ifile.ZZMassFSR_4e.Draw("same") #redraw on top

c3 = TCanvas( 'FSRMass_4mu', '4mu Candidates with FSR',)
ifile.ZZMassFSR_4mu.SetLineColor(kBlack)
ifile.ZZMassNoFSR4mu.SetLineColor(kOrange)
ifile.ZZMassFSRTrue4mu.SetLineColor(kRed)
ifile.ZZMassFSRFake4mu.SetLineColor(kGreen)
ifile.ZZMassFSR_4mu.SetLineWidth(2)
ifile.ZZMassNoFSR4mu.SetLineWidth(2)
ifile.ZZMassFSRTrue4mu.SetLineWidth(2)
ifile.ZZMassFSRFake4mu.SetLineWidth(2)
ifile.ZZMassFSR_4mu.Draw()
ifile.ZZMassNoFSR4mu.Draw("same")
ifile.ZZMassFSRTrue4mu.Draw("same")
ifile.ZZMassFSRFake4mu.Draw("same")
ifile.ZZMassFSR_4mu.Draw("same") #redraw on top


def plotTrueFake(hFSR,hNoFSR,hTrue,hFalse) :
    hFSR.SetLineColor(kBlack)
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
    
c4 = TCanvas("FSRMass_2e2mu","FSRMass_2e2mu")
plotTrueFake(ifile.ZZMassFSR_2e2mu,ifile.ZZMassNoFSR2e2mu,ifile.ZZMassFSRTrue2e2mu,ifile.ZZMassFSRFake2e2mu)


purity = 1.-ifile.fsrET2vsDR_fake.GetEntries()/(ifile.fsrET2vsDR_fake.GetEntries()+ifile.fsrET2vsDR_true.GetEntries())
print("FSR purity: = {:.3g}".format(purity))
