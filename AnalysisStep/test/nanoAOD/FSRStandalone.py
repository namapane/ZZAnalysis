# Plot files produced with FSRStandaloneAnalyzer
# Extract FSR candidates with:
# grep -e "^FSR:" out_FSR.txt -a | sed -e s/FSR:// -e s/:/" "/g > FSR_tree.txt
#
from __future__ import print_function
import math
from ROOT import TAttLine,TCanvas, TFile, TColor, gStyle, TH1F, TH1, gROOT, TNtuple
from ROOT import kBlack, kBlue, kRed, kOrange, kGreen, kPink,kCyan

def printCanvases(type="png", path=".") :
    canvases = gROOT.GetListOfCanvases()
    for c in canvases :
        c.Print(path+"/"+c.GetTitle()+"."+type)

gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetHistLineWidth(4)
gStyle.SetPadGridX(True);
gStyle.SetPadGridY(True);
gStyle.SetGridColor(0);
gStyle.SetGridStyle(3);
gStyle.SetGridWidth(1);

TH1.SetDefaultSumw2()


n=TNtuple("test", "test", "run:ls:event:pT:eta:phi:gRelIso:SCVetoTight:lID:lTight:dR:dRET2:isLoose:isTight:maskLoose:maskTight:isSelStd:isSelNano:muEleOverlap:pTGen:etaGen:phiGen")
n.ReadFile("FSR_tree2.txt");


p_ele = n.Draw("pT","isTight&&abs(lID)==11&&pTGen>0")/float(n.Draw("pT","isTight&&abs(lID)==11"))
p_mu  = n.Draw("pT","isTight&&abs(lID)==13&&pTGen>0")/float(n.Draw("pT","isTight&&abs(lID)==13"))

print("purity e/mu: ", p_ele, p_mu)

n.Draw("pTGen>0:SCVetoTight>>h(2,0,2,2,0,2)","isTight","textcolz")
