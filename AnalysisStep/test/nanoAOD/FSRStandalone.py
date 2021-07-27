# Plot files produced with FSRStandaloneAnalyzer
# Extract FSR candidates with:
# grep -e "^FSR:" out_FSR.txt -a | sed -e s/FSR:// -e s/:/" "/g > FSR_tree.txt
# zgrep -ah "^FSR:" AAAOK/ggH125_Chunk*/log.txt.gz | sed -e s/FSR:// -e s/:/" "/g > FSR_tree.txt
from __future__ import print_function
import math
from ROOT import TAttLine,TCanvas, TFile, TColor, gStyle, TH1F, TH1, gROOT, TNtuple
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


# isLoose : gamma passes nanoAOD cuts (preselected->always true)
# isTight : gamma passes H4l cuts
# SCVetoTight : True = would be vetoed by tighter SC veto (as in nanoAOD)
# ltight  : lepton matched to gamma passes tight cuts
# maskLoose : masked by another photon passing loose sel
# maskTight : masked by another photon passing tight sel
# isSelStd  : selected by standard sel: (isTight&&!maskTight)
# isSelNano : tight sel applied after matching loose photons to leptons (emulate tight on top of nanoAOD)
# muEleOverlap : photon can be associated (within the max DR cut) to both a mu and an ele
# pTGen : pT of the gen FSR photon matched to the reco photon. Only FSR photons from W,Z->l are considered.
# genFS: final state (4mu=0, 4e=1, 2e2mu=2)
# genMother: 25 = l from H->Z; 23 = Z, etc.


n=TNtuple("test", "test", "run:ls:event:pT:eta:phi:gRelIso:SCVetoTight:lID:lTight:dR:dRET2:isLoose:isTight:maskLoose:maskTight:isSelStd:isSelNano:muEleOverlap:pTGen:etaGen:phiGen:genFS:genMother")
n.ReadFile("FSR_tree4.txt");

selection = "isSelStd&&lTight&&genFS<3"
#selection = "isSelStd&&lTight&&genFS<3&&genMother==25" # only for matched lepton from H
#selection = "isSelStd&&lTight&&(!SCVetoTight)" # Add tighter SC veto
p_ele = n.Draw("pT",selection+"&&abs(lID)==11&&pTGen>0")/float(n.Draw("pT",selection+"&&abs(lID)==11"))
p_mu  = n.Draw("pT",selection+"&&abs(lID)==13&&pTGen>0")/float(n.Draw("pT",selection+"&&abs(lID)==13"))

print("purity e/mu: ", p_ele, p_mu)

# Effect of tight SC veto
n.Draw("pTGen>0:SCVetoTight>>h(2,0,2,2,0,2)",selection,"textcolz")
h = gROOT.FindObject('h')
fakesRej=h.GetBinContent(2,1)
truesRej=h.GetBinContent(2,2)
allFakes=h.GetBinContent(1,1)+fakesRej
allTrues=h.GetBinContent(1,2)+truesRej
deltaTrues=truesRej/allTrues
deltaPur=fakesRej/(allTrues+allFakes)
print("tight SC veto: Delta_eff = ", deltaTrues, "Delta_fakes", fakesRej/allFakes, "Delta_purity", deltaPur)

# Effect of matching loose photons in nanoAOD+tighter cut at analysis level
# ...

# Effect of matching muons and electrons independently
# ...

# Effect of matching FSR to loose leptons (US's comment). Note: nanoAOD matches to even looser lepts than those considered by maskTight.
# ...
