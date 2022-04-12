from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import *
import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

class weightFiller(Module):
    def __init__(self, XS, APPLY_K_NNLOQCD_ZZGG, APPLY_K_NNLOQCD_ZZQQB, APPLY_K_NNLOEW_ZZQQB, APPLY_QCD_GGF_UNCERT):
        self.writeHistFile = False
        self.XS = XS
        self.APPLY_K_NNLOQCD_ZZGG = APPLY_K_NNLOQCD_ZZGG
        self.APPLY_K_NNLOQCD_ZZQQB = APPLY_K_NNLOQCD_ZZQQB
        self.APPLY_K_NNLOEW_ZZQQB = APPLY_K_NNLOEW_ZZQQB
        self.APPLY_QCD_GGF_UNCERT = APPLY_QCD_GGF_UNCERT

        basePath='%s/src/ZZAnalysis/AnalysisStep/' % os.environ['CMSSW_BASE']

        ## ggZZ QCD k-factors
        self.spkfactor_ggzz_nnlo = [None]*9
        self.spkfactor_ggzz_nlo  = [None]*9
        if self.APPLY_K_NNLOQCD_ZZGG != 0 :            
            strZZGGKFVar = ["Nominal", "PDFScaleDn", "PDFScaleUp", "QCDScaleDn", "QCDScaleUp", "AsDn", "AsUp", "PDFReplicaDn", "PDFReplicaUp"]
            # NNLO
            ggZZKFactorFile = ROOT.TFile.Open(basePath+'data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root')
            for i in range(0,9) :
                 for i in range(0,9) :
                     self.spkfactor_ggzz_nnlo[i] = (ggZZKFactorFile.Get('sp_kfactor_%s' % strZZGGKFVar[i])).Clone('sp_kfactor_%s_NNLO' % strZZGGKFVar[i])
            ggZZKFactorFile.Close()
            # NLO
            ggZZKFactorFile = ROOT.TFile.Open(basePath+'data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NLO_NNPDF_NarrowWidth_13TeV.root')
            for i in range(0,9) :
                 for i in range(0,9) :
                     self.spkfactor_ggzz_nlo[i] = (ggZZKFactorFile.Get('sp_kfactor_%s' % strZZGGKFVar[i])).Clone('sp_kfactor_%s_NNLO' % strZZGGKFVar[i])
            ggZZKFactorFile.Close()            

        # qqZZ QCD K-fatcors (for qqZZ)
        # We include these using plain root, but the code could be called directly from the .so lib
        if self.APPLY_K_NNLOQCD_ZZQQB :
            if "/kFactors_C.so" not in ROOT.gSystem.GetLibraries() :
                p_helper = basePath+'src/kFactors.C+'
                print('Loading C++ helper from ' + p_helper)
                ROOT.gROOT.ProcessLine('.L ' + p_helper)
            self.kfactor_qqZZ_qcd_M = ROOT.kfactor_qqZZ_qcd_M

        if self.APPLY_K_NNLOEW_ZZQQB :
            #ewkTable = EwkCorrections::readFile_and_loadEwkTable(basePath+'data/kfactors/ZZ_EwkCorrections.dat')
            pass

        # ggH NNLOPS weights
        if self.APPLY_QCD_GGF_UNCERT :        
            NNLOPS_weight_file = ROOT.TFile.Open(basePath+'data/ggH_NNLOPS_Weights/NNLOPS_reweight.root')
            self.gr_NNLOPSratio_pt_powheg_0jet = NNLOPS_weight_file.Get("gr_NNLOPSratio_pt_powheg_0jet")
            self.gr_NNLOPSratio_pt_powheg_1jet = NNLOPS_weight_file.Get("gr_NNLOPSratio_pt_powheg_1jet")
            self.gr_NNLOPSratio_pt_powheg_2jet = NNLOPS_weight_file.Get("gr_NNLOPSratio_pt_powheg_2jet")
            self.gr_NNLOPSratio_pt_powheg_3jet = NNLOPS_weight_file.Get("gr_NNLOPSratio_pt_powheg_3jet")


    ### Imported from: https://github.com/CJLST/ZZAnalysis/blob/Run2_CutBased_UL/AnalysisStep/test/Ntuplizers/HZZ4lNtupleMaker.cc#L1950
    def evalSpline(self, sp, xval):
        xmin = sp.GetXmin()
        xmax = sp.GetXmax()
        res = 0.
        if xval<xmin :
            res=sp.Eval(xmin)
            deriv=sp.Derivative(xmin)
            res += deriv*(xval-xmin)
        elif (xval>xmax) :
            res=sp.Eval(xmax)
            deriv=sp.Derivative(xmax)
            res += deriv*(xval-xmax)
        else :
            res=sp.Eval(xval)
        return res
        

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.APPLY_K_NNLOQCD_ZZQQB :
           self.out.branch("KFactor_QCD_qqZZ_M_Weight", "F")            
        if self.APPLY_K_NNLOQCD_ZZGG > 0 : 
            self.out.branch("KFactor_QCD_ggZZ_Nominal_Weight", "F")
        if self.APPLY_K_NNLOEW_ZZQQB :
            self.out.branch("KFactor_EW_qqZZ_Weight", "F")
        if self.APPLY_QCD_GGF_UNCERT :
            self.out.branch("ggH_NNLOPS_Weight", "F")

            
        self.out.branch("overallEventWeight", "F") #Overall weight, including relevant k-factors and corrections


    def analyze(self, event):
        KFactor_EW_qqZZ =1.        
        KFactor_QCD_qqZZ_M = 1.
        KFactor_QCD_ggZZ_Nominal = 1.
        ggH_NNLOPS_Weight = 1.

        ### QCD weights for ggH, ggZZ
        if self.APPLY_K_NNLOQCD_ZZGG == 1 or self.APPLY_K_NNLOQCD_ZZGG == 2 :
            KFactor_QCD_ggZZ_Nominal = self.evalSpline(self.spkfactor_ggzz_nnlo[0], event.GenZZMass) #1: NNLO/LO
            if self.APPLY_K_NNLOQCD_ZZGG == 2 : # 2:NNLO/NLO
                denominator = self.evalSpline(self.spkfactor_ggzz_nlo[0], event.GenZZMass)
                KFactor_QCD_ggZZ_Nominal /= denominator
        elif self.APPLY_K_NNLOQCD_ZZGG == 3 :# 2:NLO/LO
            KFactor_QCD_ggZZ_Nominal = self.evalSpline(self.spkfactor_ggzz_nlo[0], event.GenZZMass)
        elif self.APPLY_K_NNLOQCD_ZZGG !=0 :
            print ("Unsupported: APPLY_K_NNLOQCD_ZZGG=", self.APPLY_K_NNLOQCD_ZZGG) 
            exit(1)

        ### QCD and EW weights for qqZZ
        if self.APPLY_K_NNLOQCD_ZZQQB :
            flavor = 1 # same Z flavors
            if event.GenZZFinalState == 14641 or event.GenZZFinalState == 28561 or event.GenZZFinalState == 50625 :
                flavor = 2 # differenr Z flavors
            KFactor_QCD_qqZZ_M = self.kfactor_qqZZ_qcd_M(event.GenZZMass, flavor, 2)/self.kfactor_qqZZ_qcd_M(event.GenZZMass, flavor, 1)

        if self.APPLY_K_NNLOEW_ZZQQB :
            KFactor_EW_qqZZ = 1 #FIXME: this is 1 up to 2*mZ. To be implemented for higher masses.

        if self.APPLY_QCD_GGF_UNCERT :
            htxsNJets = event.HTXS_njets30
            htxsHPt = event.HTXS_Higgs_pt
            if htxsNJets==0 :
                ggH_NNLOPS_Weight = self.gr_NNLOPSratio_pt_powheg_0jet.Eval(min(htxsHPt, 125.0))
            elif htxsNJets==1 :
                ggH_NNLOPS_Weight = self.gr_NNLOPSratio_pt_powheg_1jet.Eval(min(htxsHPt,625.0))
            elif htxsNJets==2 :
                ggH_NNLOPS_Weight = self.gr_NNLOPSratio_pt_powheg_2jet.Eval(min(htxsHPt,800.0))
            elif htxsNJets>=3 :
                ggH_NNLOPS_Weight = self.gr_NNLOPSratio_pt_powheg_3jet.Eval(min(htxsHPt,925.0))


        # L1 pre-firing weights. FIXME: not available in Nano02Apr2020
        # L1prefiringWeight = event.L1PreFiringWeight_Nom
        # L1prefiringWeightUp = event.L1PreFiringWeight_Up
        # L1prefiringWeightDn = event.L1PreFiringWeight_Dn
        
        #FIXME: event.ZZ_dataMCWeight is not included, since that can be stored per-candidate if storeAllCands=True.
        w_total = self.XS * event.Generator_weight * event.puWeight * KFactor_EW_qqZZ * KFactor_QCD_ggZZ_Nominal * KFactor_QCD_qqZZ_M * ggH_NNLOPS_Weight

        if self.APPLY_K_NNLOQCD_ZZQQB :
            self.out.fillBranch("KFactor_QCD_qqZZ_M_Weight", KFactor_QCD_qqZZ_M)
        if self.APPLY_K_NNLOQCD_ZZGG > 0 : 
            self.out.fillBranch("KFactor_QCD_ggZZ_Nominal_Weight", KFactor_QCD_ggZZ_Nominal)
        if self.APPLY_K_NNLOEW_ZZQQB :
            self.out.fillBranch("KFactor_EW_qqZZ_Weight", KFactor_EW_qqZZ)
        if self.APPLY_QCD_GGF_UNCERT :
            self.out.fillBranch("ggH_NNLOPS_Weight", ggH_NNLOPS_Weight)

        self.out.fillBranch("overallEventWeight", w_total)

        return True


