from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR


def passEleBDT(pt, SCeta, BDT, era = "2018") :
    fSCeta = abs(SCeta)  
    # FIXME: This should be era-dependent.
    # Using 2017 WP for ElectronMVAEstimatorRun2Fall17IsoV2Values since other BDTs are not available in nanoAOD
    return (pt<=10. and ((fSCeta<0.8 and BDT > 0.85216885148) or (fSCeta>=0.8 and fSCeta<1.479 and BDT >  0.82684550976) or (fSCeta>=1.479 and BDT >  0.86937630022))) or (pt>10. and  ((fSCeta<0.8 and BDT >  0.98248928759) or (fSCeta>=0.8 and fSCeta<1.479 and BDT >  0.96919224579) or (fSCeta>=1.479 and BDT >  0.79349796445)))


class lepFiller(Module):
    def __init__(self, cuts):
        self.writeHistFile=False
        self.cuts = cuts
        self.muLooseId = cuts["muLooseId"]
        self.muTightId = cuts["muTightId"]
        self.eleLooseId = cuts["eleLooseId"]
        self.eleTightId = cuts["eleTightId"]
        

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree        
        
        self.out.branch("FsrPhoton_mass", "F", lenVar="nFsrPhoton") # Hack so that photon.p4() works

        self.out.branch("Electron_isLoose", "O", lenVar="nElectron")
        self.out.branch("Electron_isTightIso", "O", lenVar="nElectron")
        self.out.branch("Electron_myFsrPhotonIdx", "I", lenVar="nElectron")
        self.out.branch("Electron_pfRelIso03FsrCorr", "F", lenVar="nElectron")

        self.out.branch("Muon_isLoose", "O", lenVar="nMuon")
        self.out.branch("Muon_isTightIso", "O", lenVar="nMuon")
        self.out.branch("Muon_myFsrPhotonIdx", "I", lenVar="nMuon")
        self.out.branch("Muon_pfRelIso03FsrCorr", "F", lenVar="nMuon")

#    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#        pass

    # Recompute lepton isolation subtracting energy of the given FSR photons
    def isoFsrCorr(self, l, selectedFSR) : 
        combRelIsoPFFSRCorr = l.pfRelIso03_all
        for f in selectedFSR :
            dR = deltaR(l.eta, l.phi, f.eta, f.phi)
            if dR >0.01 and dR < 0.3 :
                combRelIsoPFFSRCorr = max(0., l.pfRelIso03_all-f.pt/l.pt)
        return combRelIsoPFFSRCorr


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")

        # IDs (no iso)
        eleLoose = list(self.eleLooseId(e) for e in electrons)
        muLoose = list(self.muLooseId(m) for m in muons)
        eleTight = list(self.eleTightId(e) for e in electrons)
        muTight = list(self.muTightId(m) for m in muons)  

        # Skip events that do not contain enough tight leptons (note: for CRs, this should be modified)
        if (len(eleTight)+len(muTight)<4) : return False

        
        # Re-match FSR with Loose leptons, as default associati1on may be wrong due to different cuts, masking/stealing.
        # Store new FSR idxs, dRET2, and commpute FSR-corrected isolation
        fsrPhoton_myMuonIdx=[-1]*len(fsrPhotons)
        fsrPhoton_myElectronIdx=[-1]*len(fsrPhotons)
        fsrPhoton_mydROverEt2=[-1.]*len(fsrPhotons)
        muFsrPhotonIdx=[-1]*len(muons)
        eleFsrPhotonIdx=[-1]*len(electrons)

        muPhotondREt2 = [1e9]*len(muons)
        elePhotondREt2 = [1e9]*len(electrons)
        
        for ifsr, fsr in enumerate (fsrPhotons):
            if (fsr.pt < 2. or abs(fsr.eta) > 2.4 or fsr.relIso03 > self.cuts["fsr_Iso"] ) : continue
            dRmin = 0.5 # max dR for association
            closestMu=-1
            closestEle=-1
            for ilep, lep in enumerate(muons):
                if muLoose[ilep] :
                    dR =  deltaR(lep.eta, lep.phi, fsr.eta, fsr.phi)
                    if dR < dRmin and dR > 0.001 and dR/fsr.pt/fsr.pt < self.cuts["fsr_dRET2"]:
                        dRmin = dR
                        closestMu=ilep

            for ilep, lep in enumerate(electrons):
                if eleLoose[ilep] :
                    dR =  deltaR(lep.eta, lep.phi, fsr.eta, fsr.phi)
                    if dR < dRmin and dR > 0.001 and dR/fsr.pt/fsr.pt < self.cuts["fsr_dRET2"]:
                        dRmin = dR
                        closestMu=-1
                        closestEle=ilep

            if (closestMu>=0 or closestEle>=0) :
                dREt2 = dRmin/fsr.pt/fsr.pt
                fsrPhoton_mydROverEt2[ifsr] = dREt2
                if (closestMu>=0) :
                    fsrPhoton_myMuonIdx[ifsr] = closestMu
                    if dREt2<muPhotondREt2[closestMu] : #keep only the lowest-dRET2 for each lepton
                        muPhotondREt2[closestMu] = dREt2
                        muFsrPhotonIdx[closestMu] = ifsr
                if (closestEle>=0) :
                    fsrPhoton_myElectronIdx[ifsr] = closestEle
                    if dREt2<elePhotondREt2[closestEle] :
                        elePhotondREt2[closestEle] = dREt2
                        eleFsrPhotonIdx[closestEle] = ifsr


        # Recompute isolation removing all selected FSRs
        ele_isoFsrCorr = [-1.]*len(electrons)
        mu_isoFsrCorr = [-1.]*len(muons)
        ele_tightIso = [False]*len(electrons)
        mu_tightIso = [False]*len(muons)

        selectedFSR = []
        for ifsr in muFsrPhotonIdx + eleFsrPhotonIdx :
            if ifsr>=0 : selectedFSR.append(fsrPhotons[ifsr])

        for ilep,lep in enumerate(muons) :
            mu_isoFsrCorr[ilep] = self.isoFsrCorr(lep, selectedFSR)
            mu_tightIso[ilep] = muTight[ilep] and mu_isoFsrCorr[ilep] < self.cuts["relIso"]
            
        for ilep,lep in enumerate(electrons) :
            ele_isoFsrCorr[ilep] = self.isoFsrCorr(lep, selectedFSR) 
            ele_tightIso[ilep] = eleTight[ilep]# no iso cut for electrons
 
        # Tight ID, including FSR-corrected iso

        fsrM = [0.]*len(fsrPhotons)
        self.out.fillBranch("FsrPhoton_mass", fsrM)

        self.out.fillBranch("Electron_isLoose", eleLoose)
        self.out.fillBranch("Electron_isTightIso", ele_tightIso)
        self.out.fillBranch("Electron_myFsrPhotonIdx", eleFsrPhotonIdx)
        self.out.fillBranch("Electron_pfRelIso03FsrCorr", ele_isoFsrCorr)

        self.out.fillBranch("Muon_isLoose", muLoose)
        self.out.fillBranch("Muon_isTightIso", mu_tightIso)
        self.out.fillBranch("Muon_myFsrPhotonIdx", muFsrPhotonIdx)
        self.out.fillBranch("Muon_pfRelIso03FsrCorr", mu_isoFsrCorr)

        return True
