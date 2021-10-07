/** \class PhotonFiller
 *
 *  Monolithic clone of lepton-level FSR selection workflow, for performance studies
 *
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include <ZZAnalysis/AnalysisStep/interface/PhotonFwd.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>

#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
#include <Math/VectorUtil.h>
#include <TMath.h>

#include <vector>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;

struct FSRCandidate {
  edm::Ptr<pat::PackedCandidate> g;
  double gRelIso;
  bool SCVetoTight;
  float vetoer_pT;
  double dRMinEle; // distance to closest ele
  double dRMinMu;  // distance to closest mu
  const pat::Electron* closestEle = nullptr;
  const pat::Muon * closestMu = nullptr;
  const reco::Candidate* closestLep = nullptr; // closest of the above
  int igen; // index of matched gen FSR photon (if any)
  
  // Parameters, w.r.t. closest lepton
  float dR() {return min(dRMinMu,dRMinEle);}
  float dRET2() {return dR()/g->pt()/g->pt();}
  bool isFSRLoose() {return gRelIso<2   && dRET2()<0.05;}  // nanoAOD selection
  bool isFSRTight() {return gRelIso<1.8 && dRET2()<0.012;} // H4l selection
  bool isLepTight() {
    if (dRMinEle<dRMinMu) return closestEle->userFloat("isGood");
    else return closestMu->userFloat("isGood");
  }
  const reco::GenParticle* genLepton() {
    if (dRMinEle<dRMinMu) return closestEle->genLepton();
    else return closestMu->genLepton();
  }  

  // w.r.t. closest mu
  bool isFSRLoose_mu() {return closestMu!=nullptr && gRelIso<2   && dRMinMu/g->pt()/g->pt()<0.05;}
  bool isFSRTight_mu() {return closestMu!=nullptr && gRelIso<1.8 && dRMinMu/g->pt()/g->pt()<0.012;}
  
  // w.r.t. closest ele
  bool isFSRLoose_ele() {return closestEle!=nullptr && gRelIso<2   && dRMinEle/g->pt()/g->pt()<0.05;}
  bool isFSRTight_ele() {return closestEle!=nullptr && gRelIso<1.8 && dRMinEle/g->pt()/g->pt()<0.012;}

};

  


class FSRStandaloneAnalyzer : public edm::EDAnalyzer {
 public:
  /// Constructor
  explicit FSRStandaloneAnalyzer(const edm::ParameterSet&);
    
  /// Destructor
  ~FSRStandaloneAnalyzer(){};  

 private:
  virtual void beginJob(){};  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > pfCandToken;
  edm::EDGetTokenT<pat::MuonCollection> muonToken;
  edm::EDGetTokenT<pat::MuonRefVector> muonAllToken;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken;
  edm::EDGetTokenT<pat::ElectronCollection> electronAllToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > genParticleToken;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
//   edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PupInfoToken;
//   edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken;

  bool debug;
};


FSRStandaloneAnalyzer::FSRStandaloneAnalyzer(const edm::ParameterSet& iConfig) :
  pfCandToken(consumes<edm::View<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates"))), // input for photons
  muonToken(consumes<pat::MuonCollection>(edm::InputTag("softMuons"))), 
  muonAllToken(consumes<pat::MuonRefVector>(edm::InputTag("finalMuons"))), // selected as in nanoAOD
  electronToken(consumes<pat::ElectronCollection>(edm::InputTag("cleanSoftElectrons"))),
  electronAllToken(consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"))), // all
  genParticleToken(consumes<edm::View<reco::Candidate> >( edm::InputTag("prunedGenParticles"))),
  genInfoToken(consumes<GenEventInfoProduct>( edm::InputTag("generator"))),
  debug(iConfig.getUntrackedParameter<bool>("debug",false)) {}


void
FSRStandaloneAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //--- PF cands, to look for photon candidates
  edm::Handle<edm::View<pat::PackedCandidate> > pfCands; 
  iEvent.getByToken(pfCandToken, pfCands);

  edm::Handle<pat::MuonCollection> muonHandle;
  iEvent.getByToken(muonToken, muonHandle);

  edm::Handle<pat::ElectronCollection> electronHandle;
  iEvent.getByToken(electronToken, electronHandle);

  edm::Handle<pat::MuonRefVector> muonAllHandle;
  iEvent.getByToken(muonAllToken, muonAllHandle);

  edm::Handle<pat::ElectronCollection> electronAllHandle;
  iEvent.getByToken(electronAllToken, electronAllHandle);

  // MC Truth
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByToken(genParticleToken, genParticles);
  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(genInfoToken, genInfo);
  MCHistoryTools mch(iEvent, "", genParticles, genInfo);

  int genFinalState = mch.genFinalState();

  vector<const reco::Candidate *> genFSR = mch.genFSR();

  vector<int> nRecoFSRMatchedToGen(genFSR.size(),0); // check: how many photons are matched to each gen FSR?

  //----------------------
  // Loop on photons
  //----------------------

  vector<FSRCandidate> photons;
  
  typedef map<const reco::Candidate*, vector<size_t> > PhIdxsByLep;
  PhIdxsByLep photonsByLep; // indexes of photon candidates associated to each lep


  for (unsigned int i=0;i<pfCands->size();++i) {
      
    // Get the candidate as edm::Ptr
    edm::Ptr<pat::PackedCandidate> g = pfCands->ptrAt(i);

    // We only want photons
    if (g->pdgId()!=22) continue;

    // Photon preselection (is currently already applied on pat::PackedCandidate collection?)
    if (!(g->pt()>2. && fabs(g->eta())<2.4)) continue;

    //---------------------
    // // Supercluster veto
    //---------------------
    bool SCVeto=false;      // H4l veto (veto only against loose ID+SIP)
    bool SCVetoTight=false; // nanoAOD-style (veto any slimmedElectrons, without even a pT cut)


    // H4l-style veto: vs loose ID + SIP leptons 
    if (electronHandle->size()>0) {
      for (unsigned int j = 0; j< electronHandle->size(); ++j){
	const pat::Electron* e = &((*electronHandle)[j]);
	if ((e->associatedPackedPFCandidates()).size()) {
	  edm::RefVector < pat::PackedCandidateCollection > pfcands = e->associatedPackedPFCandidates();
	  for ( auto itr: pfcands ) {
	    if (g.get()==&(*itr)) {
	      if (e->userFloat("isSIP")) SCVeto=true;
// 	      if (debug) cout << "SC veto (loose only=" << SCVeto << ") " << itr->eta() << " " << itr->phi() << " " 
// 			      << fabs(g->eta() - itr->eta()) << " " << reco::deltaPhi(g->phi(), itr->phi()) << endl;
	    }
	  }
	}
      }
    }
    // nano-AOD style: vs any slimmedElectrons, without even a pT cut
    float vetoer_pT = -1;
    if (electronAllHandle->size()>0) {
      for (unsigned int j = 0; j< electronAllHandle->size(); ++j){
	const pat::Electron* e = &((*electronAllHandle)[j]);
	if ((e->associatedPackedPFCandidates()).size()) {
	  edm::RefVector < pat::PackedCandidateCollection > pfcands = e->associatedPackedPFCandidates();
	  for ( auto itr: pfcands ) {
	    if (g.get()==&(*itr)) {
	      SCVetoTight=true;
	      vetoer_pT = e->pt();
	    }
	  }
	}
      }
    }
    if (debug) cout << "FSRStandaloneAnalyzer: gamma:" << g->pt() << " " << g->eta() << " " << g->phi() << " SCVeto L/T: " << SCVeto << "/" << SCVetoTight << endl;

    if (SCVeto) continue;
      
    reco::PFCandidate gamma(0, g->p4(), reco::PFCandidate::gamma);
    double neu, chg, chgByWorstPV;
    LeptonIsoHelper::fsrIso(&gamma, pfCands, neu, chg, chgByWorstPV);
    double gRelIso = (neu + chg)/gamma.pt();    

    double dRGenVsReco = -1.;
    int igen = MCHistoryTools::fsrMatch(&gamma,genFSR);
    if (igen>=0) {
      dRGenVsReco = ROOT::Math::VectorUtil::DeltaR(genFSR[igen]->momentum(),gamma.momentum());
      if (dRGenVsReco<0.3) { //Matching cut (DR<0.3 and pt_gen>2 already applied in fsrMatch; add other criteria?)
	nRecoFSRMatchedToGen[igen]++; // FIXME
      } else {
	igen=-1;
      }
    }

    //------------------------------------------------------
    // Get the closest lepton among those satisfying loose ID + SIP
    //------------------------------------------------------
    double dRMinEle(10e9);
    double dRMinMu(10e9);
    const pat::Electron* closestEle = 0;
    const pat::Muon* closestMu = 0;
    
    // Loop over pat::Muon
    for (unsigned int j = 0; j< muonHandle->size(); ++j){
      const pat::Muon* m = &((*muonHandle)[j]);
      if (! m->userFloat("isSIP")) continue;
      double dR = ROOT::Math::VectorUtil::DeltaR(m->momentum(),g->momentum());
      if(debug) cout << "muon pt = " << m->pt() << " photon pt = " << g->pt() << " dR = " << dR 
		     << " " << m->simType()
		     << " " << m->simExtType()
		     << " " << m->simFlavour()
		     << " " << m->simHeaviestMotherFlavour()
		     << endl;
      if (dR>=0.5) continue;
      if (dR<dRMinMu) {
	dRMinMu = dR;
	closestMu = m;
      }
    }//end loop over muon collection
    
    //---------------------
    // Loop over pat::Electron
    //---------------------
    for (unsigned int j = 0; j< electronHandle->size(); ++j){
      const pat::Electron* e = &((*electronHandle)[j]);
      if (! e->userFloat("isSIP")) continue;
      double dR = ROOT::Math::VectorUtil::DeltaR(e->momentum(),g->momentum());
      
      if(debug) cout << "ele pt = " << e->pt() << " photon pt = " << g->pt() << " dR = " << dR
		     << " " << (e->genLepton()==nullptr?0:e->genLepton()->pdgId())
		     << endl;
      if (dR>=0.5) continue;
      if (dR<dRMinEle) {
	dRMinEle = dR;
	closestEle = e;	
      }
    } //end loop over electron collection


    if(closestEle!=nullptr || closestMu!=nullptr) {
      size_t iPhoton = photons.size();

      const reco::Candidate* closestLep = closestMu;
      if (dRMinEle<dRMinMu) {
	closestLep = closestEle;
      }
      photonsByLep[closestLep].push_back(iPhoton); // store photon candidates by lepton
      
      FSRCandidate gc;
      gc.g=g;
      gc.gRelIso=gRelIso;
      gc.SCVetoTight=SCVetoTight;
      gc.vetoer_pT=vetoer_pT;
      gc.dRMinEle=dRMinEle;
      gc.dRMinMu=dRMinMu;
      gc.closestLep=closestLep;
      gc.closestMu=closestMu;
      gc.closestEle=closestEle;
      gc.igen = igen;

      photons.push_back(gc);
    }
  } // end of loop over photon collection


  // Loop on all photon candidates to print info
  for (size_t iPhoton=0; iPhoton<photons.size(); ++iPhoton) {	
    FSRCandidate& gc = photons[iPhoton];
    if (!gc.isFSRLoose()) continue;

    int lID = gc.closestLep->pdgId(); // ID of the closest lepton
    float lpT = gc.closestLep->pt();

    // Are there other photons associated to the same lepton?
    vector<size_t> allPhotons = photonsByLep[gc.closestLep];
    bool maskLoose = false; // this photon is masked by another photon passing loose cuts
    bool maskTight = false; // this photon is masked by another  photon passing tight cuts
    for (const size_t& iOther : allPhotons){
      FSRCandidate& gOther = photons[iOther];
      if (iOther!=iPhoton && gOther.dRET2()<gc.dRET2()){ // Another photon would be selected for this lepton
	if (gOther.isFSRLoose()) maskLoose = true;
	if (gOther.isFSRTight()) maskTight = true;
      }
    }
      
    bool isSelStd = gc.isFSRTight() && !maskTight;  // H4L selection
    bool isSelNano = gc.isFSRTight() && !maskLoose; // nanoAOD selection + a posteriori tighter cut (may be masked by another loose)
      
    bool mu_ele_overlap = false;
    if (gc.closestMu!=0 && gc.closestEle!=0 && gc.isFSRLoose_mu() && gc.isFSRLoose_ele()) { // potential overlap if FSR are collected by mu and ele independently
      mu_ele_overlap = true;
    }


    float stealer_pT = -1; // pT of lepton that would steal this FSR
    int stealer_ID = 0;  // ID of lepton that would steal this FSR
    float stealer_DR = -1.; // DR of lepton that would steal this FSR
    // STEALING: loop over all muons and electrons with nanoAOD presel. 
    // check if muon that is not the one we are considering (by delta R, delta-pT) is closer
    for (unsigned int j = 0; j< muonAllHandle->size(); ++j){
      const pat::Muon*  m = &(*((*muonAllHandle)[j].get()));
      if (m->pdgId()==lID && ROOT::Math::VectorUtil::DeltaR(m->momentum(), gc.closestLep->momentum())<0.1 && abs(m->pt()-lpT)<1.) {
	// FIXME check
	continue;
      }
      float this_DR = ROOT::Math::VectorUtil::DeltaR(m->momentum(), gc.g->momentum());
      if ( this_DR < gc.dR()) {
	stealer_ID = m->pdgId();
	stealer_pT = m->pt();
	stealer_DR = this_DR;
	cout << "STEAL1 " << stealer_ID << " " << stealer_pT << " " << stealer_DR << endl;
      }
    }
    
    //FIXME do the same for electrons;
    //all slimmed electrions above 5GeV, cf: https://github.com/cms-sw/cmssw/blob/CMSSW_11_2_X/PhysicsTools/NanoAOD/python/electrons_cff.py#L339    
    for (unsigned int j = 0; j< electronAllHandle->size(); ++j){
      const pat::Electron* e = &((*electronAllHandle)[j]);
      if (e->pdgId()==lID && ROOT::Math::VectorUtil::DeltaR(e->momentum(), gc.closestLep->momentum())<0.1 && abs(e->pt()-lpT)<1.) {
	// FIXME check
	continue;
      }
      float this_DR = ROOT::Math::VectorUtil::DeltaR(e->momentum(), gc.g->momentum());
      if ( this_DR < gc.dR()) {
	stealer_ID = e->pdgId();
	stealer_pT = e->pt();
	stealer_DR = this_DR;
	cout << "STEAL2 " << stealer_ID << " " << stealer_pT << " " << stealer_DR << endl;
      }
    }
    



    double pTGen = -1.;
    double etaGen = 0;
    double phiGen = 0.;
    if (gc.igen>=0) {
      pTGen = genFSR[gc.igen]->pt();
      etaGen = genFSR[gc.igen]->eta();
      phiGen = genFSR[gc.igen]->phi();
    }

    const reco::GenParticle* gp =  gc.genLepton();
    int lParentType = mch.getParentCode(gp);


    cout << "FSR:" 
	 << iEvent.id().run() << ":"
	 << iEvent.id().luminosityBlock() << ":"
	 << iEvent.id().event() << " "
	 << gc.g->pt() << " "      // reco FSR pT
	 << gc.g->eta() << " " 
	 << gc.g->phi() << " "
	 << gc.gRelIso << " "      // photon iso
	 << gc.SCVetoTight << " "  // pass tighter SC veto (nanoAOD style)
	 << gc.vetoer_pT << " "       // pT of the electron causing a SC tight veto
	 << lID << " "             // ID of the reco l the photon is associated to 
	 << lpT << " "             // pT of the reco l the photon is associated to      
	 << gc.isLepTight() << " " // lepton passes tight selection
	 << gc.dR() << " "         // DR reco FSR vs reco lep
	 << gc.dRET2() << " "      // DRET2 reco FSR vs reco lep      
	 << gc.isFSRLoose() << " " // FSR passes loose (nanHoAOD) sel
	 << gc.isFSRTight() << " " // FSR passes H4l sel
	 << maskLoose << " "       // this photon is masked by another photon passing loose cuts 
	 << maskTight << " "       // this photon is masked by another  photon passing tight cuts
	 << isSelStd << " "        // Selected by H4l sel (= FSRTight && !maskTight)
	 << isSelNano << " "       // Would be selected with tight cuts over nanoAOD (= FSRTight && !maskLoose)
	 << mu_ele_overlap << " "  // Would be associated to both an ele and a mu
	 << stealer_ID << " " 
	 << stealer_pT << " " 
	 << stealer_DR << " " 
	 << pTGen << " "           // pT of gen FSR, or -1 if fake
	 << etaGen << " "
	 << phiGen << " "
	 << genFinalState << " "   
	 << lParentType            // ID of the lepton's original ancestor
      // 	 << dRGenVsReco << " " // dR of gen FSR to reco FSR, if not fake
      // 	 << Nvtx << " "
      // 	 << NObsInt << " "
      // 	 << NTrueInt << " "
	 << endl; 
  }
    



}







#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(FSRStandaloneAnalyzer);

