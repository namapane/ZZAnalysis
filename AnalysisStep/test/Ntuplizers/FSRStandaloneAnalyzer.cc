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
  double dRMinEle; // distance to closest ele
  double dRMinMu;  // distance to closest mu
  const reco::Candidate* closestEle = nullptr;
  const reco::Candidate* closestMu = nullptr;
  const reco::Candidate* closestLep = nullptr; // closest of the above
  int igen; // index of matched gen FSR photon (if any)
  
  // Parameters, w.r.t. closest lepton
  float dR() {return min(dRMinMu,dRMinEle);}
  float dRET2() {return dR()/g->pt()/g->pt();}
  int lepId() {return closestLep->pdgId();}
  bool isLoose() {return gRelIso<2   && dRET2()<0.05;}  // nanoAOD selection
  bool isTight() {return gRelIso<1.8 && dRET2()<0.012;} // H4l selection

  // w.r.t. closest mu
  bool isLoose_mu() {return closestMu!=nullptr && gRelIso<2   && dRMinMu/g->pt()/g->pt()<0.05;}
  bool isTight_mu() {return closestMu!=nullptr && gRelIso<1.8 && dRMinMu/g->pt()/g->pt()<0.012;}
  
  // w.r.t. closest ele
  bool isLoose_ele() {return closestEle!=nullptr && gRelIso<2   && dRMinEle/g->pt()/g->pt()<0.05;}
  bool isTight_ele() {return closestEle!=nullptr && gRelIso<1.8 && dRMinEle/g->pt()/g->pt()<0.012;}

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
  edm::EDGetTokenT<pat::ElectronCollection> electronToken;
  edm::EDGetTokenT<pat::ElectronCollection> electronVetoToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > genParticleToken;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
//   edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PupInfoToken;
//   edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken;

  bool debug;
};


FSRStandaloneAnalyzer::FSRStandaloneAnalyzer(const edm::ParameterSet& iConfig) :
  pfCandToken(consumes<edm::View<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates"))), // input for photons
  muonToken(consumes<pat::MuonCollection>(edm::InputTag("softMuons"))), 
  electronToken(consumes<pat::ElectronCollection>(edm::InputTag("cleanSoftElectrons"))),
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

  // MC Truth
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByToken(genParticleToken, genParticles);
  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(genInfoToken, genInfo);
  MCHistoryTools mch(iEvent, "", genParticles, genInfo);
  vector<const reco::Candidate *> genFSR = mch.genFSR();
  vector<int> nRecoFSRMatchedToGen(genFSR.size(),0);

  //----------------------
  // Loop on photons
  //----------------------

  vector<FSRCandidate> photons;
  
  typedef map<const reco::Candidate*, vector<size_t> > PhIdxsByLep;
  PhIdxsByLep photonsByLep;


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
    bool SCVeto=false;
    bool SCVetoTight=false; // Tighter (veto also non-SIP electrons)
    if (electronHandle->size()>0) {
      for (unsigned int j = 0; j< electronHandle->size(); ++j){
	const pat::Electron* e = &((*electronHandle)[j]);
	if ((e->associatedPackedPFCandidates()).size()) {
	  edm::RefVector < pat::PackedCandidateCollection > pfcands = e->associatedPackedPFCandidates();
	  for ( auto itr: pfcands ) {
	    if (g.get()==&(*itr)) {
	      SCVetoTight=true;
	      if (e->userFloat("isSIP")) SCVeto=true;
// 	      if (debug) cout << "SC veto (loose only=" << SCVeto << ") " << itr->eta() << " " << itr->phi() << " " 
// 			      << fabs(g->eta() - itr->eta()) << " " << reco::deltaPhi(g->phi(), itr->phi()) << endl;
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
    const reco::Candidate* closestMu = 0;
    const reco::Candidate* closestEle = 0;
    
    // Loop over pat::Muon
    for (unsigned int j = 0; j< muonHandle->size(); ++j){
      const pat::Muon* m = &((*muonHandle)[j]);
      if (! m->userFloat("isSIP")) continue;
      double dR = ROOT::Math::VectorUtil::DeltaR(m->momentum(),g->momentum());
      if(debug) cout << "muon pt = " << m->pt() << " photon pt = " << g->pt() << " dR = " << dR << endl;
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
      if(debug) cout << "ele pt = " << e->pt() << " photon pt = " << g->pt() << " dR = " << dR << endl;
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
      
    int lID = gc.closestLep->pdgId(); // ID of the closest lepton

    // Are there other photons associated to the same lepton?
    vector<size_t> allPhotons = photonsByLep[gc.closestLep];
    bool maskLoose = false; // this photon is masked by another, loose
    bool maskTight = false; // this photon is masked by another, tight
    for (const size_t& iOther : allPhotons){
      FSRCandidate& gOther = photons[iOther];
      if (iOther!=iPhoton && gOther.dRET2()<gc.dRET2()){ // Another photon could be selected for this lepton
	if (gOther.isLoose()) maskLoose = true;
	if (gOther.isTight()) maskTight = true;
      }
    }
      
    bool isSelStd = gc.isTight() && !maskTight;  // H4L selection
    bool isSelNano = gc.isTight() && !maskLoose; // nanoAOD selection + a posteriori tighter cut (may be masked by another loose)
      
    bool mu_ele_overlap = false;
    if (gc.closestMu!=0 && gc.closestEle!=0 && gc.isLoose_mu() && gc.isLoose_ele()) { // potential overlap if FSR are collected by mu and ele independently
      mu_ele_overlap = true;
      cout << "OVERLAP" << endl;
    }
      
    double pTGen = -1.;
    double etaGen = 0;
    double phiGen = 0.;
    if (gc.igen>=0) {
      pTGen = genFSR[gc.igen]->pt();
      etaGen = genFSR[gc.igen]->eta();
      phiGen = genFSR[gc.igen]->phi();
    }

    if (!gc.isLoose()) continue;

    cout << "FSR:" 
	 << iEvent.id().run() << ":"
	 << iEvent.id().luminosityBlock() << ":"
	 << iEvent.id().event() << " "
	 << gc.g->pt() << " "  // reco FSR pT
	 << gc.g->eta() << " " 
	 << gc.g->phi() << " "
	 << gc.gRelIso << " " // photon iso
	 << gc.SCVetoTight << " " // pass tighter SC veto 
	 << lID << " " // this is the ID of the reco l the photon is associated to 
	 << gc.dR() << " "  // DR-reco FSR vs reco lep
	 << gc.dRET2() << " "  // DR-reco FSR vs reco lep      
	 << gc.isLoose() << " "
	 << gc.isTight() << " "
	 << maskLoose << " " // masked by another passing loose sel
	 << maskTight << " " // masked by another passing loose sel
	 << isSelStd << " " // Selected by standard sel
	 << isSelNano << " " // nanoAOD + tight sel
	 << mu_ele_overlap << " " 
      // 	 << dRGenVsReco << " " // dR of gen FSR to reco FSR, if not fake 
	 << pTGen << " "       // pT of gen FSR, or -1 if fake
	 << etaGen << " "
	 << phiGen << " "
      // 	 << Nvtx << " "
      // 	 << NObsInt << " "
      // 	 << NTrueInt << " "
	 << endl; 
  }
    



}







#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(FSRStandaloneAnalyzer);

