#ifndef TtSemiLepHitFitProducer_h
#define TtSemiLepHitFitProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/JetMCUtils/interface/combination.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLepEvtPartons.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h"

#include "TopQuarkAnalysis/TopHitFit/interface/RunHitFit.h"
#include "TopQuarkAnalysis/TopHitFit/interface/LeptonTranslatorBase.h"
#include "TopQuarkAnalysis/TopHitFit/interface/JetTranslatorBase.h"
#include "TopQuarkAnalysis/TopHitFit/interface/METTranslatorBase.h"

using std::string;
using std::vector;

template <typename LeptonCollection>
class TtSemiLepHitFitProducer : public edm::EDProducer {

 public:

  explicit TtSemiLepHitFitProducer(const edm::ParameterSet&);
  ~TtSemiLepHitFitProducer() override;

 private:
  // produce
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<vector<pat::Jet> > jetsToken_;
  edm::EDGetTokenT<LeptonCollection>  lepsToken_;
  edm::EDGetTokenT<vector<pat::MET> > metsToken_;

  /// maximal number of jets (-1 possible to indicate 'all')
  int maxNJets_;
  /// maximal number of combinations to be written to the event
  int maxNComb_;

  /// maximum eta value for muons, needed to limited range in which resolutions are provided
  double maxEtaMu_;
  /// maximum eta value for electrons, needed to limited range in which resolutions are provided
  double maxEtaEle_;
  /// maximum eta value for jets, needed to limited range in which resolutions are provided
  double maxEtaJet_;

  /// input tag for b-tagging algorithm
  string bTagAlgo_;
  /// min value of bTag for a b-jet
  double minBTagValueBJet_;
  /// max value of bTag for a non-b-jet
  double maxBTagValueNonBJet_;
  /// switch to tell whether to use b-tagging or not
  bool useBTag_;

  /// constraints
  double mW_;
  double mTop_;

  /// jet correction level
  string jetCorrectionLevel_;

  /// jet energy scale
  double jes_;
  double jesB_;

  struct FitResult {
    int Status;
    double Chi2;
    double Prob;
    double MT;
    double SigMT;
    pat::Particle HadB;
    pat::Particle HadP;
    pat::Particle HadQ;
    pat::Particle LepB;
    pat::Particle LepL;
    pat::Particle LepN;
    vector<int> JetCombi;
    bool operator< (const FitResult& rhs) { return Chi2 < rhs.Chi2; };
  };

  typedef hitfit::RunHitFit<pat::Electron,pat::Muon,pat::Jet,pat::MET> PatHitFit;

  edm::FileInPath hitfitDefault_;
  edm::FileInPath hitfitElectronResolution_;
  edm::FileInPath hitfitMuonResolution_;
  edm::FileInPath hitfitUdscJetResolution_;
  edm::FileInPath hitfitBJetResolution_;
  edm::FileInPath hitfitMETResolution_;

  hitfit::LeptonTranslatorBase<pat::Electron> electronTranslator_;
  hitfit::LeptonTranslatorBase<pat::Muon>     muonTranslator_;
  hitfit::JetTranslatorBase<pat::Jet>         jetTranslator_;
  hitfit::METTranslatorBase<pat::MET>         metTranslator_;

  PatHitFit* HitFit;
};

template<typename LeptonCollection>
TtSemiLepHitFitProducer<LeptonCollection>::TtSemiLepHitFitProducer(const edm::ParameterSet& cfg):
  jetsToken_               (consumes<vector<pat::Jet> >    (cfg.getParameter<edm::InputTag>("jets"))),
  lepsToken_               (consumes<LeptonCollection>     (cfg.getParameter<edm::InputTag>("leps"))),
  metsToken_               (consumes<vector<pat::MET> >    (cfg.getParameter<edm::InputTag>("mets"))),
  maxNJets_                (cfg.getParameter<int>          ("maxNJets"            )),
  maxNComb_                (cfg.getParameter<int>          ("maxNComb"            )),
  bTagAlgo_                (cfg.getParameter<string>       ("bTagAlgo"            )),
  minBTagValueBJet_        (cfg.getParameter<double>       ("minBDiscBJets"       )),
  maxBTagValueNonBJet_     (cfg.getParameter<double>       ("maxBDiscLightJets"   )),
  useBTag_                 (cfg.getParameter<bool>         ("useBTagging"         )),
  mW_                      (cfg.getParameter<double>       ("mW"                  )),
  mTop_                    (cfg.getParameter<double>       ("mTop"                )),
  jetCorrectionLevel_      (cfg.getParameter<string>       ("jetCorrectionLevel"  )),
  jes_                     (cfg.getParameter<double>       ("jes"                 )),
  jesB_                    (cfg.getParameter<double>       ("jesB"                )),

  // The following five initializers read the config parameters for the
  // ASCII text files which contains the physics object resolutions.
  hitfitDefault_           (cfg.getUntrackedParameter<edm::FileInPath>(string("hitfitDefault")           , edm::FileInPath(string("TopQuarkAnalysis/TopHitFit/data/setting/RunHitFitConfiguration.txt")))),
  hitfitElectronResolution_(cfg.getUntrackedParameter<edm::FileInPath>(string("hitfitElectronResolution"), edm::FileInPath(string("TopQuarkAnalysis/TopHitFit/data/resolution/tqafElectronResolution.txt")))),
  hitfitMuonResolution_    (cfg.getUntrackedParameter<edm::FileInPath>(string("hitfitMuonResolution")    , edm::FileInPath(string("TopQuarkAnalysis/TopHitFit/data/resolution/tqafMuonResolution.txt")))),
  hitfitUdscJetResolution_ (cfg.getUntrackedParameter<edm::FileInPath>(string("hitfitUdscJetResolution") , edm::FileInPath(string("TopQuarkAnalysis/TopHitFit/data/resolution/tqafUdscJetResolution.txt")))),
  hitfitBJetResolution_    (cfg.getUntrackedParameter<edm::FileInPath>(string("hitfitBJetResolution")    , edm::FileInPath(string("TopQuarkAnalysis/TopHitFit/data/resolution/tqafBJetResolution.txt")))),
  hitfitMETResolution_     (cfg.getUntrackedParameter<edm::FileInPath>(string("hitfitMETResolution")     , edm::FileInPath(string("TopQuarkAnalysis/TopHitFit/data/resolution/tqafKtResolution.txt")))),

  // The following four initializers instantiate the translator between PAT objects
  // and HitFit objects using the ASCII text files which contains the resolutions.
  electronTranslator_(hitfitElectronResolution_.fullPath()),
  muonTranslator_    (hitfitMuonResolution_.fullPath()),
  jetTranslator_     (hitfitUdscJetResolution_.fullPath(), hitfitBJetResolution_.fullPath(), jetCorrectionLevel_, jes_, jesB_),
  metTranslator_     (hitfitMETResolution_.fullPath())
{
  // Create an instance of RunHitFit and initialize it.
  HitFit = new PatHitFit(electronTranslator_,
                         muonTranslator_,
                         jetTranslator_,
                         metTranslator_,
                         hitfitDefault_.fullPath(),
                         mW_,
                         mW_,
                         mTop_);

  maxEtaMu_  = 2.4;
  maxEtaEle_ = 2.5;
  maxEtaJet_ = 2.5;

  edm::LogVerbatim( "TopHitFit" )
    << "\n"
    << "+++++++++++ TtSemiLepHitFitProducer ++++++++++++ \n"
    << " Due to the eta ranges for which resolutions     \n"
    << " are provided in                                 \n"
    << " TopQuarkAnalysis/TopHitFit/data/resolution/     \n"
    << " so far, the following cuts are currently        \n"
    << " implemented in the TtSemiLepHitFitProducer:     \n"
    << " |eta(muons    )| <= " << maxEtaMu_  <<        " \n"
    << " |eta(electrons)| <= " << maxEtaEle_ <<        " \n"
    << " |eta(jets     )| <= " << maxEtaJet_ <<        " \n"
    << "++++++++++++++++++++++++++++++++++++++++++++++++ \n";

  produces< vector<pat::Particle> >("PartonsHadP");
  produces< vector<pat::Particle> >("PartonsHadQ");
  produces< vector<pat::Particle> >("PartonsHadB");
  produces< vector<pat::Particle> >("PartonsLepB");
  produces< vector<pat::Particle> >("Leptons");
  produces< vector<pat::Particle> >("Neutrinos");

  produces< vector<vector<int> > >();
  produces< vector<double> >("Chi2");
  produces< vector<double> >("Prob");
  produces< vector<double> >("MT");
  produces< vector<double> >("SigMT");
  produces< vector<int> >("Status");
  produces< int >("NumberOfConsideredJets");
}

template<typename LeptonCollection>
TtSemiLepHitFitProducer<LeptonCollection>::~TtSemiLepHitFitProducer()
{
  delete HitFit;
}

template<typename LeptonCollection>
void TtSemiLepHitFitProducer<LeptonCollection>::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  std::unique_ptr< vector<pat::Particle> > pPartonsHadP( new vector<pat::Particle> );
  std::unique_ptr< vector<pat::Particle> > pPartonsHadQ( new vector<pat::Particle> );
  std::unique_ptr< vector<pat::Particle> > pPartonsHadB( new vector<pat::Particle> );
  std::unique_ptr< vector<pat::Particle> > pPartonsLepB( new vector<pat::Particle> );
  std::unique_ptr< vector<pat::Particle> > pLeptons    ( new vector<pat::Particle> );
  std::unique_ptr< vector<pat::Particle> > pNeutrinos  ( new vector<pat::Particle> );

  std::unique_ptr< vector<vector<int> > > pCombi ( new vector<vector<int> > );
  std::unique_ptr< vector<double>       > pChi2  ( new vector<double> );
  std::unique_ptr< vector<double>       > pProb  ( new vector<double> );
  std::unique_ptr< vector<double>       > pMT    ( new vector<double> );
  std::unique_ptr< vector<double>       > pSigMT ( new vector<double> );
  std::unique_ptr< vector<int>          > pStatus( new vector<int> );
  std::unique_ptr< int > pJetsConsidered( new int );

  edm::Handle<vector<pat::Jet> > jets;
  evt.getByToken(jetsToken_, jets);

  edm::Handle<vector<pat::MET> > mets;
  evt.getByToken(metsToken_, mets);

  edm::Handle<LeptonCollection> leps;
  evt.getByToken(lepsToken_, leps);

  // -----------------------------------------------------
  // skip events with no appropriate lepton candidate in
  // or empty MET or less jets than partons
  // -----------------------------------------------------

  const unsigned int nPartons = 4;

  // Clear the internal state
  HitFit->clear();

  // Add lepton into HitFit
  bool foundLepton = false;
  if (!leps->empty()) {
    double maxEtaLep = maxEtaMu_;
    if ( !dynamic_cast<const reco::Muon*>(&((*leps)[0])) ) // assume electron if it is not a muon
      maxEtaLep = maxEtaEle_;
    for (unsigned iLep=0; iLep<(*leps).size(); ++iLep) {
      if (std::abs((*leps)[iLep].eta()) <= maxEtaLep) {
        HitFit->AddLepton((*leps)[iLep]);
        foundLepton = true;
        break;
      }
    }
  }

  // Add jets into HitFit
  unsigned int nJetsFound = 0;
  for(unsigned iJet=0; iJet<(*jets).size() && (int)nJetsFound!=maxNJets_; ++iJet) {
    double jet_abseta = std::abs((*jets)[iJet].eta());
    if ( (*jets)[iJet].isCaloJet() )
      jet_abseta = std::abs(((reco::CaloJet*) ((*jets)[iJet]).originalObject())->detectorP4().eta());
    if (jet_abseta <= maxEtaJet_) {
      HitFit->AddJet((*jets)[iJet]);
      ++nJetsFound;
    }
  }
  *pJetsConsidered = nJetsFound;

  //
  // R U N   H I T F I T
  //
  // Run the kinematic fit and get how many permutations are possible
  // in the fit
  std::list<FitResult> FitResultList;
  if (foundLepton and !mets->empty() and nJetsFound>=nPartons) {
    // Add missing transverse energy into HitFit
    HitFit->SetMet((*mets)[0]);

    // Number of all permutations of the event
    HitFit->FitAllPermutation();

    // Get the fit results for all permutations
    vector<hitfit::Fit_Result> hitFitResult = HitFit->GetFitAllPermutation();

    //
    // BEGIN PART WHICH EXTRACTS INFORMATION FROM HITFIT
    //

    // Loop over all permutations and extract the information, save into a vector that is later sorted according to chi2
    for (const auto &fit : hitFitResult) {
      /*
        Get jet permutation according to TQAF convention
        11 : leptonic b
        12 : hadronic b
        13 : hadronic W
        14 : hadronic W
      */
      vector<int> hitCombi(4);

      // Get the number of jets and loop over the jets
      const size_t nHitFitJet = fit.njets();
      for (size_t jet = 0, nJets = fit.njets(); jet < ; ++jet) {
        int jet_type = fit.jet(jet).type();

        switch (jet_type) {
          case 11: hitCombi[TtSemiLepEvtPartons::LepB     ] = jet;
            break;
          case 12: hitCombi[TtSemiLepEvtPartons::HadB     ] = jet;
            break;
          case 13: hitCombi[TtSemiLepEvtPartons::LightQ   ] = jet;
            break;
          case 14: hitCombi[TtSemiLepEvtPartons::LightQBar] = jet;
            break;
        }
      }

      // Store the kinematic quantities in the corresponding containers.
      hitfit::Lepjets_Event_Jet hadP_ = fit.jet(hitCombi[TtSemiLepEvtPartons::LightQ   ]);
      hitfit::Lepjets_Event_Jet hadQ_ = fit.jet(hitCombi[TtSemiLepEvtPartons::LightQBar]);
      hitfit::Lepjets_Event_Jet hadB_ = fit.jet(hitCombi[TtSemiLepEvtPartons::HadB     ]);
      hitfit::Lepjets_Event_Jet lepB_ = fit.jet(hitCombi[TtSemiLepEvtPartons::LepB     ]);
      hitfit::Lepjets_Event_Lep lepL_ = fit.lep(0);

      if (hitFitResult[fit].chisq() > 0 and (!useBTag_ // only take into account converged fits use btag information if chosen
                           or (jets->at(hitCombi[TtSemiLepEvtPartons::LightQ   ]).bDiscriminator(bTagAlgo_) < maxBTagValueNonBJet_
                           and jets->at(hitCombi[TtSemiLepEvtPartons::LightQBar]).bDiscriminator(bTagAlgo_) < maxBTagValueNonBJet_
                           and jets->at(hitCombi[TtSemiLepEvtPartons::HadB     ]).bDiscriminator(bTagAlgo_) > minBTagValueBJet_
                           and jets->at(hitCombi[TtSemiLepEvtPartons::LepB     ]).bDiscriminator(bTagAlgo_) > minBTagValueBJet_)))
      {
        FitResult result;
        result.Status = 0;
        result.Chi2 = hitFitResult[fit].chisq();
        result.Prob = exp(-1.0*(hitFitResult[fit].chisq())/2.0);
        result.MT   = hitFitResult[fit].mt();
        result.SigMT= hitFitResult[fit].sigmt();
        result.HadB = pat::Particle(reco::LeafCandidate(0, math::XYZTLorentzVector(hadB_.p().x(), hadB_.p().y(), hadB_.p().z(), hadB_.p().t()), math::XYZPoint()));
        result.HadP = pat::Particle(reco::LeafCandidate(0, math::XYZTLorentzVector(hadP_.p().x(), hadP_.p().y(), hadP_.p().z(), hadP_.p().t()), math::XYZPoint()));
        result.HadQ = pat::Particle(reco::LeafCandidate(0, math::XYZTLorentzVector(hadQ_.p().x(), hadQ_.p().y(), hadQ_.p().z(), hadQ_.p().t()), math::XYZPoint()));
        result.LepB = pat::Particle(reco::LeafCandidate(0, math::XYZTLorentzVector(lepB_.p().x(), lepB_.p().y(), lepB_.p().z(), lepB_.p().t()), math::XYZPoint()));
        result.LepL = pat::Particle(reco::LeafCandidate(0, math::XYZTLorentzVector(lepL_.p().x(), lepL_.p().y(), lepL_.p().z(), lepL_.p().t()), math::XYZPoint()));
        result.LepN = pat::Particle(reco::LeafCandidate(0, math::XYZTLorentzVector(fittedEvent.met().x(), fittedEvent.met().y(), fittedEvent.met().z(), fittedEvent.met().t()), math::XYZPoint()));
        result.JetCombi = hitCombi;

        FitResultList.emplace_back(result);
      }
    }
  }

  // -----------------------------------------------------
  // feed out result
  // starting with the JetComb having the smallest chi2
  // -----------------------------------------------------
  if (FitResultList.size()==0) { // in case no fit results were stored in the list
    // the kinFit getters return empty objects here
    pPartonsHadP->push_back( pat::Particle() );
    pPartonsHadQ->push_back( pat::Particle() );
    pPartonsHadB->push_back( pat::Particle() );
    pPartonsLepB->push_back( pat::Particle() );
    pLeptons    ->push_back( pat::Particle() );
    pNeutrinos  ->push_back( pat::Particle() );
    // indices referring to the jet combination
    vector<int> invalidCombi;
    for (unsigned int i = 0; i < nPartons; ++i) invalidCombi.push_back( -1 );
    pCombi->push_back( invalidCombi );
    // chi2
    pChi2->push_back( -1. );
    // chi2 probability
    pProb->push_back( -1. );
    // fitted top mass
    pMT->push_back( -1. );
    pSigMT->push_back( -1. );
    // status of the fitter
    pStatus->push_back( -1 );
  } else {
    // sort results w.r.t. chi2 values
    FitResultList.sort();

    int iComb = 0;
    for (auto result = FitResultList.begin(); result != FitResultList.end(); ++result) {
      if (maxNComb_>=0 and ++iComb>maxNComb_) break;
      // partons
      pPartonsHadP->push_back( result->HadP );
      pPartonsHadQ->push_back( result->HadQ );
      pPartonsHadB->push_back( result->HadB );
      pPartonsLepB->push_back( result->LepB );
      // lepton
      pLeptons->push_back( result->LepL );
      // neutrino
      pNeutrinos->push_back( result->LepN );
      // indices referring to the jet combination
      pCombi->push_back( result->JetCombi );
      // chi2
      pChi2->push_back( result->Chi2 );
      // chi2 probability
      pProb->push_back( result->Prob );
      // fitted top mass
      pMT->push_back( result->MT );
      pSigMT->push_back( result->SigMT );
      // status of the fitter
      pStatus->push_back( result->Status );
    }
  }
  evt.put(std::move(pCombi));
  evt.put(std::move(pPartonsHadP), "PartonsHadP");
  evt.put(std::move(pPartonsHadQ), "PartonsHadQ");
  evt.put(std::move(pPartonsHadB), "PartonsHadB");
  evt.put(std::move(pPartonsLepB), "PartonsLepB");
  evt.put(std::move(pLeptons    ), "Leptons"    );
  evt.put(std::move(pNeutrinos  ), "Neutrinos"  );
  evt.put(std::move(pChi2       ), "Chi2"       );
  evt.put(std::move(pProb       ), "Prob"       );
  evt.put(std::move(pMT         ), "MT"         );
  evt.put(std::move(pSigMT      ), "SigMT"      );
  evt.put(std::move(pStatus     ), "Status"     );
  evt.put(std::move(pJetsConsidered), "NumberOfConsideredJets");
}

#endif
