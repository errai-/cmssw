/**
  @file PatJetHitFitTranslator.cc

  @brief Specialization of template class JetTranslatorBase in the
  package HitFit

  @author Haryo Sumowidagdo <Suharyo.Sumowidagdo@cern.ch>

  @par Created
  Sat Jun 27 17:49:21 2009 UTC

  @version $Id: PatJetHitFitTranslator.cc,v 1.2 2011/05/26 12:57:18 mseidel Exp $
 */


#include <TopQuarkAnalysis/TopHitFit/interface/JetTranslatorBase.h>
#include <DataFormats/PatCandidates/interface/Jet.h>

namespace hitfit {

using std::string;

template<>
JetTranslatorBase<pat::Jet>::JetTranslatorBase()
{
  string CMSSW_BASE(getenv("CMSSW_BASE"));
  string resolution_filename = CMSSW_BASE + string("/src/TopQuarkAnalysis/PatHitFit/data/exampleJetResolution.txt");
  udscResolution_ = EtaDepResolution(resolution_filename);
  bResolution_    = EtaDepResolution(resolution_filename);
  jetCorrectionLevel_ = "L7Parton";
  jes_            = 1.0;
  jesB_           = 1.0;
} // JetTranslatorBase<pat::Jet>::JetTranslatorBase()


template<>
JetTranslatorBase<pat::Jet>::JetTranslatorBase(const string& udscFile, const string& bFile)
{
  string CMSSW_BASE(getenv("CMSSW_BASE"));

  string udscResolution_filename = udscFile;
  if (udscFile.empty()) udscResolution_filename = CMSSW_BASE + string("/src/TopQuarkAnalysis/PatHitFit/data/exampleJetResolution.txt");

  string bResolution_filename = bFile;
  if (bFile.empty()) bResolution_filename = CMSSW_BASE + string("/src/TopQuarkAnalysis/PatHitFit/data/exampleJetResolution.txt");

  udscResolution_ = EtaDepResolution(udscResolution_filename);
  bResolution_    = EtaDepResolution(bResolution_filename);
  jetCorrectionLevel_ = "L7Parton";
  jes_            = 1.0;
  jesB_           = 1.0;
} // JetTranslatorBase<pat::Jet>::JetTranslatorBase(const string& ifile)


template<>
JetTranslatorBase<pat::Jet>::JetTranslatorBase(const string& udscFile, const string& bFile, const string& jetCorrectionLevel,
                                               double jes, double jesB)
{
  string CMSSW_BASE(getenv("CMSSW_BASE"));

  string udscResolution_filename = udscFile;
  if (udscFile.empty()) udscResolution_filename = CMSSW_BASE + string("/src/TopQuarkAnalysis/TopHitFit/data/resolution/tqafUdscJetResolution.txt");

  string bResolution_filename = bFile;
  if (bFile.empty()) bResolution_filename = CMSSW_BASE + string("/src/TopQuarkAnalysis/TopHitFit/data/resolution/tqafBJetResolution.txt");

  udscResolution_ = EtaDepResolution(udscResolution_filename);
  bResolution_    = EtaDepResolution(bResolution_filename);
  jetCorrectionLevel_ = jetCorrectionLevel;
  jes_            = jes;
  jesB_           = jesB;
} // JetTranslatorBase<pat::Jet>::JetTranslatorBase(const string& ifile)


template<>
JetTranslatorBase<pat::Jet>::~JetTranslatorBase()
{
} // JetTranslatorBase<pat::Jet>::~JetTranslatorBase()


template<>
Lepjets_Event_Jet
JetTranslatorBase<pat::Jet>::operator()(const pat::Jet& jet, int type /*= hitfit::unknown_label */)
{
  Fourvec p;

  double jet_eta = jet.eta();

  //if (jet.isPFJet()) do nothing at the moment!
  if (jet.isCaloJet()) jet_eta = static_cast<const reco::CaloJet*>(jet.originalObject())->detectorP4().eta();

  Vector_Resolution jet_resolution;

  if (type == hitfit::hadb_label || type == hitfit::lepb_label || type == hitfit::higgs_label) {
    jet_resolution = bResolution_.GetResolution(jet_eta);
    pat::Jet bPartonCorrJet(jet.correctedJet(jetCorrectionLevel_,"BOTTOM"));
    bPartonCorrJet.scaleEnergy(jesB_);
    p = Fourvec(bPartonCorrJet.px(),bPartonCorrJet.py(),bPartonCorrJet.pz(),bPartonCorrJet.energy());
  } else {
    jet_resolution = udscResolution_.GetResolution(jet_eta);
    pat::Jet udsPartonCorrJet(jet.correctedJet(jetCorrectionLevel_,"UDS"));
    udsPartonCorrJet.scaleEnergy(jes_);
    p = Fourvec(udsPartonCorrJet.px(),udsPartonCorrJet.py(),udsPartonCorrJet.pz(),udsPartonCorrJet.energy());
  }

  return Lepjets_Event_Jet(p, type, jet_resolution);
} // Lepjets_Event_Jet JetTranslatorBase<pat::Jet>::operator()(const pat::Jet& j,int type)


template<>
const EtaDepResolution&
JetTranslatorBase<pat::Jet>::udscResolution() const
{
  return udscResolution_;
}


template<>
const EtaDepResolution&
JetTranslatorBase<pat::Jet>::bResolution() const
{
  return bResolution_;
}


template<>
bool
JetTranslatorBase<pat::Jet>::CheckEta(const pat::Jet& jet) const
{
  double jet_eta = jet.eta();

  //if (jet.isPFJet()) do nothing at the moment!
  if (jet.isCaloJet()) jet_eta = static_cast<const reco::CaloJet*>(jet.originalObject())->detectorP4().eta();
  return bResolution_.CheckEta(jet_eta) and udscResolution_.CheckEta(jet_eta);
}

} // namespace hitfit
