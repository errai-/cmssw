/**
  @file PatMuonHitFitTranslator.cc

  @brief Specialization of template class LeptonTranslatorBase in the
  package HitFit for pat::Muon.

  @author Haryo Sumowidagdo <Suharyo.Sumowidagdo@cern.ch>

  @par Created
  Sat Jun 27 17:49:15 2009 UTC

  @version $Id: PatMuonHitFitTranslator.cc,v 1.8 2010/08/06 22:03:03 haryo Exp $
 */

#include "TopQuarkAnalysis/TopHitFit/interface/LeptonTranslatorBase.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

namespace hitfit {

using std::string;

template<>
LeptonTranslatorBase<pat::Muon>::LeptonTranslatorBase()
{
  string resolution_filename = string(getenv("CMSSW_BASE")) + string("/src/TopQuarkAnalysis/PatHitFit/data/exampleMuonResolution.txt");
  resolution_ = EtaDepResolution(resolution_filename);
} // LeptonTranslatorBase<pat::Muon>::LeptonTranslatorBase()


template<>
LeptonTranslatorBase<pat::Muon>::LeptonTranslatorBase(const string& ifile)
{
  string resolution_filename = ifile ;
  if (ifile.empty()) resolution_filename = string(getenv("CMSSW_BASE")) + string("/src/TopQuarkAnalysis/PatHitFit/data/exampleMuonResolution.txt");

  resolution_ = EtaDepResolution(resolution_filename);
} // LeptonTranslatorBase<pat::Muon>::LeptonTranslatorBase(const string& s)


template<>
LeptonTranslatorBase<pat::Muon>::~LeptonTranslatorBase()
{
} // LeptonTranslatorBase<pat::Muon>::~LeptonTranslatorBase()


template<>
Lepjets_Event_Lep
LeptonTranslatorBase<pat::Muon>::operator()(const pat::Muon& lepton, int type /*= hitfit::lepton_label */)
{
  Fourvec p(lepton.px(),lepton.py(),lepton.pz(),lepton.energy());

  double            muon_eta        = lepton.eta();
  Vector_Resolution muon_resolution = resolution_.GetResolution(muon_eta);

  return Lepjets_Event_Lep(p, muon_label, muon_resolution);
} // Lepjets_Event_Lep LeptonTranslatorBase<pat::Muon>::operator()


template<>
const EtaDepResolution&
LeptonTranslatorBase<pat::Muon>::resolution() const
{
  return resolution_;
}

template<>
bool
LeptonTranslatorBase<pat::Muon>::CheckEta(const pat::Muon& lepton) const
{
  return resolution_.CheckEta(lepton.eta());
}

} // namespace hitfit
