/**
  @file PatElectronHitFitTranslator.cc

  @brief Specialization of template class LeptonTranslatorBase in the
  package HitFit for pat::Electron.

  @author Haryo Sumowidagdo <Suharyo.Sumowidagdo@cern.ch>

  @par Created
  Sat Jun 27 17:49:06 2009 UTC

  @version $Id: PatElectronHitFitTranslator.cc,v 1.8 2010/08/06 22:02:52 haryo Exp $
 */

#include "TopQuarkAnalysis/TopHitFit/interface/LeptonTranslatorBase.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

namespace hitfit {

using std::string;

template<>
LeptonTranslatorBase<pat::Electron>::LeptonTranslatorBase()
{
  string resolution_filename = string(getenv("CMSSW_BASE")) + string("/src/TopQuarkAnalysis/PatHitFit/data/exampleElectronResolution.txt");
  resolution_ = EtaDepResolution(resolution_filename);
} // LeptonTranslatorBase<pat::Electron>::LeptonTranslatorBase()


template<>
LeptonTranslatorBase<pat::Electron>::LeptonTranslatorBase(const string& ifile)
{
  string resolution_filename = ifile ;
  if (ifile.empty()) resolution_filename = string(getenv("CMSSW_BASE")) + string("/src/TopQuarkAnalysis/PatHitFit/data/exampleElectronResolution.txt");

  resolution_ = EtaDepResolution(resolution_filename);
} // LeptonTranslatorBase<pat::Electron>::LeptonTranslatorBase(const string& ifile)


template<>
LeptonTranslatorBase<pat::Electron>::~LeptonTranslatorBase()
{
}


template<>
Lepjets_Event_Lep
LeptonTranslatorBase<pat::Electron>::operator()(const pat::Electron& lepton, int type /* = hitfit::lepton_label */)
{
  Fourvec p(lepton.px(),lepton.py(),lepton.pz(),lepton.energy());

  double            electron_eta        = lepton.superCluster()->eta();
  Vector_Resolution electron_resolution = resolution_.GetResolution(electron_eta);

  return Lepjets_Event_Lep(p, electron_label, electron_resolution);
} // Lepjets_Event_Lep LeptonTranslatorBase<pat::Electron>::operator()

template<>
const EtaDepResolution&
LeptonTranslatorBase<pat::Electron>::resolution() const
{
  return resolution_;
}

template<>
bool
LeptonTranslatorBase<pat::Electron>::CheckEta(const pat::Electron& lepton) const
{
  double            electron_eta        = lepton.superCluster()->eta();
  return resolution_.CheckEta(electron_eta);
}

} // namespace hitfit
