//
//

/**
  @file CLHEPHitFitTranslator.cc

  @brief Translator for four-vector classes in CLHEP.
  Contains specific implementation for the following classes:
  - LeptonTranslatorBase<CLHEP::HepLorenzVector>
  - JetTranslatorBase<CLHEP::HepLorenzVector>
  - METTranslatorBase<CLHEP::HepLorenzVector>

  @author Haryo Sumowidagdo <Suharyo.Sumowidagdo@cern.ch>

  @date Thu Jul 29 17:54:01 CEST 2010

  @version $Id: CLHEPHitFitTranslator.cc,v 1.4 2010/08/06 22:00:44 haryo Exp $
 */

#include <cstdlib>
#include "CLHEP/Vector/LorentzVector.h"

#include "TopQuarkAnalysis/TopHitFit/interface/LeptonTranslatorBase.h"
#include "TopQuarkAnalysis/TopHitFit/interface/JetTranslatorBase.h"
#include "TopQuarkAnalysis/TopHitFit/interface/METTranslatorBase.h"

namespace hitfit {

using std::string;

template<>
LeptonTranslatorBase<CLHEP::HepLorentzVector>::LeptonTranslatorBase()
{
  string resolution_filename = string(getenv("CMSSW_BASE")) + string("/src/TopQuarkAnalysis/TopHitFit/data/exampleElectronResolution.txt");
  resolution_ = EtaDepResolution(resolution_filename);
} // LeptonTranslatorBase<CLHEP::HepLorentzVector>::LeptonTranslatorBase()


template<>
LeptonTranslatorBase<CLHEP::HepLorentzVector>::LeptonTranslatorBase(const string& ifile)
{
  string resolution_filename = ifile;
  if (ifile.empty()) resolution_filename = string(getenv("CMSSW_BASE")) + string("/src/TopQuarkAnalysis/TopHitFit/data/exampleElectronResolution.txt");

  resolution_ = EtaDepResolution(resolution_filename);
} // LeptonTranslatorBase<CLHEP::HepLorentzVector>::LeptonTranslatorBase(const string& ifile)


template<>
LeptonTranslatorBase<CLHEP::HepLorentzVector>::~LeptonTranslatorBase()
{
}


template<>
Lepjets_Event_Lep
LeptonTranslatorBase<CLHEP::HepLorentzVector>::operator()(const CLHEP::HepLorentzVector& lepton, int type /* = hitfit::lepton_label */)
{
  Fourvec p(lepton.px(),lepton.py(),lepton.pz(),lepton.e());

  double lepton_eta(lepton.eta());
  Vector_Resolution lepton_resolution = resolution_.GetResolution(lepton_eta);
  return Lepjets_Event_Lep(p, type, lepton_resolution);
} // Lepjets_Event_Lep LeptonTranslatorBase<CLHEP::HepLorentzVector>::operator()(const CLHEP::HepLorentzVector& lepton)


template<>
const EtaDepResolution&
LeptonTranslatorBase<CLHEP::HepLorentzVector>::resolution() const
{
  return resolution_;
}


template<>
bool
LeptonTranslatorBase<CLHEP::HepLorentzVector>::CheckEta(const CLHEP::HepLorentzVector& lepton) const
{
  return resolution_.CheckEta(lepton.eta());
}


template<>
JetTranslatorBase<CLHEP::HepLorentzVector>::JetTranslatorBase()
{
  string CMSSW_BASE(getenv("CMSSW_BASE"));

  string udsc_resolution_filename = CMSSW_BASE + string("/src/TopQuarkAnalysis/TopHitFit/data/exampleJetResolution.txt");
  const string& b_resolution_filename = udsc_resolution_filename;

  udscResolution_ = EtaDepResolution(udsc_resolution_filename);
  bResolution_    = EtaDepResolution(b_resolution_filename);
} // JetTranslatorBase<CLHEP::HepLorentzVector>::JetTranslatorBase()


template<>
JetTranslatorBase<CLHEP::HepLorentzVector>::JetTranslatorBase(const string& udscFile, const string& bFile)
{
  string CMSSW_BASE(getenv("CMSSW_BASE"));

  string udsc_resolution_filename = udscFile;
  if (udscFile.empty()) udsc_resolution_filename = CMSSW_BASE + string("/src/TopQuarkAnalysis/TopHitFit/data/exampleJetResolution.txt");

  string b_resolution_filename = bFile;
  if (bFile.empty()) b_resolution_filename = CMSSW_BASE + string("/src/TopQuarkAnalysis/TopHitFit/data/exampleJetResolution.txt");

  udscResolution_ = EtaDepResolution(udsc_resolution_filename);
  bResolution_    = EtaDepResolution(b_resolution_filename);
} // JetTranslatorBase<CLHEP::HepLorentzVector>::JetTranslatorBase(const string& udscFile,const string& bFile)


template<>
JetTranslatorBase<CLHEP::HepLorentzVector>::~JetTranslatorBase()
{
} // JetTranslatorBase<CLHEP::HepLorentzVector>::~JetTranslatorBase()


template<>
Lepjets_Event_Jet
JetTranslatorBase<CLHEP::HepLorentzVector>::operator()(const CLHEP::HepLorentzVector& jet, int type /* = hitfit::unknown_label */)
{
  Fourvec p(jet.px(),jet.py(),jet.pz(),jet.e());

  double            jet_eta        = jet.eta();

  Vector_Resolution jet_resolution;
  if (type == hitfit::lepb_label||type == hitfit::hadb_label||type== hitfit::higgs_label)
    jet_resolution = bResolution_.GetResolution(jet_eta);
  else
    jet_resolution = udscResolution_.GetResolution(jet_eta);

  return Lepjets_Event_Jet(p, type, jet_resolution);
} // Lepjets_Event_Jet JetTranslatorBase<CLHEP::HepLorentzVector>::operator()(const CLHEP::HepLorentzVector& j,int type)


template<>
const EtaDepResolution&
JetTranslatorBase<CLHEP::HepLorentzVector>::udscResolution() const
{
  return udscResolution_;
}


template<>
const EtaDepResolution&
JetTranslatorBase<CLHEP::HepLorentzVector>::bResolution() const
{
  return bResolution_;
}


template<>
bool
JetTranslatorBase<CLHEP::HepLorentzVector>::CheckEta(const CLHEP::HepLorentzVector& jet) const
{
  return udscResolution_.CheckEta(jet.eta()) && bResolution_.CheckEta(jet.eta());
}


template<>
METTranslatorBase<CLHEP::HepLorentzVector>::METTranslatorBase()
{
  resolution_ = Resolution(string("0,0,12"));
} // METTranslatorBase<CLHEP::HepLorentzVector>::METTranslatorBase()


template<>
METTranslatorBase<CLHEP::HepLorentzVector>::METTranslatorBase(const string& ifile)
{
  const Defaults_Text defs(ifile);
  string resolution_string(defs.get_string("met_resolution"));
  resolution_ = Resolution(resolution_string);
} // METTranslatorBase<CLHEP::HepLorentzVector>::METTranslatorBase(const string& ifile)


template<>
METTranslatorBase<CLHEP::HepLorentzVector>::~METTranslatorBase()
{
} // METTranslatorBase<CLHEP::HepLorentzVector>::~METTranslatorBase()


template<>
Fourvec
METTranslatorBase<CLHEP::HepLorentzVector>::operator()(const CLHEP::HepLorentzVector& m)
{
  return Fourvec (m.px(),m.py(),0.0,m.e());
} // Fourvec METTranslatorBase<CLHEP::HepLorentzVector>::operator()(const CLHEP::HepLorentzVector& m)



template<>
Resolution
METTranslatorBase<CLHEP::HepLorentzVector>::KtResolution(const CLHEP::HepLorentzVector& m) const
{
  return resolution_;
} // Resolution METTranslatorBase<CLHEP::HepLorentzVector>::KtResolution(const CLHEP::HepLorentzVector& m)



template<>
Resolution
METTranslatorBase<CLHEP::HepLorentzVector>::METResolution(const CLHEP::HepLorentzVector& m) const
{
  return KtResolution(m);
} // Resolution METTranslatorBase<CLHEP::HepLorentzVector>::METResolution(const CLHEP::HepLorentzVector& m)

} // namespace hitfit
