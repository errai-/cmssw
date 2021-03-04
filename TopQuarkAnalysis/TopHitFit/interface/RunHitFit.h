/**
  @file RunHitFit.h

  @brief Template class of experiment-independent interface to HitFit.

  @author Haryo Sumowidagdo <Suharyo.Sumowidagdo@cern.ch>

  @par Creation date:
  May 2009.

  @par Modification History:
  Nov 2009: Haryo Sumowidagdo <Suharyo.Sumowidagdo@cern.ch>:
  Add doxygen tags for automatic generation of documentation.
 */

#ifndef HITFIT_RUNHITFIT_H
#define HITFIT_RUNHITFIT_H

#include <algorithm>

#include "TopQuarkAnalysis/TopHitFit/interface/Defaults_Text.h"
#include "TopQuarkAnalysis/TopHitFit/interface/LeptonTranslatorBase.h"
#include "TopQuarkAnalysis/TopHitFit/interface/JetTranslatorBase.h"
#include "TopQuarkAnalysis/TopHitFit/interface/METTranslatorBase.h"
#include "TopQuarkAnalysis/TopHitFit/interface/Fit_Result.h"
#include "TopQuarkAnalysis/TopHitFit/interface/Lepjets_Event.h"
#include "TopQuarkAnalysis/TopHitFit/interface/Top_Fit.h"

// Explanation about the MIN/MAX definitions:
//
// For a given number of jets, there is a corresponding number of
// permutations how to assign each jet in the event to the corresponding
// parton-level jet.
// The number of permutations up to 10 jets are given below for Tt and
// TtH events.
//
// NJet         Npermutation (Tt)       Npermutation (TtH)
// 4            24                      --
// 5            120                     --
// 6            360                     360
// 7            840                     2520
// 8            1680                    10080
// 9            3024                    30240
// 10           5040                    75600
//
// The formulas for the number of permutations for Tt and TtH events
// given n jets in the event are
//
//         n!
// Tt:  -------- ; n >= 4
//      (n - 4)!
//
//          n!
// TtH: ---------- ; n >= 6
//      (n - 6)!2!
//
// The current MAX settings are chosen for a maximum number of 8 jets
// Increasing this limit should be done with caution, as it will
// increase the number of permutations rapidly.
//

namespace hitfit{

using std::vector;
using std::string;

/**
  @brief Template class of experiment-independent interface to HitFit.
  This class is intended to be used inside the programming environment of
  a specific experiment, where each type of physics objects has its
  own class/type.  For using HitFit with generic four-vector classes,
  user can't use this class and have to use the Top_Fit class directly.
  The reason is: this class is designed assuming electron and muon are
  represented by different object type, a situation which is guaranteed to
  happen in any experiments.
  The class contains some static integer constants to limit the maximum
  amount of jets in an event before fitting. See the description of those
  constants for details.  The numbers of permutations for \f$t\bar{t}\f$ and
  \f$t\bar{t}H\f$ as a function of the number of jets
  \f$N_{\mathrm{jet}}\f$ in the event for a few values of are

  <center>
  <table border="2">

  <tr>
  <th>\f$N_{\mathrm{jet}}\f$</th>
  <th>\f$N_{t\bar{t}}\f$</th>
  <th>\f$N_{t\bar{t}H}\f$</th>
  </tr>

  <tr>
  <td align="right">4</td>
  <td align="right">24</td>
  <td align="right">N/A</td>
  </tr>

  <tr>
  <td align="right">5</td>
  <td align="right">120</td>
  <td align="right">N/A</td>
  </tr>

  <tr>
  <td align="right">6</td>
  <td align="right">360</td>
  <td align="right">360</td>
  </tr>

  <tr>
  <td align="right">7</td>
  <td align="right">840</td>
  <td align="right">2520</td>
  </tr>

  <tr>
  <td align="right">8</td>
  <td align="right">1680</td>
  <td align="right">20160</td>
  </tr>

  </table>
  </center>

  If adjusting the limits defined by the static constants is desired, then
  please the following formulas.

  The number for possible
  permutations, \f$N_{t\bar{t}}\f$, as a function of number of jets,
  \f$n\f$, for \f$t\bar{t}\f$ event is given by:
  \f[
  N_{t\bar{t}}(n) = \frac{n!}{(n-4)!};~ n \ge 4
  \f]
  The number for possible permutations, \f$N_{t\bar{t}H}\f$, as a
  function of number of jets, \f$n\f$, for \f$t\bar{t}H\f$ is given by:
  \f[
  N_{t\bar{t}}(n) = \frac{n!}{(n-6)!2!};~ n \ge 6
  \f]

  @param AElectron The typename of the electron physics object class to
  be translated into HitFit's Lepjets_Event_Lep.

  @param AMuon The typename of the muon physics object class to
  be translated into HitFit's Lepjets_Event_Lep.

  @param AJet The typename of the jet physics object class to
  be translated into HitFit's Lepjets_Event_Jet.

  @param AMet The typename of the missing transverse energy physics
  object class be translated into HitFit's Fourvec.

 */
template <class AElectron,
          class AMuon,
          class AJet,
          class AMet>
class RunHitFit {
private:
  /**
     The translator from AElectron to Lepjets_Event_Lep.
   */
  LeptonTranslatorBase<AElectron> _ElectronTranslator;

  /**
     The translator from AMuon to Lepjets_Event_Lep.
   */
  LeptonTranslatorBase<AMuon>     _MuonTranslator;

  /**
     The translator from AJet to Lepjets_Event_Jet.
   */
  JetTranslatorBase<AJet>         _JetTranslator;

  /**
     The translator from AMet to Fourvec.
   */
  METTranslatorBase<AMet>         _METTranslator;

  /**
     The internal event, which only contains lepton and missing transverse energy.
   */
  Lepjets_Event                   _event;

  /**
     The internal array of jets.
     Jets are kept in this array and not added into the internal
     event. The reason is: the jet energy correction applied to a jet
     is dependent on the assumed jet type (b or light) in the
     permutation.  Therefore the decision is to store
     jets in their original format/data type.

     Before a fit to a particular permutation is done,
     this class convert the jets in this array into
     Lepjets_Event_Jet format, taking into consideration
     the assumed jet type and applying the appropriate jet energy
     correction.
   */
  vector<AJet>                    _jets;
  /**
     Indicator for b jets
   */
  vector<bool>                    _isBJet;
  /**
     Indicator for light quark jets
   */
  vector<bool>                    _isLJet;

  /**
     The interface between the event and the fitting algorithm.
   */
  Top_Fit                         _Top_Fit;

  /**
     The results of the kinematic fit.
   */
  vector<Fit_Result>              _Fit_Results;

public:
  /**
     @brief Constructor.

     @param el The function object to translate from AElectron to
     Lepjets_Event_Lep.

     @param mu The function object to translate from AMuon to
     Lepjets_Event_Lep.

     @param jet The function object to translate from AJet to
     Lepjets_Event_Jet.

     @param met The function object to translate from AMet to
     Fourvec.

     @param default_file The path of ASCII text files which contains the
     parameter settings for this instance of RunHitFit.

     @param lepw_mass The mass to which the leptonic \f$ W- \f$ boson should be
     constrained to.  A value of zero means this constraint will be removed.

     @param hadw_mass The mass to which the hadronic \f$ W- \f$ boson should be
     constrained to.  A value of zero means this constraint will be removed.

     @param top_mass The mass to which the top quark should be constrained
     to.  A value of zero means this constraint will be removed.
   */
  RunHitFit(const LeptonTranslatorBase<AElectron>& el,
            const LeptonTranslatorBase<AMuon>&     mu,
            const JetTranslatorBase<AJet>&         jet,
            const METTranslatorBase<AMet>&         met,
            const string                           default_file,
            double                                 lepw_mass,
            double                                 hadw_mass,
            double                                 top_mass):
      _ElectronTranslator(el),
      _MuonTranslator    (mu),
      _JetTranslator     (jet),
      _METTranslator     (met),
      _event             (0,0),
      _Top_Fit(Top_Fit_Args(Defaults_Text(default_file)),lepw_mass,hadw_mass,top_mass)
  {}

  /**
     @brief Destructor.
   */
  ~RunHitFit()
  {
  }

  /**
     @brief Clear the internal event, fit results, and jets.
   */
  void
  clear()
  {
    _event = Lepjets_Event(0,0);
    _jets.clear();
    _isBJet.clear();
    _isLJet.clear();
    _Fit_Results.clear();
  }

  /**
     @brief Add one electron into the internal event.

     @param electron The electron to be added into the internal event.
   */
  void
  AddLepton(const AElectron& electron)
  {
    _event.add_lep(_ElectronTranslator(electron,electron_label));
  }

  /**
     @brief Add one muon into the internal event.

     @param muon The muon to be added into the internal event.
  */
  void
  AddLepton(const AMuon& muon)
  {
    _event.add_lep(_MuonTranslator(muon,muon_label));
  }

  /**
     @brief Add one jet into the internal event.  This function will
     do nothing if the internal event has already contained the maximally
     allowed number of jets.

     Explanation about this function: This function does not directly add
     the jet into the internal event.  Rather, this function store the
     jet in an internal array.
     The reason is: jet energy correction
     and resolution depends on the jet type in the permutation.
     Therefore RunHitFit will only add jet into the event after a specific
     jet permutation has been determined.
     This is done in the FitAllPermutation function().

     @param jet The jet to be added into the internal event.
  */
  void
  AddJet(const AJet& jet, bool isB, bool isL)
  {
    if (_jets.size() < MAX_HITFIT_JET) {
      _jets.emplace_back(jet);
      _isBJet.push_back(isB);
      _isLJet.push_back(isL);
      //if (fabs(jet.eta())>2.4) std::cout << fabs(jet.eta()) << " " << isB << " " << isL << std::endl;
    }
  }

  /**
     @brief Set the missing transverse energy of the internal event.
   */
  void
  SetMet(const AMet& met)
  {
    _event.met()    = _METTranslator(met);
    _event.kt_res() = _METTranslator.KtResolution(met);
    _event.sumEt()  = met.sumEt();
  }

  /**
     @brief Set the \f$k_{T}\f$ resolution of the internal event.

     @param res The resolution.
   */
  void
  SetKtResolution(const Resolution& res)
  {
    _event.kt_res() = res;
  }

  /**
     @brief Set the \f$E_{T}\!\!\!\!/\f$ resolution of the internal event.

     @param res The \f$E_{T}\!\!\!\!/\f$ resolution, same as \f$k_{T}\f$
     resolution.
   */
  void
  SetMETResolution(const Resolution& res)
  {
    SetKtResolution(res);
  }

  /**
     @brief Return a constant reference to the underlying Top_Fit object.
   */
  const Top_Fit&
  GetTopFit() const
  {
    return _Top_Fit;
  }

  /**
     @brief Fit all permutations of the internal event.  Returns the
     number of permutations.
   */
  vector<Fit_Result>::size_type
  FitAllPermutation()
  {
    if (_jets.size() < MIN_HITFIT_JET) {
      // For ttbar lepton+jets, a minimum of MIN_HITFIT_JETS jets
      // is required
      return 0;
    }

    if (_jets.size() > MAX_HITFIT_JET) {
      // Restrict the maximum number of jets in the fit
      // to prevent loop overflow
      return 0;
    }

    _Fit_Results.clear();

    // Prepare the array of jet types for permutation
    // Both light quark jets are initially given hadw1_label, so that next_permutation behaves nicely.
    vector<int> jet_types (_jets.size(), unknown_label);
    jet_types[0] = lepb_label;
    jet_types[1] = hadb_label;
    jet_types[2] = hadw1_label;
    jet_types[3] = hadw1_label;

    if (_Top_Fit.args().do_higgs_flag() and _jets.size() >= MIN_HITFIT_TTH) {
      jet_types[4] = higgs_label;
      jet_types[5] = higgs_label;
    }

    std::stable_sort(jet_types.begin(),jet_types.end());

    // begin loop over unique jet permutations
    do {
      // Don't fit if the jet types don't match => saves much time.
      bool typeAgreement = true;
      // Look exactly at the four leading jets.
      for (size_t ijet = 0; ijet < 4; ++ijet) {
        int jT = jet_types[ijet];
        if (jT==hadw1_label) {
          if (!_isLJet[ijet]) {
            typeAgreement = false;
            break;
          }
        } else if (jT==hadb_label or jT==lepb_label) {
          if (!_isBJet[ijet]) {
            typeAgreement = false;
            break;
          }
        }
      }
      if (typeAgreement) {
        // Variables to be calculated only once per two neutrino solutions.
        double umwhad = 0, umthad = 0, nuz_store = 0;

        // loop over two neutrino solution
        for (int nusol = 0 ; nusol < 2 ; ++nusol) {
          // Copy the event
          Lepjets_Event fev = _event;

          // Add jets into the event, with the assumed type in accord with the permutation.
          // The translator _JetTranslator will correctly return object of Lepjets_Event_Jet with
          // jet energy correction applied in accord with the assumed jet type (b or light).
          for (size_t j = 0; j != _jets.size(); ++j) fev.add_jet(_JetTranslator(_jets[j],jet_types[j]));

          // Set jet types.
          fev.set_jet_types(jet_types);

          // Prepare the placeholders for unique kinematic quantities.
          double umtlep = 0;
          double mt, sigmt;
          Column_Vector pullx, pully;

          // Do the fit
          double chisq = _Top_Fit.fit_one_perm(fev,
                                               umwhad, umthad, umtlep,
                                               nuz_store,
                                               mt, sigmt,
                                               pullx, pully);
          // Store output of the fit
          _Fit_Results.emplace_back(chisq,
                                    fev,
                                    pullx, pully,
                                    umwhad,
                                    (umthad + umtlep)/2.,
                                    mt, sigmt);
        } // end loop over two neutrino solution
      }
    } while (std::next_permutation (jet_types.begin(), jet_types.end()));
    // end loop over unique jet permutations

    return _Fit_Results.size();
  }

  /**
      @brief Return the results of fitting all permutations of the
      internal event.
   */
  vector<Fit_Result>
  GetFitAllPermutation()
  {
    return _Fit_Results;
  }

  /**
     Minimum number of jet as input to HitFit in Tt event
   */
  static const unsigned int MIN_HITFIT_JET =   4 ;

  /**
     Minimum number of jet as input to HitFit in TtH event
   */
  static const unsigned int MIN_HITFIT_TTH =   6 ;

  /**
     Maximum number of jet as input to HitFit in each event
   */
  static const unsigned int MAX_HITFIT_JET =   8 ;

  /**
     Maximum number of HitFit permutation in each event.
   */
  static const unsigned int MAX_HITFIT     = 1680;

  /**
     Maximum number of fitted variables in HitFit in each event
   */
  static const unsigned int MAX_HITFIT_VAR =  32 ;
};

} // namespace hitfit

#endif // #ifndef RUNHITFIT_H
