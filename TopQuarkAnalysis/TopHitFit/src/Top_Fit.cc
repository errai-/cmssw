// File: src/Top_Fit.cc
// Purpose: Handle jet permutations.
// Created: Jul, 2000, sss, based on run 1 mass analysis code.
//
// XXX handle merging jets.
// XXX btagging for ttH.
//
// CMSSW File      : src/Top_Fit.cc
// Original Author : Scott Stuart Snyder <snyder@bnl.gov> for D0
// Imported to CMSSW by Haryo Sumowidagdo <Suharyo.Sumowidagdo@cern.ch>
//


/**
    @file Top_Fit.cc

    @brief Handle and fit jet permutations of an event.  This is the
    primary interface between user's Lepjets_Event and HitFit kinematic
    fitting algorithm.  See the documentation for the header file
    Top_Fit.h for details.

    @author Scott Stuart Snyder <snyder@bnl.gov>

    @par Creation date:
    Jul 2000.

    @par Modification History:
    Apr 2009: Haryo Sumowidagdo <Suharyo.Sumowidagdo@cern.ch>:
    Imported to CMSSW.<br>
    Nov 2009: Haryo Sumowidagdo <Suharyo.Sumowidagdo@cern.ch>:
    Added doxygen tags for automatic generation of documentation.

    @par Terms of Usage:
    With consent for the original author (Scott Snyder).

 */

#include "TopQuarkAnalysis/TopHitFit/interface/Top_Fit.h"
#include "TopQuarkAnalysis/TopHitFit/interface/Lepjets_Event.h"
#include "TopQuarkAnalysis/TopHitFit/interface/Top_Decaykin.h"
#include "TopQuarkAnalysis/TopHitFit/interface/Defaults.h"
#include "TopQuarkAnalysis/TopHitFit/interface/Fit_Results.h"
#include "TopQuarkAnalysis/TopHitFit/interface/fourvec.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>

using std::cout;
using std::endl;
using std::abs;
using std::next_permutation;
using std::stable_sort;
using std::vector;
using std::ostream;


namespace hitfit {

//*************************************************************************
// Argument handling.
//


Top_Fit_Args::Top_Fit_Args (const Defaults& defs)
//
// Purpose: Constructor.
//
// Inputs:
//   defs -        The Defaults instance from which to initialize.
//
  : _print_event_flag (defs.get_bool ("print_event_flag")),
    _do_higgs_flag (defs.get_bool ("do_higgs_flag")),
    _jet_mass_cut (defs.get_float ("jet_mass_cut")),
    _mwhad_min_cut (defs.get_float ("mwhad_min_cut")),
    _mwhad_max_cut (defs.get_float ("mwhad_max_cut")),
    _mtdiff_max_cut (defs.get_float ("mtdiff_max_cut")),
    _nkeep (defs.get_int ("nkeep")),
    _solve_nu_tmass (defs.get_bool ("solve_nu_tmass")),
    _args (defs)
   {
}


bool Top_Fit_Args::print_event_flag () const
//
// Purpose: Return the print_event_flag parameter.
//          See the header for documentation.
//
{
  return _print_event_flag;
}


bool Top_Fit_Args::do_higgs_flag () const
//
// Purpose: Return the do_higgs_flag parameter.
//          See the header for documentation.
//
{
  return _do_higgs_flag;
}


double Top_Fit_Args::jet_mass_cut () const
//
// Purpose: Return the jet_mass_cut parameter.
//          See the header for documentation.
//
{
  return _jet_mass_cut;
}


double Top_Fit_Args::mwhad_min_cut () const
//
// Purpose: Return the mwhad_min_cut parameter.
//          See the header for documentation.
//
{
  return _mwhad_min_cut;
}


double Top_Fit_Args::mwhad_max_cut () const
//
// Purpose: Return the mwhad_max_cut parameter.
//          See the header for documentation.
//
{
  return _mwhad_max_cut;
}


double Top_Fit_Args::mtdiff_max_cut () const
//
// Purpose: Return the mtdiff_max_cut parameter.
//          See the header for documentation.
//
{
  return _mtdiff_max_cut;
}


int Top_Fit_Args::nkeep () const
//
// Purpose: Return the nkeep parameter.
//          See the header for documentation.
//
{
  return _nkeep;
}


bool Top_Fit_Args::solve_nu_tmass () const
//
// Purpose: Return the solve_nu_tmass parameter
//          See the header for documentation.
//
{
  return _solve_nu_tmass;
}


const Constrained_Top_Args& Top_Fit_Args::constrainer_args () const
//
// Purpose: Return the contained subobject parameters.
//
{
  return _args;
}



//*************************************************************************
// Helper functions.
//

namespace {

/**
    @brief Helper function: apply mass cuts to see if this
    event should be rejected before fitting.

    @param ev The event to test.

    @param args The parameter settings.

    @param mwhad  The hadronic  \f$ W- \f$ boson mass.

    @param umthad The mass of the hadronic top quark before fit.

    @param umtlep The mass of the leptonic top quark before fit.
 */
bool test_for_bad_masses (const Lepjets_Event& ev,
                          const Top_Fit_Args& args,
                          double mwhad,
                          double umthad,
                          double umtlep)
//
// Purpose: Apply mass cuts to see if this event should be rejected
//          without fitting.
//
// Inputs:
//   ev -          The event to test.
//   args -        Parameter setting.
//   mwhad -       The hadronic W mass.
//   umthad -      The hadronic top mass.
//   umtlep -      The leptonic top mass.
//
// Returns:
//   True if the event should be rejected.
//
{
  // Reject the event if any jet's mass is too large.
  if (ev.sum (lepb_label).m()  > args.jet_mass_cut() ||
      ev.sum (hadb_label).m()  > args.jet_mass_cut() ||
      ev.sum (hadw1_label).m() > args.jet_mass_cut() ||
      ev.sum (hadw2_label).m() > args.jet_mass_cut())
      return true;

  // Reject if if the hadronic W mass is outside the window.
  if (mwhad < args.mwhad_min_cut()) return true;

  // Reject if if the hadronic W mass is outside the window.
  if (mwhad > args.mwhad_max_cut()) return true;

  // And if the two top masses are too far apart.
  if (abs (umthad - umtlep) > args.mtdiff_max_cut()) return true;

  // It's ok.
  return false;
}

} // unnamed namespace


//*************************************************************************


Top_Fit::Top_Fit (const Top_Fit_Args& args,
                  double lepw_mass,
                  double hadw_mass,
                  double top_mass)
//
// Purpose: Constructor.
//
// Inputs:
//   args -        The parameter settings for this instance.
//   lepw_mass -   The mass to which the leptonic W should be constrained,
//                 or 0 to skip this constraint.
//   hadw_mass -   The mass to which the hadronic W should be constrained,
//                 or 0 to skip this constraint.
//   top_mass -    The mass to which the top quarks should be constrained,
//                 or 0 to skip this constraint.
//
  : _args (args),
    _constrainer(args.constrainer_args(), lepw_mass, hadw_mass, top_mass),
    _lepw_mass(lepw_mass),
    _hadw_mass(hadw_mass)
{
}


double Top_Fit::fit_one_perm (Lepjets_Event& ev,
                              double& umwhad,
                              double& umthad,
                              double& umtlep,
                              double& nuz_store,
                              double& mt,
                              double& sigmt,
                              Column_Vector& pullx,
                              Column_Vector& pully)
//
// Purpose: Fit a single jet permutation.
//
// Inputs:
//   ev -          The event to fit.
//                 The object labels must have already been assigned.
//
// Outputs:
//   ev-           The event after the fit.
//   umwhad -      Hadronic W mass before fitting.
//   umthad -      Top mass before fitting, hadronic.
//   umtlep -      Top mass before fitting, leptonic.
//   nuz_store -   Store for the second neutrino solution.
//   mt -          Top mass after fitting.
//   sigmt -       Top mass uncertainty after fitting.
//   pullx -       Vector of pull quantities for well-measured variables.
//   pully -       Vector of pull quantities for poorly-measured variables.
//
// Returns:
//   The fit chisq, or < 0 if the fit didn't converge.
//
// Adaptation note by Haryo Sumowidagdo:
//   This function is rewritten in order to make its purpose reflects
//   the function's name.  The function nows only fit one jet permutation
//   with one neutrino solution only.
//
//
{
  mt = 0;
  sigmt = 0;

  double nuz = 0;
  // The hadronic results and the neutrino solutions need to be calculated only once per two neutrino solutions.
  if (umwhad==0 and umthad==0 and nuz_store==0) {
    umwhad = Top_Decaykin::hadw(ev).m();
    umthad = Top_Decaykin::hadt(ev).m();

    // Find the neutrino solutions by requiring either:
    // 1) that the leptonic top have the same mass as the hadronic top.
    // 2) that the mass of the lepton and neutrino is equal to the W mass
    if (_args.solve_nu_tmass())
      Top_Decaykin::solve_nu_tmass (ev, umthad, nuz, nuz_store);
    else
      Top_Decaykin::solve_nu (ev, _lepw_mass, nuz, nuz_store);
  } else {
    // The neutrino solutions have already been found; use the stored value.
    nuz = nuz_store;
  }

  // Set up to use the selected neutrino solution
  ev.met().setZ(nuz);

  // Note: We have set the neutrino Pz, but we haven't set the neutrino energy.
  // Remember that originally the neutrino energy was equal to
  // sqrt(nu_px*nu_px + nu_py*nu_py).  Calculating the invariant mass squared
  // for the neutrino will give negative mass squared.
  // Therefore we need to adjust (increase) the neutrino energy in order to
  // make its mass remain zero.
  adjust_e_for_mass(ev.met(),0);

  // The leptonic top mass can be calculated only after the neutrino has been set.
  umtlep = Top_Decaykin::lept(ev).m();

  // Trace, if requested.
  if (_args.print_event_flag()) {
    cout << "Top_Fit::fit_one_perm() : Before fit:\n";
    Top_Decaykin::dump_ev (cout, ev);
  }

  // Maybe reject this event.
  umwhad = Top_Decaykin::hadw(ev).m();

  if (_hadw_mass > 0 and test_for_bad_masses(ev, _args, umwhad, umthad, umtlep)) {
    cout << "Top_Fit: bad mass comb.\n";
    return -999;
  }

  // Do the fit.
  double chisq = _constrainer.constrain(ev, mt, sigmt, pullx, pully);

  // Trace, if requested.
  if (_args.print_event_flag()) {
    cout << "Top_Fit::fit_one_perm() : After fit:\n";
    cout << "chisq: " << chisq << " mt: " << mt << " ";
    Top_Decaykin::dump_ev(cout, ev);
  }

  // Done!
  return chisq;
}

/**
  @brief Output stream operator, print the content of this Top_Fit object
  to an output stream.

  @param s The output stream to which to write.

  @param fitter The instance of Top_Fit to be printed.
 */
std::ostream& operator<< (std::ostream& s, const Top_Fit& fitter)
//
// Purpose: Print the object to S.
//
// Inputs:
//   s -           The stream to which to write.
//   fitter -      The object to write.
//
// Returns:
//   The stream S.
//
{
  return s << fitter._constrainer;
}

const Top_Fit_Args& Top_Fit::args() const
{
  return _args;
}

} // namespace hitfit
