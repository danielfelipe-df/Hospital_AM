/**
 * @file constants.h
 * @author Daniel Felipe
 * @date 2020
 * @brief Header containing the Constants class
 */


#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

/**
 * @brief This class contains all important constants variables that can be
 * modified by the external CSV file.
 *
 * @class Constants
 */

class Constants{
public:

  /* Attributes */

  /// Variable keeps the name of CSV file.
  std::string name;

  /// @f$(\xi)@f$ Test sensitivity
  double xi;
  /// @f$(\theta)@f$ Massive test covering
  double theta;
  /// @f$(\iota)@f$ Slight symptomatics ratio tested continually
  double iota;

  /// N95 mask effectiveness
  double N95;
  /// Surgical mask effectiveness
  double TBQ;
  /// Hands washing effectiveness
  double HW;
  /// Physical distance effectiveness
  double SDP;

  /// @f$(\alpha)@f$ Adherence to isolation
  double alpha;
  /// @f$(\mu)@f$ Low - High risk staff crossed interaction parameter
  double mu;
  /// @f$(\chi)@f$ Low risk staff internal interaction parameter
  double chi;
  /// @f$(\phi_1)@f$ High risk staff internal interaction parameter
  double phi1;
  /// @f$(\eta)@f$ Interaction parameter between patients and High risk staff
  double eta;
  /// @f$(\lambda)@f$ Interaction parameter between familiars and staff
  double lambda;

  /// Number of agents (staff) traced given a positive test result
  unsigned int trace;
  /// Number of agents (staff) in the positive test result's contact network
  unsigned int trace_net;

  /// @f$(\nu)@f$ Number of days with massive testing strategy
  unsigned int nu;
  /// @f$(\delta)@f$ Number of days without massive testing strategy
  unsigned int delta;

  /// Number of gaussian functions used in the external prevalence
  size_t N_gauss;
  /// Height of the curves' peak
  double *A_gauss;
  /// Position of the center of curves' peak
  double *Mu_gauss;
  /// Gaussian RMS width
  double *Sigma_gauss;

  /// Number of linear sections in the laboral @f$\beta@f$
  size_t N_betaL;
  /// Laboral @f$\beta@f$ slope
  double *m_betaL;
  /// Laboral @f$\beta@f$ y-intercept
  double *b_betaL;
  /// Time limits for the laboral @f$\beta@f$ linear section
  double *lim_betaL;

  /// Number of linear sections in the familiar @f$\beta@f$
  size_t N_betaF;
  /// Familiar @f$\beta@f$ slope
  double *m_betaF;
  /// Familiar @f$\beta@f$ y-intercept
  double *b_betaF;
  /// Time limits for the familiar @f$\beta@f$ linear section
  double *lim_betaF;

  /// Boolean for isolation of Slight symptomatics
  bool AisLev;


  /* Constructor */
  Constants();

  /* Destructor */
  ~Constants();

};

#endif /* CONSTANTS_H */
