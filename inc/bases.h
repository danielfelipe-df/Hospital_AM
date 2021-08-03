/**
 * @file bases.h
 * @author Daniel Felipe
 * @date 2020
 * @brief File containing the Constants variable, aliases of different objects
 * and other constants aren't in Constants class.
 */


#ifndef BASES_H
#define BASES_H

/* Constants Class ****************************************************/

#include <constants.h>

/// Create the Constants variable
const Constants MyCons;


/* Agents container *****************************************************/

#include <vector>

/// Index agents container alias
typedef std::vector<unsigned int> group; 


/* Distributions and alias ************************************************/

#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/normal.hpp>

/// Weibull distribution alias
typedef boost::math::weibull_distribution<> weib_d;
/// Gamma distribution alias
typedef boost::math::gamma_distribution<> gamma_d;
/// Lognormal distribution alias
typedef boost::math::lognormal_distribution<> lognormal_d;
/// Normal distribution alias
typedef boost::math::normal_distribution<> normal_d;


/* Size of populations **********************************************/

/// Number of people in the system
const int N = 1000;
/// Number of high risk agents
const int Na = N*0.18;
/// Number of low risk agents
const int Nb = N*0.82;


/* Dynamic periods ***************************************************/

/// Latent period
const double De = 2.9;
/// Pre-symptomatic infectious period for slight agents
const double Dpl = 1.3;
/// Pre-symptomatic infectious period for severe agents
const double Dpg = 2.3;
/// Infectious period for slight agents
const double Dil = 1.7;
/// Infectious period for severe agents
const double Dig = 2.9;
/// Reciprocal of latent period
const double USDe = 1.0/De;
/// Reciprocal of pre-symptomatic infectious period for slight agents
const double USDpl = 1.0/Dpl;
/// Reciprocal of pre-symptomatic infectious period for severe agents
const double USDpg = 1.0/Dpg;
/// Reciprocal of infectious period for slight agents
const double USDil = 1.0/Dil;
/// Reciprocal of infectious period for slight agents
const double USDig = 1.0/Dig;
/// Reciprocal of complete infectious period
const double USDplil = 1.0/(Dpl+Dil);


/* Complementary ratios **************************************************/

/// Asymptomatic infected ratio
const double kappa = 0.2;
/// Slight symptomatics ratio
const double psi = 0.95;
/// Minimum time the agents must keep in each state
const double TM = 0.5;
/// Isolation time for traced agents
const double TtraceMax = 15;
/// Time delivery results distribution
const lognormal_d dist_Tt(0.51192, 0.41694);

/**************************************************************************/

#endif /* BASES_H */
