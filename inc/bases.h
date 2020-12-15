#ifndef BASES_H
#define BASES_H

/* Índice de personas *****************************************************/

#include <vector>

typedef std::vector<int> grupo; //Redefino el vector de índices de personas como grupo

/* Distribuciones *********************************************************/

#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/normal.hpp>

typedef boost::math::weibull_distribution<> weib_d;
typedef boost::math::gamma_distribution<> gamma_d;
typedef boost::math::lognormal_distribution<> lognormal_d;
typedef boost::math::normal_distribution<> normal_d;

/* Tamaño de las poblaciones **********************************************/

const int N = 1000; //Número de personas en el sistema
const int Na = N*0.18; //Número de personas de riesgo alto
const int Nb = N*0.82; //Número de personas de riesgo bajo


#endif
