#ifndef BASES_H
#define BASES_H

/* Clase de constantes ****************************************************/

#include <constants.h>
#include <string>

const constants MyCons;


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

/* Tasas de la dinámica ***************************************************/

const double De = 2.9; //Periodo latente
const double Dpl = 1.3; //Periodo infeccioso pre-sintomático para leves
const double Dpg = 2.3; //Periodo infecciosos pre-sintomático para graves
const double Dil = 1.7; //Periodo infeccioso de leves
const double Dig = 2.9; //Periodo infeccioso de graves
const double USDe = 1.0/De; //Inverso periodo latente
const double USDpl = 1.0/Dpl; //Inverso del periodo infeccioso pre-sintomático para leves
const double USDpg = 1.0/Dpg; //Inverso del periodo infeccioso pre-sintomático para graves
const double USDil = 1.0/Dil; //Inverso del periodo infeccioso de leves
const double USDig = 1.0/Dig; //Inverso del periodo infeccioso de graves
const double USDplil = 1.0/(Dpl+Dil); //Inverso del periodo infeccioso total


/* Tasas complementarias **************************************************/

const double kappa = 0.2; //Proporción de contagiados que son asintomáticos
const double psi = 0.95; //Proporción de sintomáticos que son leves
const double TM = 0.5; //Tiempo mínimo que deben pasar los agentes en cada estado
const double TtraceMax = 15; //Tiempo que los agentes pasaran aislados debido al rastreo
const lognormal_d dist_Tt(0.51192, 0.41694); //Distribución de los tiempos de entrega de examen

/**************************************************************************/

#endif
