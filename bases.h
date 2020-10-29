#ifndef BASES_H
#define BASES_H

/**************************************************************************/


#include <vector>

typedef std::vector<int> grupo; //Redefino el vector de índices de personas como grupo


/**************************************************************************/


#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/lognormal.hpp>

typedef boost::math::weibull_distribution<> weib_d;
typedef boost::math::gamma_distribution<> gamma_d;
typedef boost::math::lognormal_distribution<> lognormal_d;


/**************************************************************************/


const int N = 1000; //Número de personas en el sistema
const int Na = N*0.18; //Número de personas de riesgo alto
const int Nb = N*0.82; //Número de personas de riesgo bajo

const double beta = 1.3; //Beta de infección
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

const double kappa = 0.8; //Proporción de contagiados que son asintomáticos
const double psi = 0.95; //Proporción de sintomáticos que son leves
const double xi = 0.65; //Sensibilidad de la prueba
const double theta = 0.5; //Cobertura
const double iota = 0.6; //Proporción de sintomáticos leves testeados continuamente
const double Tt = 0.5; //Tiempo de reporte de la prueba
const double TM = 0.5; //Tiempo mínimo que deben pasar los agentes en cada estado

const double alpha = 0.95; //Adherencia al aislamiento
const double mu = 0.1; //Tasa de contacto cruzada
const double chi = 0.7; //Tasa de contacto de tipo bajo
const double eta = 0.7; //Tasa de contacto con externos
const double phi1 = 0.6; //Tasa de contacto de tipo alto

const double rho = Dpl + Dil - Tt; //Tiempo de aislamiento de asintomáticos
const double epsilon = Dig - Tt; //Tiempo de aislamiento de infecciosos graves
const double lambda = Dil - Tt; //Tiempo de aislamiento de infecciosos leves (continuo)
const double tau = lambda; //Tiempo de aislamiento de infecciosos leves (masivo)
const double USrho = 1.0/rho; //Inverso del tiempo de aislamiento de asintomáticos
const double USepsilon = 1.0/epsilon; //Inverso del tiempo de aislamiento de infecciosos graves
const double USlambda = 1.0/lambda; //Inverso del tiempo de aislamiento de infecciosos leves (continuo)
const double UStau = 1.0/tau; //Inverso del tiempo de aislamiento de infecciosos leves (masivo)

const int trace = 3; //Número de personas rastreadas cada vez que haya contagio

#endif
