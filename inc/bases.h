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
const double xi = 0.9; //Sensibilidad de la prueba
const double theta = 0.5; //Cobertura
const double iota = 0.6; //Proporción de sintomáticos leves testeados continuamente
const lognormal_d dist_Tt(0.51192, 0.41694); //Distribución del delay de entrega de resultados
const double TM = 0.5; //Tiempo mínimo que deben pasar los agentes en cada estado
const double N95 = 0.31; //Reducción en la probabilidad de contagiarse usando tapabocas N95
const double TBQ = 0.84; //Reducción en el probabilidad de contagiarse usando tapabocas quirúrgico
const double HW = 0.64; //Reducción en la probabilidad de contagiarse si se lava las manos
const double SDP = 0.89; //Reducción en la probabilidad de contagiarse si mantiene distanciamiento
const size_t N_beta = 3; //Número de funciones en las que se divide beta
const double m_beta[N_beta] = {9.396, -4.698, 0.0}; //Pendiente de la función que define el beta
const double b_beta[N_beta] = {0.001, 3.5245, 0.001}; //Corte en y de la función que define el beta
const double lim_beta[N_beta+1] = {0.0, 0.25, 0.75, 1.0}; //Límites de los tramos

/* Tasas de contacto ******************************************************/

const double alpha = 0.9; //Adherencia al aislamiento
const double mu = 0.1; //Tasa de contacto cruzada
const double chi = 0.7; //Tasa de contacto de tipo bajo
const double eta = 0.7; //Tasa de contacto con externos
const double phi1 = 0.6; //Tasa de contacto de tipo alto

/* Rastreo ****************************************************************/

const int trace = 3; //Número de personas rastreadas cada vez que haya contagio
const unsigned int trace_net = 6; //Número de contactos en promedio por cada persona del staff


/* Parámetros de gaussianas ***********************************************/

/* Prevalencia externa alta */
const size_t N_gauss = 2;
const double A_gauss[N_gauss] = {1.52e-3, 1.39e-2};
const double Mu_gauss[N_gauss] = {8.51e1, 1.27e2};
const double Sigma_gauss[N_gauss] = {8.25, 3.52e1};


/* Prevalencia externa media
const size_t N_gauss = 4;
const double A_gauss[N_gauss] = {8.04e-3, 1.11e-3, 1.18e-3, 4.4e-3};
const double Mu_gauss[N_gauss] = {1.07e2, 9.23e1, 1.85e2, 2.32e2};
const double Sigma_gauss[N_gauss] = {3.45e1, 6.13, 1.35e1, 4.46e1};
*/

/* Prevalencia externa baja
const size_t N_gauss = 2;
const double A_gauss[N_gauss] = {1.33e-3, 6.87e-4};
const double Mu_gauss[N_gauss] = {1.41e2, 7.84e1};
const double Sigma_gauss[N_gauss] = {5.5e1, 4.15e1};
*/

/**************************************************************************/

#endif
