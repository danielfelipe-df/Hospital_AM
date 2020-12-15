#ifndef CONS_DYNAMICS_H
#define CONS_DYNAMICS_H

/* Periodos de los comportamentos ***************************************************/

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


/* Variables de la tasa de transmisión ***************************************************/

const size_t N_beta = 3; //Número de funciones en las que se divide beta
const double m_beta[N_beta] = {9.396, -4.698, 0.0}; //Pendiente de la función que define el beta
const double b_beta[N_beta] = {0.001, 3.5245, 0.001}; //Corte en y de la función que define el beta
const double lim_beta[N_beta+1] = {0.0, 0.25, 0.75, 1.0}; //Límites de los tramos

const double N95 = 0.31; //Reducción en la probabilidad de contagiarse usando tapabocas N95
const double TBQ = 0.84; //Reducción en el probabilidad de contagiarse usando tapabocas quirúrgico
const double HW = 0.64; //Reducción en la probabilidad de contagiarse si se lava las manos
const double SDP = 0.89; //Reducción en la probabilidad de contagiarse si mantiene distanciamiento


/* Proporciones de compartimentos infecciosos ******************************************************/

const double kappa = 0.2; //Proporción de contagiados que son asintomáticos
const double psi = 0.95; //Proporción de sintomáticos que son leves


/* Prevalencia externa ***********************************************/

/* Alta */
const size_t N_gauss = 2;
const double A_gauss[N_gauss] = {1.52e-3, 1.39e-2};
const double Mu_gauss[N_gauss] = {8.51e1, 1.27e2};
const double Sigma_gauss[N_gauss] = {8.25, 3.52e1};


/* Media
const size_t N_gauss = 4;
const double A_gauss[N_gauss] = {8.04e-3, 1.11e-3, 1.18e-3, 4.4e-3};
const double Mu_gauss[N_gauss] = {1.07e2, 9.23e1, 1.85e2, 2.32e2};
const double Sigma_gauss[N_gauss] = {3.45e1, 6.13, 1.35e1, 4.46e1};
*/

/* Baja
const size_t N_gauss = 2;
const double A_gauss[N_gauss] = {1.33e-3, 6.87e-4};
const double Mu_gauss[N_gauss] = {1.41e2, 7.84e1};
const double Sigma_gauss[N_gauss] = {5.5e1, 4.15e1};
*/

/**************************************************************************/

#endif
