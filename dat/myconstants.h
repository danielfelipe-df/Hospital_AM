#ifndef MYCONSTANTS_H
#define MYCONSTANTS_H

/* Prevalencias externas *******************************/

//Prevalencia externa alta
const size_t N_gaussA = 2;
const double A_gaussA[N_gaussA] = {1.52e-3, 1.39e-2};
const double Mu_gaussA[N_gaussA] = {8.51e1, 1.27e2};
const double Sigma_gaussA[N_gaussA] = {8.25, 3.52e1};

//Prevalencia externa media
const size_t N_gaussM = 2;
const double A_gaussM[N_gaussM] = {8.04e-3, 1.11e-3};//, 1.18e-3, 4.4e-3};
const double Mu_gaussM[N_gaussM] = {1.07e2, 9.23e1};//, 1.85e2, 2.32e2};
const double Sigma_gaussM[N_gaussM] = {3.45e1, 6.13};//, 1.35e1, 4.46e1};

//Prevalencia externa baja
const size_t N_gaussB = 2;
const double A_gaussB[N_gaussB] = {1.33e-3, 6.87e-4};
const double Mu_gaussB[N_gaussB] = {1.41e2, 7.84e1};
const double Sigma_gaussB[N_gaussB] = {5.5e1, 4.15e1};

/******************************************************/


/* Tasas de contagio **********************************/

// Tasa de contagio laboral
const size_t N_betaL = 2;
const double m_betaL[N_betaL] = {0.0, 0.0};
const double b_betaL[N_betaL] = {0.88, 0.0};
const double lim_betaL[N_betaL+1] = {0.0, 1/3.0, 1.0};

// Tasa de contagio familiar
const size_t N_betaF = 2;
const double m_betaF[N_betaF] = {0.0, 0.0};
const double b_betaF[N_betaF] = {0.0, 0.88};
const double lim_betaF[N_betaF+1] = {0.0, 1/3.0, 1.0};

/******************************************************/

#endif
