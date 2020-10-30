#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <vector>
#include "Random64.h"
#include "bases.h"
#include "trabajadores.h"


/* Esta función me dice cuánto tiempo se demora en hacerse una reacción y qué reacción es. */
std::vector<double> contagio(std::vector<grupo> &Val, std::vector<grupo> &Vba, double prev, Crandom &ran, double t, double* tj);


/* Esta función implementa el método de bisección para hallar el tiempo */
double biseccion(double* A, double prom, double sigma, double t, double B, double ranr, double* tj, int n,  std::vector<lognormal_d> &dist);


/* Es la función integrada de la prevalencia externa */
double function(double A, double prom, double sigma, double t, double tau);


/* Es la función phi */
double phi(double* A, double* tj, unsigned int n, double prom, double sigma, double B, double deltat, double t,  std::vector<lognormal_d> &dist);


#endif
