#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <vector>
#include <Random64.h>
#include <bases.h>
#include <Constants/cons_dynamics.h>
#include <Constants/cons_contact.h>
#include <trabajadores.h>


/* Esta función me dice cuánto tiempo se demora en hacerse una reacción y qué reacción es. */
std::vector<double> contagio(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, double t, double* tj);


/* Esta función implementa el método de bisección para hallar el tiempo */
double biseccion(double* A, double t, double B, double ranr, double* tj, int n, std::vector<lognormal_d> &dist);


/* Es la función integrada de la prevalencia externa */
double function(double A, double prom, double sigma, double t, double tau);


/* Es la función phi */
double phi(double* A, double* tj, unsigned int n, double B, double deltat, double t, std::vector<lognormal_d> &dist);


/* Es la función auxiliar de phi para cuando los límites de la integral no se diferencian por más de 1.0 de diferencia */
void aux_phi_function(double t0, double diff, double &value_gnum, double &value_bnum);


/* Es la función auxiliar principal de phi para hacer toda la integración */
void main_aux_phi_function(double t0, double diff, double &value_gnum, double &value_bnum, double t);


/* Con esta función identifico a la persona que hico la infección */
int who_infected(grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, double cons1, double cons2, Crandom &ran, int alti, double t, double TBa, double TBb);


/* Esta es la función madre para la selección del que infectó */
int selection_infectious(grupo &Ga, grupo &Gb, grupo &Gc, grupo &Gd, Crandom &ran);


/* Esta es la función gaussiana */
double function_gauss(double x, double A, double prom, double sigma);


/* Esta es la integral, entre t0 y t1, de una gaussiana por una ecuación lineal */
double int_beta_gauss(double t0, double t1, double prom, double sigma, double A, double m, double b);


/* Esta es la integral de una ecuación lineal */
double int_beta(double t0, double t1, double m, double b);


/* Esta es la función para el beta */
double function_beta(double x);


/* Esta es la función que me devuelve el índice de la función en beta a usar */
size_t index_beta(double x);


#endif
