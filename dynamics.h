#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <vector>
#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/gamma.hpp>
#include "Random64.h"
#include "bases.h"
#include "trabajadores.h"


typedef boost::math::weibull_distribution<> weib_d;
typedef boost::math::gamma_distribution<> gamma_d;


/* Esta función me dice cuánto tiempo se demora en hacerse una reacción
 * y qué reacción es.
 */
std::vector<double> contagio(double Sa, double Sb, double STa, double STb, double Ea, double Eb, double ETa, double ETb, double Pa, double Pb, double PTa, double PTb, double PTAa, double PTAb, double La, double Lb, double LTa, double LTb, double LTAa, double LTAb, double IAa, double IAb, double Na, double Nb, double prev, Crandom &ran, double t, double* tj);


/* Esta función es la función plantilla para hacer las reacciones */
void mother_reaction(grupo &Out, grupo &In, int index, trabajadores *family, int typeout, int typein);


/* Esta función es la que genera la reacción en el testeo masivo */
void massive_reaction(grupo &S, grupo &ST, grupo &E, grupo &ET, grupo &P, grupo &PT, grupo &L, grupo &LT, grupo &RI, grupo &RT, Crandom &ran, trabajadores *family);


/* Esta función es la que genera la reacción en el testeo continuo */
void continue_reaction(grupo &L, grupo &LT, trabajadores *family, Crandom &ran);


/* Esta función me actualiza los tiempos de los testeados y los aisla si ya cumplieron el tiempo */
void tested_isolated(grupo &T, grupo &TA, grupo &G, trabajadores *family, double time, int typeout, int typein1, int typein2, bool infected, Crandom &ran);


/* Esta función implementa el método de bisección para hallar el tiempo */
double biseccion(double* A, double prom, double sigma, double t, double B, double ranr, double* tj, int n,  std::vector<weib_d> &dist);


/* Es la función integrada de la prevalencia externa */
double function(double A, double prom, double sigma, double t, double tau);


/* Es la función phi */
double phi(double* A, double* tj, unsigned int n, double prom, double sigma, double B, double deltat, double t,  std::vector<weib_d> &dist);


/* Con esta función hallo el índice de la persona que voy a hacer reaccionar */
int index_time(grupo &Out, trabajadores *family, weib_d &dist, double value);


/* Con esta función actualizo los tiempos de estado de ese grupo */
void update_times(grupo &G, trabajadores *family, double time);


/* Esta función me hace la reacción de contagio para un susceptible alto */
void reaction0(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de contagio para un susceptible bajo */
void reaction1(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de tránsito de un expuesto alto */
void reaction2(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de tránsito de un expuesto bajo */
void reaction3(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a leve alto */
void reaction4(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a leve bajo */
void reaction5(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a infeccioso grave aislado alto */
void reaction6(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a infeccioso grave aislado bajo */
void reaction7(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a recuperado alto*/
void reaction8(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a recuperado bajo*/
void reaction9(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de leve a recuperado alto */
void reaction10(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de leve a recuperado bajo */
void reaction11(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de infeccioso grave a recuperado alto */
void reaction12(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


/* Esta función me hace la reacción de infeccioso grave a recuperado bajo */
void reaction13(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);


#endif
