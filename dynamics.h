#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <vector>
#include "Random64.h"
#include "bases.h"


/* Esta función me dice cuánto tiempo se demora en hacerse una reacción
 * y qué reacción es.
 */
std::vector<double> contagio(double Sa, double Sb, double Ea, double Eb, double Pa, double Pb, double PTa, double PTb, double La, double Lb, double LTa, double LTb, double Ia, double Ib, double ITa, double ITb, double Na, double Nb, double prev, Crandom &ran, double r, double *IT, double *FT, double t);


/* Esta función es la función plantilla para hacer las reacciones */
void mother_reaction(grupo &Out, grupo &In, Crandom &ran);


/* Esta función implementa el método de la bisección para hallar el tau de la reacción debido a la prevalencia externa.*/
double biseccion(double A, double prom, double sigma, double T, double S, double t, double B);


/* Esta función es la que se usa para la bisección y la evolución temporal*/
double function(double A, double prom, double sigma, double t, double tau);


/* Esta función me hace la reacción de contagio para un susceptible alto */
inline void reaction0(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Sa, Ea, ran);
}


/* Esta función me hace la reacción de contagio para un susceptible bajo */
inline void reaction1(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Sb, Eb, ran);
}


/* Esta función me hace la reacción de tránsito de un expuesto alto */
inline void reaction2(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Ea, Pa, ran);
}


/* Esta función me hace la reacción de tránsito de un expuesto bajo */
inline void reaction3(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Eb, Pb, ran);
}


/* Esta función me hace la reacción de tránsito de un pre-sintomático a leve alto */
inline void reaction4(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Pa, La, ran);
}


/* Esta función me hace la reacción de tránsito de un pre-sintomático a leve bajo */
inline void reaction5(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Pb, Lb, ran);
}


/* Esta función me hace la reacción de tránsito de un pre-sintomático a infeccioso grave aislado alto */
inline void reaction6(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Pa, IAa, ran);
}


/* Esta función me hace la reacción de tránsito de un pre-sintomático a infeccioso grave aislado bajo */
inline void reaction7(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Pb, IAb, ran);
}


/* Esta función me hace la reacción de tránsito de un pre-sintomático a pre-sintomático testeado y aislado alto */
inline void reaction8(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Pa, PTAa, ran);
}


/* Esta función me hace la reacción de tránsito de un pre-sintomático a pre-sintomático testeado y aislado bajo */
inline void reaction9(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Pb, PTAb, ran);
}


/* Esta función me hace la reacción de tránsito de un leve a aislado leve alto */
inline void reaction10(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(La, LAa, ran);
}


/* Esta función me hace la reacción de tránsito de un leve a aislado leve bajo */
inline void reaction11(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Lb, LAb, ran);
}


/* Esta función me hace la reacción de tránsito de un leve a leve aislado y testeado alto */
inline void reaction12(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(La, LTAa, ran);
}


/* Esta función me hace la reacción de tránsito de un leve a leve aislado y testeado bajo */
inline void reaction13(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Lb, LTAb, ran);
}


/* Esta función me hace la reacción de pre-sintomático a recuperado alto */
inline void reaction14(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Pa, Ra, ran);
}


/* Esta función me hace la reacción de pre-sintomático a recuperado bajo */
inline void reaction15(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Pb, Rb, ran);
}


/* Esta función me hace la reacción de leve a recuperado alto */
inline void reaction16(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(La, Ra, ran);
}


/* Esta función me hace la reacción de leve a recuperado bajo */
inline void reaction17(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(Lb, Rb, ran);
}


/* Esta función me hace la reacción de pre-sintomático testeado y aislado a recuperado alto */
inline void reaction18(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(PTAa, Ra, ran);
}


/* Esta función me hace la reacción de pre-sintomático testeado y aislado a recuperado bajo */
inline void reaction19(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(PTAb, Rb, ran);
}


/* Esta función me hace la reacción de leve aislado a recuperado alto */
inline void reaction20(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(LAa, Ra, ran);
}


/* Esta función me hace la reacción de leve aislado a recuperado bajo */
inline void reaction21(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(LAb, Rb, ran);
}


/* Esta función me hace la reacción de leve aislado y testeado a recuperado alto */
inline void reaction22(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(LTAa, Ra, ran);
}


/* Esta función me hace la reacción de leve aislado y testeado a recuperado bajo */
inline void reaction23(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(LTAb, Rb, ran);
}


/* Esta función me hace la reacción de tránsito de infeccioso grave a recuperado alto */
inline void reaction24(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(IAa, Ra, ran);
}


/* Esta función me hace la reacción de tránsito de infeccioso grave a recuperado bajo*/
inline void reaction25(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran){
  mother_reaction(IAb, Rb, ran);
}



#endif
