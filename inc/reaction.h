/**
 * @file reaction.h
 * @author Daniel Felipe
 * @date 2020
 * @brief Header containing the mother reaction and reactions functions
 */


#ifndef REACTION_H
#define REACTION_H

#include <Random64.h>
#include <bases.h>
#include <workers.h>

/* Esta función es la función plantilla para hacer las reacciones */
void mother_reaction(grupo &Out, grupo &In, int index, Workers *family, int typeout, int typein);


/* Con esta función hallo el índice de la persona que voy a hacer reaccionar */
int index_time(grupo &Out, Workers *family, lognormal_d &dist, double value);


/* Esta función me escoje el infectado, si el infeccioso es un externo */
int infected(Workers *family, int max);


/* Esta función me hace la reacción de contagio para un susceptible alto */
void reaction0(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de contagio para un susceptible bajo */
void reaction1(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de tránsito de un expuesto alto */
void reaction2(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de tránsito de un expuesto bajo */
void reaction3(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a leve alto */
void reaction4(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a leve bajo */
void reaction5(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a infeccioso grave aislado alto */
void reaction6(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a infeccioso grave aislado bajo */
void reaction7(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a recuperado alto*/
void reaction8(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a recuperado bajo*/
void reaction9(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de leve a recuperado alto */
void reaction10(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de leve a recuperado bajo */
void reaction11(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de infeccioso grave a recuperado alto */
void reaction12(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


/* Esta función me hace la reacción de infeccioso grave a recuperado bajo */
void reaction13(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);


#endif /* REACTION_H */
