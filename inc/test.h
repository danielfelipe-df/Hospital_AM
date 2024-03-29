/**
 * @file test.h
 * @author Daniel Felipe
 * @date 2020
 * @brief Header containing the testing functions and reactions
 */


#ifndef TEST_H
#define TEST_H

#include <Random64.h>
#include <bases.h>
#include <workers.h>

#include <map>


/* Esta función es la función plantilla para hacer los cambios de estado normal a testeado */
void tested_reaction(group &Out, group &In, int index, Workers *family, Stages typeout, Stages typein, double ran, bool dist);


/* Esta función complementa la reacción cuando algún leve tipo iota se aisla. Aquí se le hace el test para saber si se queda ahí o se devuelve. */
void tested_lev_ais(int agent, Workers *family, double value, bool dist);


/* Esta función es la que genera la reacción en el testeo masivo */
void massive_reaction(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos);


/* Esta función me actualiza los tiempos de los testeados y si cumplieron el tiempo los aisla o los devuelvo al estado normal.*/
int tested_isolated_inf(group &T, group &TA, group &G, Workers *family, double time, Stages typeout, Stages typein1, Stages typein2, Crandom &ran);


/* Esta función me actualiza el tiempo de los testeados masivamente, pero fuera de la zona de testeo masivo. Si ya cumplieron, los muevo. */
void tested_massive(group &T, group &G, Workers *family, double time, Stages typeout, Stages typein);


/* Con esta función muevo los testeados masivos a su respectivo group después de pasar los días de testeo */
void move_massive(group &T, group &G, Workers *family, Stages typeout, Stages typein);


/* Con esta función actualizo el tiempo de los leves y reviso si el test es positivo o negativo */
void result_lev_ais(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Workers *altos, Workers *bajos, double time, Crandom &ran);


/* Con esta función muevo todos los testeados masivos a su respectivo grupo después de pasar los días de testeo */
void move_massive_all(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Workers *altos, Workers *bajos);


/* Esta función me actualiza el tiempo de los testeados masivamente, pero fuera de la zona de testeo masivo. Si ya cumplieron, los muevo. */
void tested_massive_all(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Workers *altos, Workers *bajos, double time);

#endif /* TEST_H */
