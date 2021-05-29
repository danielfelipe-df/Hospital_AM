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


/* Esta función es la función plantilla para hacer los cambios de estado normal a testeado */
void tested_reaction(grupo &Out, grupo &In, int index, Workers *family, int typeout, int typein, double ran, bool dist);


/* Esta función complementa la reacción cuando algún leve tipo iota se aisla. Aquí se le hace el test para saber si se queda ahí o se devuelve. */
void tested_lev_ais(int agent, Workers *family, double value, bool dist);


/* Esta función es la que genera la reacción en el testeo masivo */
void massive_reaction(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos);


/* Esta función me actualiza los tiempos de los testeados y si cumplieron el tiempo los aisla o los devuelvo al estado normal.*/
int tested_isolated_inf(grupo &T, grupo &TA, grupo &G, Workers *family, double time, int typeout, int typein1, int typein2, Crandom &ran);


/* Esta función me actualiza el tiempo de los testeados masivamente, pero fuera de la zona de testeo masivo. Si ya cumplieron, los muevo. */
void tested_massive(grupo &T, grupo &G, Workers *family, double time, int typeout, int typein);


/* Con esta función muevo los testeados masivos a su respectivo grupo después de pasar los días de testeo */
void move_massive(grupo &T, grupo &G, Workers *family, unsigned int typeout, unsigned int typein);


/* Con esta función actualizo el tiempo de los leves y reviso si el test es positivo o negativo */
void result_lev_ais(std::vector<grupo> &Val, std::vector<grupo> &Vba, Workers *altos, Workers *bajos, double time, Crandom &ran);


/* Con esta función muevo todos los testeados masivos a su respectivo grupo después de pasar los días de testeo */
void move_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, Workers *altos, Workers *bajos);


/* Esta función me actualiza el tiempo de los testeados masivamente, pero fuera de la zona de testeo masivo. Si ya cumplieron, los muevo. */
void tested_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, Workers *altos, Workers *bajos, double time);

#endif /* TEST_H */
