#include <fstream>
#include "other_functions.h"


void update_times(grupo &G, trabajadores *family, double time){
  for(unsigned int i=0; i<G.size(); i++){family[G[i]].tstate += time;}
}


void update_massive(grupo &G, trabajadores *family, double time){
  for(unsigned int i=0; i<G.size(); i++){family[G[i]].time += time;}
}


void update_times_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time){
  update_times(Val[2], altos, time);	update_times(Vba[2], bajos, time); //Expuestos
  update_times(Val[3], altos, time);	update_times(Vba[3], bajos, time); //Expuestos testeados
  update_times(Val[4], altos, time);	update_times(Vba[4], bajos, time); //Expuestos aislados
  update_times(Val[5], altos, time);	update_times(Vba[5], bajos, time); //Presintomáticos
  update_times(Val[6], altos, time);	update_times(Vba[6], bajos, time); //Presintomáticos testeados
  update_times(Val[7], altos, time);	update_times(Vba[7], bajos, time); //Presintomáticos aislados
  update_times(Val[8], altos, time);	update_times(Vba[8], bajos, time); //Leves
  update_times(Val[9], altos, time);	update_times(Vba[9], bajos, time); //Leves testeados
  update_times(Val[10], altos, time);	update_times(Vba[10], bajos, time); //Leves aislados
  update_times(Val[11], altos, time);	update_times(Vba[11], bajos, time); //Graves
}


void update_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time){
  update_massive(Val[1], altos, time);	update_massive(Vba[1], bajos, time); //Susceptibles testeados
  update_massive(Val[3], altos, time);	update_massive(Vba[3], bajos, time); //Expuestos testeados
  update_massive(Val[13], altos, time);	update_massive(Vba[13], bajos, time); //Recuperados testeados
}


void print_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, double time, std::string name){
  std::ofstream fout;
  fout.open(name, std::ios_base::app);

  fout << time;
  for(unsigned int i=0; i<Val.size(); i++){fout << '\t' << Val[i].size() << '\t' << Vba[i].size();}
  fout << std::endl;

  fout.close();
}


void print_types(std::vector<grupo> &Val, std::vector<grupo> &Vba, double time, std::string name){
  unsigned int veca[6] = {}, vecb[6] = {};

  veca[4] += Val[11].size();    vecb[4] += Vba[11].size(); //Graves
  for(unsigned int i=0; i<=1; i++){veca[0] += Val[i].size();    vecb[0] += Vba[i].size();} //Susceptibles
  for(unsigned int i=0; i<=2; i++){
    veca[1] += Val[i+2].size();    vecb[1] += Vba[i+2].size(); //Expuestos
    veca[2] += Val[i+5].size();    vecb[2] += Vba[i+5].size(); //Presintomáticos
    veca[3] += Val[i+8].size();    vecb[3] += Vba[i+8].size(); //Leves
    veca[5] += Val[i+12].size();    vecb[5] += Vba[i+12].size(); //Recuperados
  }

  std::ofstream fout;
  fout.open(name, std::ios_base::app);

  fout << time;
  for(unsigned int i=0; i<6; i++){fout << '\t' << veca[i] << '\t' << vecb[i];}
  fout << std::endl;

  fout.close();
}


void print_inf(std::vector<grupo> &Val, std::vector<grupo> &Vba, double time, std::string name){
  unsigned int veca[3] = {}, vecb[3] = {};

  veca[2] += Val[11].size();    vecb[2] += Vba[11].size(); //Graves
  for(unsigned int i=0; i<=2; i++){
    veca[0] += Val[i+5].size();    vecb[0] += Vba[i+5].size(); //Presintomáticos
    veca[1] += Val[i+8].size();    vecb[1] += Vba[i+8].size(); //Leves
  }

  double value1, value2;
  value1 = veca[0] + veca[1] + veca[2];
  value2 = vecb[0] + vecb[1] + vecb[2];
  
  std::ofstream fout;
  fout.open(name, std::ios_base::app);
  fout << time << '\t' << value1 << '\t' << value2 << '\t' << value1 + value2 << std::endl;
  fout.close();
}


void print_net(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, std::string name){
  std::ofstream fout;
  fout.open(name);

  fout << "Source\tTarget" << std::endl;
  for(unsigned int i=12; i<=14; i++){
    for(unsigned int j=0; j<Val[i].size(); j++){
      for(unsigned int k=0; k<altos[Val[i][j]].my_inf.size(); k++){
	fout << Val[i][j] << '\t' << altos[Val[i][j]].my_inf[k] << std::endl;
      }
    }
    for(unsigned int j=0; j<Vba[i].size(); j++){
      for(unsigned int k=0; k<bajos[Vba[i][j]].my_inf.size(); k++){
	fout << Vba[i][j] + Na << '\t' << bajos[Vba[i][j]].my_inf[k] << std::endl;
      }
    }
  }

  fout.close();
}
