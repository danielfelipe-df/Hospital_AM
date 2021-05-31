#include <fstream>
#include <other_functions.h>


void update_times(group &G, Workers *family, double time){
  for(unsigned int i=0; i<G.size(); i++){family[G[i]].tstate += time;}
}


void update_massive(group &G, Workers *family, double time){
  for(unsigned int i=0; i<G.size(); i++){family[G[i]].time += time;}
}


void update_times_all(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Workers *altos, Workers *bajos, double time){
  update_times(Val["EXP"], altos, time);	update_times(Vba["EXP"], bajos, time); //Expuestos
  update_times(Val["EXPT"], altos, time);	update_times(Vba["EXPT"], bajos, time); //Expuestos testeados
  update_times(Val["EXPA"], altos, time);	update_times(Vba["EXPA"], bajos, time); //Expuestos aislados
  update_times(Val["PRE"], altos, time);	update_times(Vba["PRE"], bajos, time); //Presintomáticos
  update_times(Val["PRET"], altos, time);	update_times(Vba["PRET"], bajos, time); //Presintomáticos testeados
  update_times(Val["PREA"], altos, time);	update_times(Vba["PREA"], bajos, time); //Presintomáticos aislados
  update_times(Val["MSYM"], altos, time);	update_times(Vba["MSYM"], bajos, time); //Leves
  update_times(Val["MSYMT"], altos, time);	update_times(Vba["MSYMT"], bajos, time); //Leves testeados
  update_times(Val["MSYMA"], altos, time);	update_times(Vba["MSYMA"], bajos, time); //Leves aislados
  update_times(Val["SSYMA"], altos, time);	update_times(Vba["SSYMA"], bajos, time); //Graves
}


void update_massive_all(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Workers *altos, Workers *bajos, double time){
  update_massive(Val["SUST"], altos, time);	update_massive(Vba["SUST"], bajos, time); //Susceptibles testeados
  update_massive(Val["EXPT"], altos, time);	update_massive(Vba["EXPT"], bajos, time); //Expuestos testeados
  update_massive(Val["RECT"], altos, time);	update_massive(Vba["RECT"], bajos, time); //Recuperados testeados
}


void print_all(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, double time, std::string name){
  std::ofstream fout;
  fout.open(name, std::ios_base::app);

  fout << time << '\t' ;
  fout << Val["SUS"].size() << '\t' << Vba["SUS"].size() << '\t';
  fout << Val["SUST"].size() << '\t' << Vba["SUST"].size() << '\t';
  fout << Val["SUSA"].size() << '\t' << Vba["SUSA"].size() << '\t';
  fout << Val["EXP"].size() << '\t' << Vba["EXP"].size() << '\t';
  fout << Val["EXPT"].size() << '\t' << Vba["EXPT"].size() << '\t';
  fout << Val["EXPA"].size() << '\t' << Vba["EXPA"].size() << '\t';
  fout << Val["PRE"].size() << '\t' << Vba["PRE"].size() << '\t';
  fout << Val["PRET"].size() << '\t' << Vba["PRET"].size() << '\t';
  fout << Val["PREA"].size() << '\t' << Vba["PREA"].size() << '\t';
  fout << Val["MSYM"].size() << '\t' << Vba["MSYM"].size() << '\t';
  fout << Val["MSYMT"].size() << '\t' << Vba["MSYMT"].size() << '\t';
  fout << Val["MSYMA"].size() << '\t' << Vba["MSYMA"].size() << '\t';
  fout << Val["SSYMA"].size() << '\t' << Vba["SSYMA"].size() << '\t';
  fout << Val["RECI"].size() << '\t' << Vba["RECI"].size() << '\t';
  fout << Val["RECT"].size() << '\t' << Vba["RECT"].size() << '\t';
  fout << Val["RECA"].size() << '\t' << Vba["RECA"].size() << '\t';
  fout << std::endl;

  fout.close();
}


void print_types(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, double time, std::string name){
  unsigned int veca[6] = {}, vecb[6] = {};

  // Susceptibles
  veca[0] += Val["SUS"].size();  vecb[0] += Vba["SUS"].size();
  veca[0] += Val["SUST"].size();  vecb[0] += Vba["SUST"].size();
  veca[0] += Val["SUSA"].size();  vecb[0] += Vba["SUSA"].size();

  // Expuestos
  veca[1] += Val["EXP"].size();  vecb[1] += Vba["EXP"].size();
  veca[1] += Val["EXPT"].size();  vecb[1] += Vba["EXPT"].size();
  veca[1] += Val["EXPA"].size();  vecb[1] += Vba["EXPA"].size();

  // Pre-sintomáticos
  veca[2] += Val["PRE"].size();  vecb[2] += Vba["PRE"].size();
  veca[2] += Val["PRET"].size();  vecb[2] += Vba["PRET"].size();
  veca[2] += Val["PREA"].size();  vecb[2] += Vba["PREA"].size();

  // Sintomáticos leves
  veca[3] += Val["MSYM"].size();  vecb[3] += Vba["MSYM"].size();
  veca[3] += Val["MSYMT"].size();  vecb[3] += Vba["MSYMT"].size();
  veca[3] += Val["MSYMA"].size();  vecb[3] += Vba["MSYMA"].size();

  //Graves
  veca[4] += Val["SSYMA"].size();    vecb[4] += Vba["SSYMA"].size();

  // Recuperados
  veca[5] += Val["RECI"].size();  vecb[5] += Vba["RECI"].size();
  veca[5] += Val["RECT"].size();  vecb[5] += Vba["RECT"].size();
  veca[5] += Val["RECA"].size();  vecb[5] += Vba["RECA"].size();


  std::ofstream fout;
  fout.open(name, std::ios_base::app);

  fout << time;
  for(unsigned int i=0; i<6; i++){fout << '\t' << veca[i] << '\t' << vecb[i];}
  fout << std::endl;

  fout.close();
}


void print_inf(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, double time, std::string name){
  unsigned int veca[3] = {}, vecb[3] = {};

  // Pre-sintomáticos
  veca[0] += Val["PRE"].size();  vecb[0] += Vba["PRE"].size();
  veca[0] += Val["PRET"].size();  vecb[0] += Vba["PRET"].size();
  veca[0] += Val["PREA"].size();  vecb[0] += Vba["PREA"].size();

  // Sintomáticos leves
  veca[1] += Val["MSYM"].size();  vecb[1] += Vba["MSYM"].size();
  veca[1] += Val["MSYMT"].size();  vecb[1] += Vba["MSYMT"].size();
  veca[1] += Val["MSYMA"].size();  vecb[1] += Vba["MSYMA"].size();

  //Graves
  veca[2] += Val["SSYMA"].size();    vecb[2] += Vba["SSYMA"].size();

  double value1, value2;
  value1 = veca[0] + veca[1] + veca[2];
  value2 = vecb[0] + vecb[1] + vecb[2];

  std::ofstream fout;
  fout.open(name, std::ios_base::app);
  fout << time << '\t' << value1 << '\t' << value2 << '\t' << value1 + value2 << std::endl;
  fout.close();
}


void print_net(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Workers *altos, Workers *bajos, std::string name){
  std::string names[3] = {"RECI", "RECT", "RECA"};
  std::ofstream fout;
  fout.open(name);

  fout << "Source\tTarget" << std::endl;
  for(unsigned int i=0; i<3; i++){
    for(unsigned int j=0; j<Val[names[i]].size(); j++){
      for(unsigned int k=0; k<altos[Val[names[i]][j]].my_inf.size(); k++){
	fout << Val[names[i]][j] << '\t' << altos[Val[names[i]][j]].my_inf[k] << std::endl;
      }
    }
    for(unsigned int j=0; j<Vba[names[i]].size(); j++){
      for(unsigned int k=0; k<bajos[Vba[names[i]][j]].my_inf.size(); k++){
	fout << Vba[names[i]][j] + Na << '\t' << bajos[Vba[names[i]][j]].my_inf[k] << std::endl;
      }
    }
  }

  fout.close();
}
