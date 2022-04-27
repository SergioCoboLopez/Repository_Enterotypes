//25/01/2018. Este codigo proporciona un ejemplo muy simple de como se utilizan los unordered multimaps
//(los diccionarios de python) en cpp.

#include <iostream>
#include <fstream>
#include <unordered_map>
#include "armadillo"
#include <tuple>

using namespace std;
using namespace arma;

//Definimos una estructura para poder colocar dos elementos como value del multimap
//======================
struct PersonMicrobe
{int Person, Microbe;};
//======================



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main()
{
  //Cargamos Datos
  //============================================
  mat Data;
  Data.load("../Input_Data/micro_dataset2.txt");
  //============================================
  
  //Declaramos el diccionario
  //============================================
  std::unordered_multimap<float, PersonMicrobe > u;
  //============================================

  //Rellenamos el diccionario en un bucle usando "emplace"
  //============================================
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      u.emplace(Data(i,j),PersonMicrobe{i,j});
    }}
  //============================================

  //Imprimir por pantalla los contenidos del diccionario
  //==========================================================================
  for (unordered_multimap<float, PersonMicrobe>::iterator Values = u.begin();
     Values != u.end(); ++Values)
  {
    cout << "Microbe " << (*Values).first << " corresponds to person " <<
           (*Values).second.Person << " and to microbe " << (*Values).second.Microbe << endl;
  }
  //==========================================================================
  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}
