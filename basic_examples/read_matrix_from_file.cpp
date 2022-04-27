#include <iostream>
#include <fstream>
using namespace std;

//Declarar variables
//======================================================================
double distances [4][3];  
ifstream in("/export/home/shared/Projects/Microbiome/Input_Data/micro_dataset.txt", ios::in);
//======================================================================

//Funcion para comprobar nombre del archivo de entrada
//======================================================================
void Alerta_Nombre (ifstream& in)
{
  if (!in) {
    cout << "Unable to open your file. Probably a typo.\n";
  }
}
//======================================================================

int main()
 {
  int fila;
  int columna;

  Alerta_Nombre (in);


  for (fila = 1; fila < 4; fila++) {
    for (columna = 0; columna < 3; columna++) {
      in >> distances[fila][columna];
      cout << distances[fila][columna];
      cout << ' ';
    }
  }
  cout << '\n';
  cout << distances[2][2];
  cout << '\n';

  in.close();
  return 0; 
 }
