#include <iostream>
#include <fstream>
using namespace std;

//0. FUNCTIONS
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//0.1. Open File
//======================================================================
double Open_File(string& path,string& name,int& Number_of_rows,int& Number_of_columns)
{
  string file= path + name;
  cout<< file + '\n';
  ifstream Output( file.c_str());
  if (!Output) {
    cout << "Unable to open your file. Probably a typo.\n";
  }

  double Data_Matrix [Number_of_rows][Number_of_columns];
  int row_counter;int col_counter;
  for (row_counter = 0; row_counter < Number_of_rows; row_counter++) {
    for (col_counter = 0; col_counter < Number_of_columns; col_counter++) {
      Output >> Data_Matrix[row_counter][col_counter];
      //      cout << Data_Matrix[row_counter][col_counter];
    }
  }
  Output.close();
  return Data_Matrix;
  
}
//======================================================================








//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main()
{
  //Declarar variables
  //======================================================================
  int N_rows=1220; int N_cols=370;
  double Data_Matrix [N_rows][N_cols];
  int row_counter;
  int col_counter;
  
  std::string path_to_file= "/export/home/shared/Projects/Microbiome/Input_Data/";
  std::string name_of_file= "toy_data.txt";
  std::string file=path_to_file+name_of_file;
  //ifstream Output( file.c_str());
  //======================================================================

  //  ifstream Output;
  double Matrix [N_rows][N_cols];
  Matrix=Open_File(path_to_file,name_of_file,N_rows,N_cols);

  //Open_File(path_to_file,name_of_file);

  // for (row_counter = 0; row_counter < N_rows; row_counter++) {
  //   for (col_counter = 0; col_counter < N_cols; col_counter++) {
  //     Output >> Data_Matrix[row_counter][col_counter];
  //     cout << Data_Matrix[row_counter][col_counter];
  //   }
  // }
  // cout << '\n';
  // cout << Data_Matrix;
  // cout << '\n';

  // Output.close();
  
  return 0;
}
