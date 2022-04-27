#include <iostream>
#include "armadillo"

using namespace std;
using namespace arma;

struct Matrices{
  mat M1;
  int a;
};

void Fun_Matrix(mat *A,int i){
  (*A)(0,1)+=i;
  return;
    }

void Fun_Double(double *n,int i){
  *n=i;
  cout<<n<<endl;
  cout<<*n<<endl;
  return;
}

void Fun_Struct(Matrices *Matrices1){
  cout<<"Primer Elemento con * "<<endl;
  cout<<(*Matrices1).M1<<endl;
  cout<<"\n Primer Elemento con -> "<<endl;
  cout<<Matrices1->M1<<endl;

  cout<<"\n \n Segundo Elemento con * "<<endl;
  cout<<(*Matrices1).a<<endl;
  cout<<"\n Segundo Elemento con -> "<<endl;
  cout<<Matrices1->a<<endl;
return;
}

int main () {

  int filas=4;int columnas=4;
  mat A=zeros(filas,columnas);
  double n=5.0;
  double n1=4.0;

  Matrices Matrices1;
  Matrices1.M1=ones(filas,columnas);
  Matrices1.a=4;

  Fun_Struct(&Matrices1);

  // for (int i=0;i<10;i++){
    // Fun_Matrix(&A,i);
    // cout<<"Valor de la matriz para la iteracion "<<i<<endl;
    // cout<<A<<endl;

    // Fun_Double(&n,i);
    // cout<<"Valor del numero para la iteracion "<<i<<endl;
    // cout<<n<<endl;

    
  //}


  

  return 0;
}
