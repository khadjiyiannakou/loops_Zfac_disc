#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <complex>


#define NT 48
#define NS 16
#define KAPPA 0.1373
#define MU 0.006
#define NSTOCH 7000

using namespace std;

// ========================= Matrix Class ==========================//
class Matrix{
private:
  int m;
  int n;
public:
  Matrix(int,int);
  ~Matrix();
  complex<double> **data;
  int size() const { return m*n;}
  int M() const {return m;}
  int N() const {return n;}
  void print();
  //  complex<double>& operator()(int,int);
  //  complex<double> operator()(int,int);
  Matrix operator*(Matrix&);
  Matrix operator-(Matrix&);
  Matrix operator*(complex<double>);
  Matrix commutator(Matrix&);
  complex<double> trace();
};

Matrix::Matrix(int M, int N):m(M),n(M){
  try{
  data = new complex<double>*[m];
  for(int i = 0; i < m; i++)
    data[i] = new complex<double>[n];
  }
  catch (bad_alloc& ba){
    wcerr << "Error caught:" << ba.what() << endl;
    exit(EXIT_FAILURE);
  }
  for(int i = 0 ; i < m ; i++)
    memset(data[i],0,n*sizeof(complex<double>));
}

Matrix::~Matrix(){
  for(int i = 0; i < n; i++)
    delete[] data[i];
  delete[] data;
}

void Matrix::print(){
  for(int i = 0 ; i < m ; i++){
    for(int j = 0 ; j < n ; j++)
      cout << scientific << showpos << " (" << data[i][j].real() << "," << data[i][j].imag() << ") ";
    cout << endl;
  }
}

Matrix Matrix::operator *(Matrix &B){
  if(n != B.M()){
    wcerr << "Error matrix dimensions do not agree\n";
    exit(EXIT_FAILURE);
  }
  Matrix C(m,B.N());

  for(int i = 0 ; i < m ; i++)
    for(int j = 0 ; j < B.N() ; j++)
      for(int k = 0 ; k < n ; k++)
	C.data[i][j] += data[i][k] * B.data[k][j];

  return C;
}

Matrix Matrix::operator *(complex<double> a){
  Matrix C(m,n);
  for(int i = 0 ; i < m ; i++)
    for(int j = 0 ; j < n ; j++)
      C.data[i][j] = data[i][j] * a;
  return C;
}

Matrix Matrix::operator -(Matrix &B){
  Matrix C(m,n);
  for(int i = 0 ; i < m ; i++)
    for(int j = 0 ; j < n ; j++)
      C.data[i][j] = data[i][j] - B.data[i][j];
  return C;
}

Matrix Matrix::commutator(Matrix &B){
  Matrix L = (*this) * B;
  Matrix R = B * (*this);
  return L-R;
}
//complex<double>& Matrix::operator()(int i ,int j){
//  return data[i][j];
//}
//============================================================//

complex<double> Matrix::trace(){
  if(m!=n){
    wcerr << "Error, only square matrices have trace\n";
    exit(EXIT_FAILURE);
  }
  complex<double> res = complex<double>(0,0);
  for(int i = 0 ; i < m ; i++)
    res += data[i][i];
  return res;
}

int main(int argc, char* argv[]){

  if(argc != 4){
    wcerr << "Error wrong number of inputs provided\n" ;
    exit(EXIT_FAILURE);
  }
  //======================== Open files ==========================//
  fstream f_standard, f_generalized, f_output;
  f_standard.open(argv[1], fstream::in);
  if(!f_standard){
    wcerr << "Error cannot find loops for standard end trick\n";
    exit(EXIT_FAILURE);
  }

  f_generalized.open(argv[2], fstream::in);
  if(!f_generalized){
    wcerr << "Error cannot find loops for generalized end trick\n";
    exit(EXIT_FAILURE);
  }

  f_output.open(argv[3], fstream::out);
  if(!f_output){
    wcerr << "Error cannot open file for writting results\n";
    exit(EXIT_FAILURE);
  }
  //=========================== Read data from files==================================//
  int dummy;
  complex<double> *loop_std_raw = new complex<double>[NT*NS];
  complex<double> *loop_gen_raw = new complex<double>[NT*NS];
  if(!loop_std_raw || !loop_gen_raw){
    wcerr << "Error allocating memory\n";
    exit(EXIT_FAILURE);
  }

  cout.precision(7);
  for(int i = 0 ; i < NT*NS; i++){
    f_standard >> dummy >> dummy >> dummy >> dummy >> dummy >> loop_std_raw[i].real() >> loop_std_raw[i].imag();
    f_generalized >> dummy >> dummy >> dummy >> dummy >> dummy >> loop_gen_raw[i].real() >> loop_gen_raw[i].imag();
    //    cout << scientific << showpos << loop_gen_raw[i].real() << "\t" << loop_gen_raw[i].imag() << endl;
  }
  cout << "Data loaded successfully\n";
  //=============================================================//
  //  Matrix g(4,4);
  Matrix g1(4,4);
  Matrix g2(4,4);
  Matrix g3(4,4);
  Matrix g4(4,4);
  Matrix one(4,4);
  complex<double> imagUnit = complex<double>(0,1);

  g1.data[0][3].imag() = 1;
  g1.data[1][2].imag() = 1;
  g1.data[2][1].imag() = -1;
  g1.data[3][0].imag() = -1;

  g2.data[0][3].real() = 1;
  g2.data[1][2].real() = -1;
  g2.data[2][1].real() = -1;
  g2.data[3][0].real() = 1;

  g3.data[0][2].imag() = 1;
  g3.data[1][3].imag() = -1;
  g3.data[2][0].imag() = -1;
  g3.data[3][1].imag() = 1;

  g4.data[0][0].real() = 1;
  g4.data[1][1].real() = 1;
  g4.data[2][2].real() = -1;
  g4.data[3][3].real() = -1;
  
  one.data[0][0].real() = 1;
  one.data[1][1].real() = 1;
  one.data[2][2].real() = 1;
  one.data[3][3].real() = 1;
  

  Matrix g5 = g1*g2*g3*g4;
  Matrix ig5 = g5*imagUnit;
  Matrix I = one*imagUnit;
  Matrix g5g1 = g5*g1;
  Matrix g5g2 = g5*g2;
  Matrix g5g3 = g5*g3;
  Matrix g5g4 = g5*g4;
  Matrix Comm_g4g1 = g4.commutator(g1);
  Matrix Comm_g4g2 = g4.commutator(g2);
  Matrix Comm_g4g3 = g4.commutator(g3);
  Matrix Comm_g1g2 = g1.commutator(g2); 
  Matrix Comm_g1g3 = g1.commutator(g3); 
  Matrix Comm_g2g3 = g2.commutator(g3); 

  Matrix ig5Comm_g4g1 = g5*Comm_g4g1*imagUnit*0.5;
  Matrix ig5Comm_g4g2 = g5*Comm_g4g2*imagUnit*0.5;
  Matrix ig5Comm_g4g3 = g5*Comm_g4g3*imagUnit*0.5;

  Matrix ig5Comm_g1g2 = g5*Comm_g1g2*imagUnit*0.5;
  Matrix ig5Comm_g1g3 = g5*Comm_g1g3*imagUnit*0.5;
  Matrix ig5Comm_g2g3 = g5*Comm_g2g3*imagUnit*0.5;

  Matrix **m_std = new Matrix *[NT];
  for(int i = 0 ; i < NT ; i++)
    m_std[i] = new Matrix(4,4);

  Matrix **m_gen = new Matrix *[NT];
  for(int i = 0 ; i < NT ; i++)
    m_gen[i] = new Matrix(4,4);

  for(int i = 0 ; i < NT ; i++)
    for(int a = 0 ; a < 4 ; a++)
      for(int b = 0 ; b < 4 ; b++){
	m_std[i]->data[b][a] = loop_std_raw[i*16+a*4+b]; // take the transpose 
	m_gen[i]->data[b][a] = loop_gen_raw[i*16+a*4+b];
      }
  //  m_std[0]->print();

  //========================== create loops =====================//
  double results[16];
  complex<double> mul_std;
  complex<double> mul_gen;
  mul_std.real() = 0.;
  mul_std.imag() = (8.*MU*KAPPA*KAPPA)/NSTOCH;
  
  mul_gen.real() = (4.*KAPPA)/NSTOCH;
  mul_gen.imag() = 0.;

  for(int i = 0 ; i < 16 ; i++){
    results[i]=0.;
  }
  // Example) scalar "Needs ig5 in twisted basis, needs standard one end trick, mul factor i*8*mu*kappa^2, expected from real the part, 1 minus sign for the standard one end trick included in the code, 1 minus sign from anticommuting \bar{q}q, no sign from quda->tmLQCD "
  //Scalar -> real, +1
  //Pseudoscalar -> real, -1
  //vector -> imaginary, +1
  //Axial vector -> real, -1
  //tensor -> imaginary, +1

  for(int i = 0 ; i < NT ; i++){
    results[0] += real(((ig5*(*m_std[i])).trace())*mul_std); // scalar
    results[1] += real(((I*(*m_std[i])).trace())*mul_std*(-1.)); //pseudoscalar
    results[2] += imag(((g4*(*m_gen[i])).trace())*mul_gen); // vector g0
    results[3] += imag(((g1*(*m_gen[i])).trace())*mul_gen); // vector g1
    results[4] += imag(((g2*(*m_gen[i])).trace())*mul_gen); // vector g2
    results[5] += imag(((g3*(*m_gen[i])).trace())*mul_gen); // vector g3
    results[6] += real(((g5g4*(*m_gen[i])).trace())*mul_gen*(-1.)); // axial vector g5g0
    results[7] += real(((g5g1*(*m_gen[i])).trace())*mul_gen*(-1.)); // axial vector g5g1
    results[8] += real(((g5g2*(*m_gen[i])).trace())*mul_gen*(-1.)); // axial vector g5g2
    results[9] += real(((g5g3*(*m_gen[i])).trace())*mul_gen*(-1.)); // axial vector g5g3

    results[10] += imag(((ig5Comm_g4g1*(*m_std[i])).trace())*mul_std); // tensor [g0,g1]
    results[11] += imag(((ig5Comm_g4g2*(*m_std[i])).trace())*mul_std); // tensor [g0,g2]
    results[12] += imag(((ig5Comm_g1g2*(*m_std[i])).trace())*mul_std); // tensor [g1,g2]
    results[13] += imag(((ig5Comm_g4g3*(*m_std[i])).trace())*mul_std); // tensor [g0,g3]
    results[14] += imag(((ig5Comm_g1g3*(*m_std[i])).trace())*mul_std); // tensor [g1,g3]
    results[15] += imag(((ig5Comm_g2g3*(*m_std[i])).trace())*mul_std); // tensor [g2,g3]

  }

  // multiply minus sign coming from interchanging fermion fields
  for(int i = 0 ; i < 16 ; i++)
    results[i] *= (-1.);

  // dump results
  // it seems that tensor is zero and we do not print it
  for(int i = 0 ; i < 16; i++)
    cout << scientific << showpos << results[i] << endl;
  //========================== clear memory =====================//



  delete[] loop_std_raw;
  delete[] loop_gen_raw;

  for(int i = 0; i < NT; i++)
    delete m_std[i];
  delete[] m_std;
  
  for(int i = 0; i < NT; i++)
    delete m_gen[i];
  delete[] m_gen;
  
  return 0;
}

/*
  ig5.print();
  cout << endl;
  I.print();
  cout << endl;
  g1.print();
  cout << endl;
  g2.print();
  cout << endl;
  g3.print();
  cout << endl;
  g4.print();
  cout << endl;
  g5g1.print();
  cout << endl;
  g5g2.print();
  cout << endl;
  g5g3.print();
  cout << endl;
  g5g4.print();
  cout << endl;
  ig5Comm_g4g1.print();
  cout << endl;
  ig5Comm_g4g2.print();
  cout << endl;
  ig5Comm_g4g3.print();
  cout << endl;
  ig5Comm_g1g2.print();
  cout << endl;
  ig5Comm_g1g3.print();
  cout << endl;
  ig5Comm_g2g3.print();
  cout << endl;
  ig5Comm_g4g1.print();
  cout << endl;
  ig5Comm_g4g2.print();
  cout << endl;
  ig5Comm_g4g3.print();
  cout << endl;
  ig5Comm_g1g2.print();
  cout << endl;
  ig5Comm_g1g3.print();
  cout << endl;
  ig5Comm_g2g3.print();
  cout << endl;

*/
