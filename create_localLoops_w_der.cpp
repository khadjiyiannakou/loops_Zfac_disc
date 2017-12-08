#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <complex>
#include <vector>
#include <iomanip>
#include <stdio.h>

#define NT 48
#define NS 16
#define KAPPA 0.139305
#define MU 0.006
#define NSTOCH 512

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
  void operator=(Matrix&);
  void transpose();
  Matrix commutator(Matrix&);
  void signFlip();
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

void Matrix::operator =(Matrix &A){
  if( (this->M() != A.M()) || (this->N() != A.N()) ){
    wcerr << "Error matrix dimensions do not agree\n";
    exit(EXIT_FAILURE);
  }
  for(int i = 0 ; i < A.M(); i++)
    for(int j = 0 ; j < A.N(); j++)
      this->data[i][j] = A.data[i][j];
}

void Matrix::transpose(){
  Matrix tmp(m,n);
  tmp=*this;
  for(int i = 0 ; i < m ; i++)
    for(int j = 0 ; j < n ; j++)
      this->data[i][j] = tmp.data[j][i];
}

void Matrix::signFlip(){
  for(int i = 0 ; i < m ; i++)
    for(int j = 0 ; j < n ; j++)
      this->data[i][j] *= -1;
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

  if(argc != 6){
    wcerr << "Error wrong number of inputs provided\n" ;
    exit(EXIT_FAILURE);
  }
  //======================== Open files ==========================//
  fstream f_standard, f_generalized, f_output_ultra_local, f_output_vD, f_output_aD, f_output_tD;
  fstream f_standard_der, f_generalized_der;

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
  
  f_standard_der.open(argv[3], fstream::in);
  if(!f_standard_der){
    wcerr << "Error cannot find loops for standard derivative end trick\n";
    exit(EXIT_FAILURE);    
  }

  f_generalized_der.open(argv[4], fstream::in);
  if(!f_generalized_der){
    wcerr << "Error cannot find loops for generalized derivative end trick\n";
    exit(EXIT_FAILURE);    
  }

  string out_str(argv[5]);
  string ultra_local_str = out_str + "_ultra_local.dat";
  string vD_str = out_str + "_vD.dat";
  string aD_str = out_str + "_aD.dat";
  string tD_str = out_str + "_tD.dat";

  f_output_ultra_local.open(ultra_local_str.c_str(), fstream::out);
  if(!f_output_ultra_local){
    wcerr << "Error cannot open file for writting results\n";
    exit(EXIT_FAILURE);
  }

  f_output_vD.open(vD_str.c_str(), fstream::out);
  if(!f_output_vD){
    wcerr << "Error cannot open file for writting results\n";
    exit(EXIT_FAILURE);
  }

  f_output_aD.open(aD_str.c_str(), fstream::out);
  if(!f_output_aD){
    wcerr << "Error cannot open file for writting results\n";
    exit(EXIT_FAILURE);
  }

  f_output_tD.open(tD_str.c_str(), fstream::out);
  if(!f_output_tD){
    wcerr << "Error cannot open file for writting results\n";
    exit(EXIT_FAILURE);
  }


  //=========================== Read data from files==================================//
  int dummy;
  complex<double> *loop_std_raw = new complex<double>[NT*NS];
  complex<double> *loop_gen_raw = new complex<double>[NT*NS];
  complex<double> *loop_std_der_raw = new complex<double>[4*NT*NS];
  complex<double> *loop_gen_der_raw = new complex<double>[4*NT*NS];

  if(!loop_std_raw || !loop_gen_raw || !loop_std_der_raw || !loop_gen_der_raw){
    wcerr << "Error allocating memory\n";
    exit(EXIT_FAILURE);
  }

  cout.precision(7);
  for(int i = 0 ; i < NT*NS; i++){
    f_standard >> dummy >> dummy >> dummy >> dummy >> dummy >> loop_std_raw[i].real() >> loop_std_raw[i].imag();
    f_generalized >> dummy >> dummy >> dummy >> dummy >> dummy >> loop_gen_raw[i].real() >> loop_gen_raw[i].imag();
    //    cout << scientific << showpos << loop_gen_raw[i].real() << "\t" << loop_gen_raw[i].imag() << endl;
  }
  
  for(int mu = 0 ; mu < 4 ; mu++)
    for(int it = 0 ; it < NT; it++)
      for(int is = 0 ; is < NS ; is++){
	int i = mu*NT*NS + it*NS + is;
	f_standard_der >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> loop_std_der_raw[i].real() >> loop_std_der_raw[i].imag();
	f_generalized_der >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> loop_gen_der_raw[i].real() >> loop_gen_der_raw[i].imag();
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
	m_std[i]->data[b][a] = loop_std_raw[i*NS+a*4+b]; // take the transpose 
	m_gen[i]->data[b][a] = loop_gen_raw[i*NS+a*4+b];
      }
  //  m_std[0]->print();

  Matrix **m_std_der = new Matrix *[4*NT];
  Matrix **m_gen_der = new Matrix *[4*NT];

  for(int i = 0 ; i < 4*NT ; i++){
    m_std_der[i] = new Matrix(4,4);
    m_gen_der[i] = new Matrix(4,4);
  }

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int i = 0 ; i < NT ; i++)
      for(int a = 0 ; a < 4 ; a++)
	for(int b = 0 ; b < 4 ; b++){
	  m_std_der[mu*NT+i]->data[b][a] = loop_std_der_raw[mu*NT*NS + i*NS + a*4 + b]; // take the transpose
	  m_gen_der[mu*NT+i]->data[b][a] = loop_gen_der_raw[mu*NT*NS + i*NS + a*4 + b];
	}
  //========================== create loops =====================//
  double results[NS];
  complex<double> mul_std;
  complex<double> mul_gen;
  mul_std.real() = 0.;
  mul_std.imag() = (8.*MU*KAPPA*KAPPA)/NSTOCH;
  
  mul_gen.real() = (4.*KAPPA)/NSTOCH;
  mul_gen.imag() = 0.;

  for(int i = 0 ; i < NS ; i++){
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
    f_output_ultra_local << scientific << showpos << results[i] << endl;



  //////////// results for one-Der //////////////
  double results_vD[10]; // we will do symmetrization thus 10 combinations
  double results_aD[10];
  double results_tD[64]; // for the tensor I will print all the combinations even the zeros

  for(int i = 0 ; i < 10; i++){
    results_vD[i]=0.;
    results_aD[i]=0.;
  }
  for(int i = 0 ; i < 64; i++) results_tD[i]=0.;

  // for the derivative we need some kind of mapping because the meaning of the indices is different than the thomas code
  // in QUDA we have the following conventions
  // 0 -> x-direction
  // 1 -> y-direction
  // 2 -> z-direction
  // 3 -> t-direction
  // so we need the following mapping
  // 0 -> t-direction
  // 1 -> x-direction
  // 2 -> y-direction
  // 3 -> z-direction

  int lMap[4];
  lMap[0]=3;
  lMap[1]=0;
  lMap[2]=1;
  lMap[3]=2;

  // Martha needs the data in the following order
  // 00, 10, 11, 20, 21, 22, 30, 31, 32, 33
  

  // vector Derivative according to Martha is real (Warning: this is true without "i" factors)
  for(int i = 0 ; i < NT ; i++){
    results_vD[0] += real(((g4*(*m_gen_der[lMap[0]*NT+i])).trace())*mul_gen); //00
    results_vD[1] += ( real(((g1*(*m_gen_der[lMap[0]*NT+i])).trace())*mul_gen) + real(((g4*(*m_gen_der[lMap[1]*NT+i])).trace())*mul_gen) )/2.; //10
    results_vD[2] += real(((g1*(*m_gen_der[lMap[1]*NT+i])).trace())*mul_gen); //11
    results_vD[3] += ( real(((g2*(*m_gen_der[lMap[0]*NT+i])).trace())*mul_gen) + real(((g4*(*m_gen_der[lMap[2]*NT+i])).trace())*mul_gen) )/2.; //20
    results_vD[4] += ( real(((g2*(*m_gen_der[lMap[1]*NT+i])).trace())*mul_gen) + real(((g1*(*m_gen_der[lMap[2]*NT+i])).trace())*mul_gen) )/2.; //21
    results_vD[5] += real(((g2*(*m_gen_der[lMap[2]*NT+i])).trace())*mul_gen); //22
    results_vD[6] += ( real(((g3*(*m_gen_der[lMap[0]*NT+i])).trace())*mul_gen) + real(((g4*(*m_gen_der[lMap[3]*NT+i])).trace())*mul_gen) )/2.; //30
    results_vD[7] += ( real(((g3*(*m_gen_der[lMap[1]*NT+i])).trace())*mul_gen) + real(((g1*(*m_gen_der[lMap[3]*NT+i])).trace())*mul_gen) )/2.; //31
    results_vD[8] += ( real(((g3*(*m_gen_der[lMap[2]*NT+i])).trace())*mul_gen) + real(((g2*(*m_gen_der[lMap[3]*NT+i])).trace())*mul_gen) )/2.; //32
    results_vD[9] += real(((g3*(*m_gen_der[lMap[3]*NT+i])).trace())*mul_gen); //33
  }
  
  // axial Derivative according to Martha is imag (!!!!!!: Remember this also later wheb you will multiply the loop with the propagator)
  for(int i = 0 ; i < NT ; i++){
    results_aD[0] += imag(((g5g4*(*m_gen_der[lMap[0]*NT+i])).trace())*mul_gen*(-1.)); //00
    results_aD[1] += ( imag(((g5g1*(*m_gen_der[lMap[0]*NT+i])).trace())*mul_gen*(-1.)) + imag(((g5g4*(*m_gen_der[lMap[1]*NT+i])).trace())*mul_gen*(-1.)) )/2.; //10
    results_aD[2] += imag(((g5g1*(*m_gen_der[lMap[1]*NT+i])).trace())*mul_gen*(-1.)); //11
    results_aD[3] += ( imag(((g5g2*(*m_gen_der[lMap[0]*NT+i])).trace())*mul_gen*(-1.)) + imag(((g5g4*(*m_gen_der[lMap[2]*NT+i])).trace())*mul_gen*(-1.)) )/2.; //20
    results_aD[4] += ( imag(((g5g2*(*m_gen_der[lMap[1]*NT+i])).trace())*mul_gen*(-1.)) + imag(((g5g1*(*m_gen_der[lMap[2]*NT+i])).trace())*mul_gen*(-1.)) )/2.; //21
    results_aD[5] += imag(((g5g2*(*m_gen_der[lMap[2]*NT+i])).trace())*mul_gen*(-1.)); //22
    results_aD[6] += ( imag(((g5g3*(*m_gen_der[lMap[0]*NT+i])).trace())*mul_gen*(-1.)) + imag(((g5g4*(*m_gen_der[lMap[3]*NT+i])).trace())*mul_gen*(-1.)) )/2.; //30
    results_aD[7] += ( imag(((g5g3*(*m_gen_der[lMap[1]*NT+i])).trace())*mul_gen*(-1.)) + imag(((g5g1*(*m_gen_der[lMap[3]*NT+i])).trace())*mul_gen*(-1.)) )/2.; //31
    results_aD[8] += ( imag(((g5g3*(*m_gen_der[lMap[2]*NT+i])).trace())*mul_gen*(-1.)) + imag(((g5g2*(*m_gen_der[lMap[3]*NT+i])).trace())*mul_gen*(-1.)) )/2.; //32
    results_aD[9] += imag(((g5g3*(*m_gen_der[lMap[3]*NT+i])).trace())*mul_gen*(-1.)); //33
  }

  // now we want to create all 16 elements of the sigma_munu in order to do the for loops as Thomas does
  vector<Matrix*> sigma_munu;
  for(int i = 0 ; i < 16; i++) sigma_munu.push_back(new Matrix(4,4)); // the constructor allocates memory and initialize them to zero
  *sigma_munu[0*4+1]=ig5Comm_g4g1;
  *sigma_munu[0*4+2]=ig5Comm_g4g2;
  *sigma_munu[0*4+3]=ig5Comm_g4g3;
  *sigma_munu[1*4+0]=ig5Comm_g4g1; sigma_munu[1*4+0]->signFlip(); 
  *sigma_munu[1*4+2]=ig5Comm_g1g2;
  *sigma_munu[1*4+3]=ig5Comm_g1g3;
  *sigma_munu[2*4+0]=ig5Comm_g4g2; sigma_munu[2*4+0]->signFlip();
  *sigma_munu[2*4+1]=ig5Comm_g1g2; sigma_munu[2*4+1]->signFlip();
  *sigma_munu[2*4+3]=ig5Comm_g2g3;
  *sigma_munu[3*4+0]=ig5Comm_g4g3; sigma_munu[3*4+0]->signFlip();
  *sigma_munu[3*4+1]=ig5Comm_g1g3; sigma_munu[3*4+1]->signFlip();
  *sigma_munu[3*4+2]=ig5Comm_g2g3; sigma_munu[3*4+2]->signFlip();

  // accoring to Martha tensor derivative is real. For the tensor we need the standard one-end trick
  for(int i = 0 ; i < NT ; i++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu < 4 ; nu++)
	for(int rho = 0 ; rho < 4 ; rho++){
	  results_tD[mu*4*4+nu*4+rho] += ( real((((*sigma_munu[mu*4+nu])*(*m_std_der[lMap[rho]*NT+i])).trace())*mul_std) + real((((*sigma_munu[mu*4+rho])*(*m_std_der[lMap[nu]*NT+i])).trace())*mul_std) )/2.;
	}


  for(int i = 0 ; i < 10 ; i++){
    results_vD[i] *= (-1.);
    results_aD[i] *= (-1.);
  }

  for(int i = 0 ; i < 64 ; i++)
    results_tD[i] *= (-1.);


  int count=0;
  for(int mu = 0 ; mu < 4; mu++)
    for(int nu = 0 ; nu <= mu ; nu++){
      f_output_vD << scientific << showpos << results_vD[count] << endl;
      count++;
    }

  count = 0 ;
  for(int mu = 0 ; mu < 4; mu++)
    for(int nu = 0 ; nu <= mu ; nu++){
      f_output_aD << scientific << showpos << results_aD[count] << endl;
      count++;
    }

  count = 0;
  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++)
      for(int rho = 0 ; rho < 4 ; rho++){
	f_output_tD << scientific << showpos << results_tD[count] << endl;
	count++;
      }

  //========================== clear memory =====================//
   for(int i = 0 ; i < 16; i++) delete sigma_munu[i];

  delete[] loop_std_raw;
  delete[] loop_gen_raw;
  delete[] loop_std_der_raw;
  delete[] loop_gen_der_raw;

  for(int i = 0; i < NT; i++)
    delete m_std[i];
  delete[] m_std;
  
  for(int i = 0; i < NT; i++)
    delete m_gen[i];
  delete[] m_gen;
  
  for(int i = 0 ; i < 4*NT; i++)
    delete m_std_der[i];
  delete[] m_std_der;

  for(int i = 0 ; i < 4*NT; i++)
    delete m_gen_der[i];
  delete[] m_gen_der;

  return 0;
}

