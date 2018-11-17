#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include "autodiff.hpp"

using namespace Eigen;
using namespace autodiff;

class LeastSquaresFunc{
  public:
    var Helical_valley_function(const var params[3]);
    void Helical_valley_function_vectorizer(var funcvec[3], const var params[3]);
    var Freudenstein_and_Roth_function(const var params[2]);
    void Freudenstein_and_Roth_function_vectorizer(var funcvec[2], const var params[2]);
    var Kowalik_and_Osborne_function(const var params[4]);
    void Kowalik_and_Osborne_function_vectorizer(var funcvec[11], const var params[4]);
    var Beale_function(const var params[2]);
    void Beale_function_vectorizer(var funcvec[3], const var params[2]);
    var Exponential_fit_2(const var params[4]);
    void Exponential_fit_2_vectorizer(var funcvec[45], const var params[4]);
    double Error(var E, double answer);
};

class Strategy_1{
  private:
    double q1, q2, Alpha, Beta, Gamma;
  public:
    double ep1, ep2;
    Strategy_1();
    void Damper(double gf, double damp);
    void Change_x_params(double gf, var params[], VectorXd& x, VectorXd& dx);
};

VectorXd Gradient_Xd(const int n, var costfunc, const var params[]);
MatrixXd Hesse_Xd(const int n, var costfunc, const var params[]);
MatrixXd Jacobi_Xd(const int m, const int n, var costfuncXd[], const var params[]);

double Max_diagonal(MatrixXd& X);
double Gain_Factor(var E, var E_dash, VectorXd g, VectorXd dx);


// 以下詳細
//最小二乗の様々なプリセット関数をつくるメソッドをまとめたクラス
var LeastSquaresFunc::Helical_valley_function(const var params[3]){
  var er[3];
  var E;
  double th;
  if (val(params[0]) > 0){ th = 0.0; }
  else { th = 5.0; } 
  er[0] = 10*( params[2] - (5/M_PI)*atan(params[1]/params[0]) - th );
  er[1] = 10*( sqrt(params[0]*params[0] + params[1]*params[1]) - 1 );
  er[2] = params[2];
  E = ( er[0]*er[0] + er[1]*er[1] + er[2]*er[2] )/2;
  return E;
};

void LeastSquaresFunc::Helical_valley_function_vectorizer(var funcvec[3], const var params[3]){
  double th;
  if (val(params[0]) > 0){ th = 0.0; }
  else { th = 5.0; }
  funcvec[0] = 10*( params[2] - (5/M_PI)*atan(params[1]/params[0]) - th );
  funcvec[1] = 10*( sqrt(params[0]*params[0] + params[1]*params[1]) - 1 );
  funcvec[2] = params[2];
};

var LeastSquaresFunc::Freudenstein_and_Roth_function(const var params[2]){
  var er[2];
  var E = 0;
  er[0] = params[0] - params[1]*(2 - params[1]*(5 - params[1])) - 13;
  er[1] = params[0] - params[1]*(14 - params[1]*(1 + params[1])) - 29;
  E = ( er[0]*er[0] + er[1]*er[1] )/2;
  return E;
};

void LeastSquaresFunc::Freudenstein_and_Roth_function_vectorizer(var funcvec[2], const var params[2]){
  funcvec[0] = params[0] - params[1]*(2 - params[1]*(5 - params[1])) - 13;
  funcvec[1] = params[0] - params[1]*(14 - params[1]*(1 + params[1])) - 29;
};

var LeastSquaresFunc::Kowalik_and_Osborne_function(const var params[4]){
  var er[11];
  var E;
  double y[11], u[11];
  y[0]  = 0.1957; u[0]  = 4.0000;
  y[1]  = 0.1947; u[1]  = 2.0000;
  y[2]  = 0.1735; u[2]  = 1.0000;
  y[3]  = 0.1600; u[3]  = 0.5000;
  y[4]  = 0.0844; u[4]  = 0.2500;
  y[5]  = 0.0627; u[5]  = 0.1670;
  y[6]  = 0.0456; u[6]  = 0.1250;
  y[7]  = 0.0342; u[7]  = 0.1000;
  y[8]  = 0.0323; u[8]  = 0.0833;
  y[9]  = 0.0235; u[9]  = 0.0714;
  y[10] = 0.0246; u[10] = 0.0625;
  for (int i=0; i<11; i++){
    er[i] = y[i] - ( params[0] * u[i] * ( u[i] + params[1] ) ) 
                    / ( u[i] * ( u[i] + params[2] ) + params[3] );
    E += ( er[i]*er[i] )/2;
  }
  return E;
};

void LeastSquaresFunc::Kowalik_and_Osborne_function_vectorizer(var funcvec[11], const var params[4]){
  double y[11], u[11];
  y[0]  = 0.1957; u[0]  = 4.0000;
  y[1]  = 0.1947; u[1]  = 2.0000;
  y[2]  = 0.1735; u[2]  = 1.0000;
  y[3]  = 0.1600; u[3]  = 0.5000;
  y[4]  = 0.0844; u[4]  = 0.2500;
  y[5]  = 0.0627; u[5]  = 0.1670;
  y[6]  = 0.0456; u[6]  = 0.1250;
  y[7]  = 0.0342; u[7]  = 0.1000;
  y[8]  = 0.0323; u[8]  = 0.0833;
  y[9]  = 0.0235; u[9]  = 0.0714;
  y[10] = 0.0246; u[10] = 0.0625;
  for (int i=0; i<11; i++){
    funcvec[i] = y[i] - ( params[0] * u[i] * ( u[i] + params[1] ) ) 
                    / ( u[i] * ( u[i] + params[2] ) + params[3] );
  }
};

var LeastSquaresFunc::Beale_function(const var params[2]){
  var E;
  E =   (1.5 - params[0] + params[0]*params[1]) * (1.5 - params[0] + params[0]*params[1])/2
      + (2.25 - params[0] + params[0]*params[1]*params[1]) * (2.25 - params[0] + params[0]*params[1]*params[1])/2
      + (2.625 - params[0] + params[0]*params[1]*params[1]*params[1]) * (2.625 - params[0] + params[0]*params[1]*params[1]*params[1])/2;
  return E;
};

void LeastSquaresFunc::Beale_function_vectorizer(var funcvec[3], const var params[2]){
  funcvec[0] = (1.5 - params[0] + params[0]*params[1]);
  funcvec[1] = (2.25 - params[0] + params[0]*params[1]*params[1]);
  funcvec[2] = (2.625 - params[0] + params[0]*params[1]*params[1]*params[1]);
};

var LeastSquaresFunc::Exponential_fit_2(const var params[4]){
  var er[45];
  var E;
  double t = 0.02;
  double y[45];
  y[0]  = 0.090542;  y[1]  = 0.124569;  y[2]  = 0.179367;  y[3]  = 0.195654;
  y[4]  = 0.269707;  y[5]  = 0.286027;  y[6]  = 0.289892;  y[7]  = 0.317475;
  y[8]  = 0.308191;  y[9]  = 0.336995;  y[10]  = 0.348371;  y[11]  = 0.321337;
  y[12]  = 0.299423;  y[13]  = 0.338972;  y[14]  = 0.304763;  y[15]  = 0.288903;
  y[16]  = 0.300820;  y[17]  = 0.303974;  y[18]  = 0.283987;  y[19]  = 0.262078;
  y[20]  = 0.281593;  y[21]  = 0.267531;  y[22]  = 0.218926;  y[23]  = 0.225572;
  y[24]  = 0.200594;  y[25]  = 0.197375;  y[26]  = 0.182440;  y[27]  = 0.183892;
  y[28]  = 0.152285;  y[29]  = 0.174028;  y[30]  = 0.150874;  y[31]  = 0.126220;
  y[32]  = 0.126266;  y[33]  = 0.106384;  y[34]  = 0.118923;  y[35]  = 0.091868;
  y[36]  = 0.128926;  y[37]  = 0.119273;  y[38]  = 0.115997;  y[39]  = 0.105831;
  y[40]  = 0.075261;  y[41]  = 0.068387;  y[42]  = 0.090823;  y[43]  = 0.085205;
  y[44]  = 0.067203;
  for (int i=0; i<45; i++){
    er[i] =  y[i] - ( params[2] * exp(params[0]*t*(i+1)) + params[3] * exp(params[1]*t*(i+1)) );
    E += ( er[i]*er[i] )/2;
  }
  return E;
};
void LeastSquaresFunc::Exponential_fit_2_vectorizer(var funcvec[45], const var params[4]){
  double t = 0.02;
  double y[45];
  y[0]  = 0.090542;  y[1]  = 0.124569;  y[2]  = 0.179367;  y[3]  = 0.195654;
  y[4]  = 0.269707;  y[5]  = 0.286027;  y[6]  = 0.289892;  y[7]  = 0.317475;
  y[8]  = 0.308191;  y[9]  = 0.336995;  y[10]  = 0.348371;  y[11]  = 0.321337;
  y[12]  = 0.299423;  y[13]  = 0.338972;  y[14]  = 0.304763;  y[15]  = 0.288903;
  y[16]  = 0.300820;  y[17]  = 0.303974;  y[18]  = 0.283987;  y[19]  = 0.262078;
  y[20]  = 0.281593;  y[21]  = 0.267531;  y[22]  = 0.218926;  y[23]  = 0.225572;
  y[24]  = 0.200594;  y[25]  = 0.197375;  y[26]  = 0.182440;  y[27]  = 0.183892;
  y[28]  = 0.152285;  y[29]  = 0.174028;  y[30]  = 0.150874;  y[31]  = 0.126220;
  y[32]  = 0.126266;  y[33]  = 0.106384;  y[34]  = 0.118923;  y[35]  = 0.091868;
  y[36]  = 0.128926;  y[37]  = 0.119273;  y[38]  = 0.115997;  y[39]  = 0.105831;
  y[40]  = 0.075261;  y[41]  = 0.068387;  y[42]  = 0.090823;  y[43]  = 0.085205;
  y[44]  = 0.085205;
  for (int i=0; i<45; i++){
    funcvec[i] = y[i] - ( params[2] * exp(params[0]*t*(i+1)) + params[3] * exp(params[1]*t*(i+1)) );
  }
};

double LeastSquaresFunc::Error(var E, double answer){
  return sqrt(((val(E) - answer)*(val(E) - answer)));
};

//Strategy_1のメンバ関数
Strategy_1::Strategy_1(){
  ep1 = 1.0e-10;
  ep2 = 1.0e-10;
  q1 = 0.25;
  q2 = 0.75;
  Alpha = 1.0;
  Beta = 2.0;
  Gamma = 3.0;
};
 
void Strategy_1::Damper(double gf, double damp){
  if( gf < q1 ){
    damp = damp * Beta;
  }
  else if( q2 < gf ){
    damp = damp / Gamma;
  }
};

void Strategy_1::Change_x_params(double gf, var params[], VectorXd& x, VectorXd& dx){
  int n = x.size();
  if ( 0 < gf ){
    x = x + dx;
  }
  for (int i=0; i<n; i++){
    params[i] = x(i);
  }
}

//最小二乗のコスト関数のグラディエント，ヘッセ行列，ヤコビ行列の計算メソッド
VectorXd Gradient_Xd(const int n, var costfunc, const var params[]){
  VectorXd g(n);
  for(int i=0; i<n; i++){
    g(i) = grad(costfunc, params[i]);
  }

  return g;
};

MatrixXd Hesse_Xd(const int n, var costfunc, const var params[]){
  var g[n];
  MatrixXd H(n, n);
  for(int i=0; i<n; i++){
    g[i] = gradx(costfunc, params[i]); 
  }
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      H(i, j) = grad(g[j], params[i]);
    }
  }

  return H;
};

// m:costfuncXdの次元，n:paramsの次元
MatrixXd Jacobi_Xd(const int m, const int n, var costfuncXd[], const var params[]){
  MatrixXd J(m, n);
    for(int i=0; i<m; i++){
      for(int j=0; j<n; j++){
        J(i, j) = grad(costfuncXd[i], params[j]);
    }
  }

  return J;
};

double Max_diagonal(MatrixXd& X){
  double max = 0;
  MatrixXd D = X.diagonal();
  int n = D.rows();
  for (int i=0; i<n; i++){
    if ( max < D(i, 0) ){
      max = D(i, 0);
    }
  }
  return max;
};

double Gain_Factor(var E, var E_dash, double damp, VectorXd& g, VectorXd& dx){
  double gf;
  gf = ( val(E) - val(E_dash) ) / ( ( dx.dot(damp*dx - g) )/2 );
  return gf;
  };