#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include "autodiff.hpp"

using namespace Eigen;
using namespace autodiff;

//最小二乗の様々なプリセット関数をつくるメソッドをまとめたクラス
class LeastSquaresFunc{
  public:
    var Helical_valley_function(const var params[3]);
    void Helical_valley_function_vectorizer(var funcvec[3], const var params[3]);
    var Kowalik_and_Osborne_function(const var params[4]);
    void Kowalik_and_Osborne_function_vectorizer(var funcvec[11], const var params[4]);
    var Beale_function(const var params[2]);
    void Beale_function(var funcvec[3], const var params[2]);
};

//コスト関数のグラディエント，ヘッセ行列，最小二乗コスト関数のヤコビ行列の計算メソッド
// Xd
VectorXd Gradient_Xd(const int n, var costfunc, const var params[]);
MatrixXd Hesse_Xd(const int n, var costfunc, const var params[]);
MatrixXd Jacobi_Xd(const int m, const int n, var costfuncXd[], const var params[]);



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
}

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
}

var LeastSquaresFunc::Beale_function(const var params[2]){
  var E;
  E =   (1.5 - params[0] + params[0]*params[1]) * (1.5 - params[0] + params[0]*params[1])/2
      + (2.25 - params[0] + params[0]*params[1]*params[1]) * (2.25 - params[0] + params[0]*params[1]*params[1])/2
      + (2.625 - params[0] + params[0]*params[1]*params[1]*params[1]) * (2.625 - params[0] + params[0]*params[1]*params[1]*params[1])/2;
  return E;
}

void LeastSquaresFunc::Beale_function(var funcvec[3], const var params[2]){
  funcvec[0] = (1.5 - params[0] + params[0]*params[1]);
  funcvec[1] = (2.25 - params[0] + params[0]*params[1]*params[1]);
  funcvec[2] = (2.625 - params[0] + params[0]*params[1]*params[1]*params[1]);
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