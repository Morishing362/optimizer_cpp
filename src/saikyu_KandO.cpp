#include "../CostfuncgHJ.hpp"
#include <eigen3/Eigen/Dense>

using namespace std;

int main(int argc, char const *argv[]){
    const int n = 4;
    const int m = 11;

    var params[4];
    params[0] = 0.25;
    params[1] = 0.39;
    params[2] = 0.415;
    params[3] = 0.39;

    // varをeigenに移す
    VectorXd dx(4);
    VectorXd x(4);
    for(int j=0; j<4; j++){
        x(j) = val(params[j]);
    }
    // 初期値を出力
    cout << "\nx_init = \n" << x << "\n" << endl;

    LeastSquaresFunc LS;
    var E, E_plus_1;
    var funcvec[11];
    VectorXd e(11);
    VectorXd g(4);
    MatrixXd J, Jt, JtJ, L, I;
    double dxnorm, damp;

    damp = 1;
    I = MatrixXd::Identity(4, 4);

    for(int k=0; k<10000; k++){
        E = LS.Kowalik_and_Osborne_function(params);
        LS.Kowalik_and_Osborne_function_vectorizer(funcvec, params);
        if(k % 100  ==0){cout << E << endl;}

        for (int i=0; i<11; i++){
            e(i) = val(funcvec[i]);
        }

        J = Jacobi_Xd(m, n, funcvec, params);
        Jt = J.transpose();
        JtJ = Jt * J;
        g = Jt * e;

        x = x - 0.001*g;
    
        // params更新
        for(int j=0; j<4; j++){
            params[j] = x(j);
        }


        // dxnorm更新
        dxnorm = 0;
        for(int j=0; j<4; j++){
            dxnorm += dx(j)*dx(j);
        } 
        dxnorm = sqrt(dxnorm);

        // // 収束判定
        // if ( dxnorm < 1.0e-14){
        //     break;
        // }
    }
    
    cout << "\nJtJ = \n" << JtJ << endl;
    cout << "\nx = \n" << x <<endl;  
    cout << "\ndamp = \n" << damp <<endl;
    cout << "\n||dx|| = \n" << dxnorm <<endl;

    return 0;
}