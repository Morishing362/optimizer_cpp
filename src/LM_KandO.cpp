#include "../CostfuncgHJ.hpp"
#include <eigen3/Eigen/Dense>

using namespace std;

int main(int argc, char const *argv[]){
    const int n = 4;
    const int m = 11;

    var params[n];
    params[0] = 0.25;
    params[1] = 0.39;
    params[2] = 0.415;
    params[3] = 0.39;

    // varをeigenに移す
    VectorXd dx(n);
    VectorXd x(n);
    VectorXd x_dash(n);
    for(int j=0; j<n; j++){
        x(j) = val(params[j]);
    }
    // 初期値を出力
    cout << "\nx_init = \n" << x << "\n" << endl;

    LeastSquaresFunc LS;
    var E, E_plus_1;
    var funcvec[m];
    VectorXd e(m);
    VectorXd g(n);
    MatrixXd J, Jt, JtJ, L, I;
    double dxnorm, damp, F, answer;

    damp = 1;
    I = MatrixXd::Identity(n, n);
    F = val(LS.Kowalik_and_Osborne_function(params));
    answer = 1.53753e-4;

    cout << "\ndamp_init = \n" << damp << "\n" << endl;

    for(int k=0; k<100; k++){
        E = LS.Kowalik_and_Osborne_function(params);
        LS.Kowalik_and_Osborne_function_vectorizer(funcvec, params);
        cout << sqrt(((val(E) - answer)*(val(E) - answer))) << endl;

        for (int i=0; i<m; i++){
            e(i) = val(funcvec[i]);
        }

        J = Jacobi_Xd(m, n, funcvec, params);
        Jt = J.transpose();
        JtJ = Jt * J;
        g = Jt * e;
        L = JtJ + I*damp;

        dx = L.fullPivLu().solve(-g);
        x_dash = x + dx;
    
        // params更新
        for(int j=0; j<n; j++){
            params[j] = x_dash(j);
        }
        E_plus_1 = LS.Kowalik_and_Osborne_function(params);

        // ダンピングファクタ変更判定
        if( val(E_plus_1) > F ){
            damp = 10 * damp;
            for(int j=0; j<n; j++){
                params[j] = x(j); // params戻す
            }
        }
        else{
            damp = 0.1 * damp;
            F = val(E_plus_1);
            x = x + dx;
        }

        // dxnorm更新
        dxnorm = 0;
        for(int j=0; j<n; j++){
            dxnorm += dx(j)*dx(j);
        } 
        dxnorm = sqrt(dxnorm);

        // 収束判定
        if ( dxnorm < 1.0e-14){
            break;
        }
    }
    
    cout << "\nJtJ = \n" << JtJ << endl;
    cout << "\nx = \n" << x <<endl;  
    cout << "\ndamp = \n" << damp <<endl;
    cout << "\n||dx|| = \n" << dxnorm <<endl;

    return 0;
}