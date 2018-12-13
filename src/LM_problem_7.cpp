#include "../CostfuncgHJ.hpp"
#include <eigen3/Eigen/Dense>

using namespace std;

int main(int argc, char const *argv[]){
    const int n = 2;
    const int m = 2;

    var params[n];
    params[0] = 0.5;
    params[1] = -2;

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
    double dxnorm, damp, F;

    damp = 1;
    I = MatrixXd::Identity(n, n);
    F = val(LS.Freudenstein_and_Roth_function(params));

    cout << "\ndamp_init = \n" << damp << "\n" << endl;

    for(int k=0; k<1000; k++){
        E = LS.Freudenstein_and_Roth_function(params);
        LS.Freudenstein_and_Roth_function_vectorizer(funcvec, params);
        // cout << x(1) << endl;
        // cout << val(E) << endl;
        cout << "[" << x(0) << ", " << x(1) << "]," <<endl;


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
        
        // 収束判定
        if ( g.norm() < 1.0e-12 ){
            break;
        }
        if( dx.norm()/x.norm() < 1.0e-12 ){
            break;
        }
    
        // params更新
        for(int j=0; j<n; j++){
            params[j] = x_dash(j);
        }
        E_plus_1 = LS.Freudenstein_and_Roth_function(params);

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
    }
    
    cout << "\nJtJ = \n" << JtJ << endl;
    cout << "\nx = \n" << x <<endl;  
    cout << "\ndamp = \n" << damp <<endl;
    cout << "\n||dx|| = \n" << dxnorm <<endl;

    return 0;
}