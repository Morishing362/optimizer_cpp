#include "../CostfuncgHJ.hpp"
#include <eigen3/Eigen/Dense>

using namespace std;

int main(int argc, char const *argv[]){
    const int n = 2;
    const int m = 3;

    var params[n];
    params[0] = 1;
    params[1] = 1;

// varをeigenに移す
    VectorXd dx(n);
    VectorXd x(n);
    for(int j=0; j<n; j++){
        x(j) = val(params[j]);
    }

    cout << "\nx_init = \n" << x << "\n" << endl;

    LeastSquaresFunc LS;
    var E;
    var funcvec[m];
    VectorXd e(m);
    VectorXd g(n);
    MatrixXd J, Jt, JtJ;
    double dxnorm;

    for(int k=0; k<100; k++){

        E = LS.Beale_function(params);
        LS.Beale_function_vectorizer(funcvec, params);
        cout << E << endl;

        for (int i=0; i<m; i++){
            e(i) = val(funcvec[i]);
        }

        J = Jacobi_Xd(m, n, funcvec, params);
        Jt = J.transpose();
        JtJ = Jt * J;
        g = Jt * e;
        
        dx = JtJ.fullPivLu().solve(-g);
        x = x + dx;
    
        //params更新
        for(int j=0; j<n; j++){
            params[j] = x(j);
        }

        //dxnorm更新
        dxnorm = 0;
        for(int j=0; j<n; j++){
            dxnorm += dx(j)*dx(j);
        } 
        dxnorm = sqrt(dxnorm);

        if ( dxnorm < 1.0e-14){
            break;
        }
    }
    cout << "\nJtJ = \n" << JtJ << endl;
    cout << "\nx = \n" << x <<endl;  
    cout << "\n||dx|| = \n" << dxnorm <<endl;

    return 0;
}