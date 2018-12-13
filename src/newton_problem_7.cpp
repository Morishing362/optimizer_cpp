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
    for(int j=0; j<n; j++){
        x(j) = val(params[j]);
    }

    cout << "\nx_init = \n" << x << "\n" << endl;

    LeastSquaresFunc LS;
    var E;
    VectorXd g(n);
    MatrixXd H;
    double dxnorm, answer;

    for(int k=0; k<100; k++){

        E = LS.Freudenstein_and_Roth_function(params);
        // cout << val(E) << endl;
        // cout << x(0) << endl;
        cout << "[" << x(0) << ", " << x(1) << "]," <<endl;


        g = Gradient_Xd(n, E, params);
        H = Hesse_Xd(n, E, params);
        
        dx = H.fullPivLu().solve(-g);
        x = x + dx;

        // 収束判定
        if ( g.norm() < 1.0e-12 ){
            break;
        }
        if( dx.norm()/x.norm() < 1.0e-12 ){
            break;
        }
    
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
    }
    cout << "\nH = \n" << H << endl;
    cout << "\nx = \n" << x <<endl;  
    cout << "\n||dx|| = \n" << dxnorm <<endl;

    return 0;
}