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
    for(int j=0; j<n; j++){
        x(j) = val(params[j]);
    }

    cout << "\nx_init = \n" << x << "\n" << endl;

    LeastSquaresFunc LS;
    var E;
    VectorXd g(n);
    MatrixXd H;
    double dxnorm, answer;

    answer = 1.53753e-4;

    for(int k=0; k<100; k++){

        E = LS.Kowalik_and_Osborne_function(params);
        cout << sqrt(((val(E) - answer)*(val(E) - answer))) << endl;

        g = Gradient_Xd(n, E, params);
        H = Hesse_Xd(n, E, params);
        
        dx = H.fullPivLu().solve(-g);
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
    cout << "\nH = \n" << H << endl;
    cout << "\nx = \n" << x <<endl;  
    cout << "\n||dx|| = \n" << dxnorm <<endl;

    return 0;
}