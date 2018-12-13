#include "../CostfuncgHJ.hpp"
#include <eigen3/Eigen/Dense>

using namespace std;

int main(int argc, char const *argv[]){
    const int n = 2;
    const int m = 3;

    var params[n];
    params[0] = 3;
    params[1] = 4;

    // varをeigenに移す
    VectorXd dx(n);
    VectorXd x(n);
    VectorXd x_dash(n);
    for(int j=0; j<n; j++){
        x(j) = val(params[j]);
    }

    LeastSquaresFunc LS;
    Strategy_1 ST;
    var params_dash[n];
    var E, E_dash;
    var funcvec[m];
    VectorXd e(m);
    VectorXd g(n);
    MatrixXd J, Jt, JtJ, L, I;
    double tau, damp, gf, answer, F;

    I = MatrixXd::Identity(n, n);

    tau = 1.0;
    LS.Beale_function_vectorizer(funcvec, params);
    J = Jacobi_Xd(m, n, funcvec, params);
    Jt = J.transpose();
    JtJ = Jt * J;
    damp = tau * Max_diagonal(JtJ);
    answer = 0;
    cout << "\nJtJ = \n" << JtJ << "\n" << endl;
    cout << "\nx_init = \n" << x << "\n" << endl;
    cout << "\ndamp_init = \n" << damp << "\n" << endl;

    // cout << "\n" << "E          ||g||          damp       E_dash - E         gf" <<endl;
    for(int k=0; k<70; k++){
        E = LS.Beale_function(params);
        LS.Beale_function_vectorizer(funcvec, params);

        for (int i=0; i<m; i++){
            e(i) = val(funcvec[i]);
        }

        J = Jacobi_Xd(m, n, funcvec, params);
        Jt = J.transpose();
        JtJ = Jt * J;
        g = Jt * e;
        L = JtJ + I*damp;
        dx = L.colPivHouseholderQr().solve(-g);
        // 収束判定
        if ( g.norm() < 1.0e-12 ){
            break;
        }
        if( dx.norm()/x.norm() < 1.0e-12 ){
            break;
        }
    
        x_dash = x + dx;
        for(int j=0; j<n; j++){
            params[j] = x_dash(j);
        }
        E_dash = LS.Beale_function(params);
        for(int j=0; j<n; j++){
            params[j] = x(j);
        }

        //gain_factor
        gf = Gain_Factor(E, E_dash, damp, g, dx);
        ST.Change_x_params(gf, params, x, dx);
        // cout << damp <<endl;
        // cout << g.norm() <<endl;
        // cout << abs(answer - val(E)) << endl;
        // cout << gf <<endl;
        // cout << val(E) << "  " << g.norm() << "  " << damp << "  " << val(E_dash) - val(E) << "  " << gf << endl;
        // cout << val(E) << endl;
        cout << x(0) << endl;

        // ダンピングファクタ変更判定
        if( gf < 0.2 ){
            damp = damp * 2.0;
        }
        else if( 0.8 < gf ){
            damp = damp / 3.0;
        }
    }

    cout << "\nJtJ = \n" << JtJ << endl;
    cout << "\nx = \n" << x <<endl;
    cout << "\ndamp = \n" << damp <<endl;
    cout << "\n||dx|| / ||x|| = \n" << dx.norm() / x.norm() <<endl;

    return 0;
}



