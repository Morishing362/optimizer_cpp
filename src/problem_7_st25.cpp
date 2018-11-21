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

    LeastSquaresFunc LS;
    Strategy_2 ST;
    var params_dash[n];
    var E, E_dash;
    var funcvec[m];
    VectorXd e(m);
    VectorXd g(n);
    MatrixXd J, Jt, JtJ, L, I;
    double tau, damp, nu, beta, gamma, p, gf, answer, F;

    I = MatrixXd::Identity(n, n);

    tau = 1.0;
    LS.Freudenstein_and_Roth_function_vectorizer(funcvec, params);
    J = Jacobi_Xd(m, n, funcvec, params);
    Jt = J.transpose();
    JtJ = Jt * J;
    damp = tau * Max_diagonal(JtJ);
    answer = 24.4921;
    nu = 2.0;
    beta = 2.0;
    gamma = 3.0;
    p = 3.0;

    cout << "\nJtJ = \n" << JtJ << "\n" << endl;
    cout << "\nx_init = \n" << x << "\n" << endl;
    cout << "\ndamp_init = \n" << damp << "\n" << endl;

    // cout << "\n" << "E          ||g||          damp       E_dash - E         gf" <<endl;
    for(int k=0; k<70; k++){
        E = LS.Freudenstein_and_Roth_function(params);
        LS.Freudenstein_and_Roth_function_vectorizer(funcvec, params);

        for (int i=0; i<m; i++){
            e(i) = val(funcvec[i]);
        }

        J = Jacobi_Xd(m, n, funcvec, params);
        Jt = J.transpose();
        JtJ = Jt * J;
        g = Jt * e;
        L = JtJ + I*damp;
        dx = L.fullPivLu().solve(-g);
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
        E_dash = LS.Freudenstein_and_Roth_function(params);
        for(int j=0; j<n; j++){
            params[j] = x(j);
        }

        //gain_factor
        gf = Gain_Factor(E, E_dash, damp, g, dx);
        ST.Change_x_params(gf, params, x, dx);
        // cout << damp <<endl;
        cout << g.norm() <<endl;
        // cout << val(E) << endl;
        // cout << gf <<endl;
        // cout << val(E) << "  " << g.norm() << "  " << damp << "  " << val(E_dash) - val(E) << "  " << gf << endl;

        // ダンピングファクタ変更判定
        if (gf > 0){
            damp = damp * ST.max_in_2(beta, gamma, p, gf);
            nu = beta;
        }
        else {
            damp = damp * nu;
            nu = 2.0 * nu;
        }
    }

    cout << "\nJtJ = \n" << JtJ << endl;
    cout << "\nx = \n" << x <<endl;
    cout << "\ndamp = \n" << damp <<endl;
    cout << "\n||dx|| / ||x|| = \n" << dx.norm() / x.norm() <<endl;

    return 0;
}
