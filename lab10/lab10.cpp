#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<cstdlib>
#include"u_xt_exact.cpp"

using namespace std;

double Vg(double t) {
    double ni = 1.0;
    double t_0 = 7.5;
    double sigma = 0.75;
    return sin(2*M_PI * ni * t) * exp(-(t-t_0)/(2*sigma*sigma));
}

double uniform() {
    return rand()/static_cast<double>(RAND_MAX);
}
double recursive_function_b(double x, double t, double c, double lambda, double mi);

double recursive_function_f(double x, double t, double c, double lambda, double mi) {
    double s = -log(uniform())/(mi + lambda);
    double f_0 = 0.0;
    if (s < t) return lambda/ (lambda + mi) * recursive_function_b(x -c*s, t-s, c, lambda, mi);
    else return f_0;
}

double recursive_function_b(double x, double t, double c, double lambda, double mi) {
    double s = -log(uniform())/(mi + lambda);
    double b_0 = 0.0;
    if (s<t) return lambda/(lambda+mi) * recursive_function_f(x+c*s, t-s, c, lambda, mi);
    else return b_0;
}

void MC_algorithm(int n_paths, double x_start, double t_start, double c, double lambda, double mi, double zeta, double Gamma_g, double Gamma_l, double l, double& fxt, double& bxt) {
    for (int i = 1; i <= 2; i++) {
        double sum = 0;
        for (int n = 0; n < n_paths; n++) {
            double x = x_start;
            double t = t_start;
            double eta = 1.0;
            int sign = pow((-1), i);
            while (t > 0) {
                double s = -log(uniform())/(mi + lambda);
                if (sign == -1) {
                    if (x - c*s > 0) eta *= lambda/(lambda+mi);
                    else {
                        s = x / c;
                        sum += eta * zeta * Vg(t-s);
                        eta *= Gamma_g;
                    }
                    x -= c * s;
                    t -= s;
                }
                else if (sign == 1) {
                    if (x+c*s < l) eta *= lambda/(lambda+mi);
                    else {
                        s = (l-x)/c;
                        eta *= Gamma_l;
                    }
                    x += c * s;
                    t -= s;
                }
                sign *= -1;
            }
        }
        if (i == 1) fxt = sum/static_cast<double>(n_paths);
        else if (i == 2) bxt = sum/static_cast<double>(n_paths);
    }
}

int main(){
    //parametry układu
    //linia
    double L = 0.25, C = 100.0, R = 12.5, G = 0.5, l = 2.0;
    //odbiornik
    double R_l = 12.5, R_g = 75.0;
    //parametry Monte Carlo
    double c = 1.0/sqrt(L*C);
    double mi = G/C;
    double lambda = 1.0/2.0 * (R/L - mi);
    double R_0 = sqrt(L/C);
    double zeta = R_0/(R_0+R_g);
    double Gamma_g = (R_g-R_0)/(R_g+R_0);
    double Gamma_l = (R_l-R_0)/(R_l+R_0);

    //wartosc dokladna
    u_xt_exact();

    //implementacja
    int n_paths[] = {1000, 10000, 100000};
    double t[] = {10.0, 15.0, 25.0, 35.0, 50.0};

    MC_algorithm(n_paths[0], 0.0, 0.0, c, lambda, mi, zeta, Gamma_g, Gamma_l, l, fxt, bxt);

}