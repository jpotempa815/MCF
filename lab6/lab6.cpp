#include<iostream>
#include<cmath>
#include<cstdlib>
#include<fstream>
#include<string>
#include<vector>
#include <random>
#include "particle_translation.cpp"
#include <sstream>

using namespace std;

// const int N {100000};

double uniform(){
    return rand() / (double) RAND_MAX;
}

double gaussian(){
    double u1 = uniform();
    double u2 = uniform();
    return sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

void wiener() {
    double D = 1.;
    double dt = 0.01, t_max = 10.;
    double sigma = sqrt(2. * D * dt);

    //generator liczb z rozkladu normalnego
    mt19937_64 rd(42);
    normal_distribution<double> gauss(0.0, sigma);

    vector<int> Nmax_vals = {100, 1000, 10000, 100000};

    int Nt = static_cast<int>(t_max/dt);

    //wsp dyfuzji
    double Dxx = 0., Dyy = 0., Dxy = 0.;
    double var_Dxx = 0., var_Dyy = 0., var_Dxy = 0.;
    for (auto N : Nmax_vals) {
        string fname_xy   = "wiener_xy_N" + to_string(N) + ".csv";
        string fname_D    = "wiener_D_N" + to_string(N) + ".csv";
        ofstream out_xy(fname_xy);
        ofstream out_D (fname_D);
        vector<double> x(N, 0.0), y(N, 0.0);
        double Dxx_sr = 0., Dyy_sr = 0., Dxy_sr = 0.;
        double Dxx_sr2 = 0., Dyy_sr2 = 0., Dxy_sr2 = 0.;

        for (int k = 1; k <= Nt; k++) {
            double t = k * dt;
            double x_sr = 0., y_sr = 0., xy_sr = 0., x2_sr = 0., y2_sr = 0.;
            for (int i = 0; i < N; i++) {
                x[i] += gauss(rd);
                y[i] += gauss(rd);
                if (fabs(t-0.1) < 1e-9 || fabs(t-1) < 1e-9 || fabs(t-5) < 1e-9) {
                    out_xy << x[i] << "," << y[i] << endl;
                }
                x_sr += x[i];
                y_sr += y[i];
                xy_sr += x[i] * y[i];
                x2_sr += x[i] * x[i];
                y2_sr += y[i] * y[i];
            }
            x_sr /= static_cast<double>(N);
            y_sr /= static_cast<double>(N);
            xy_sr /= static_cast<double>(N);
            x2_sr /= static_cast<double>(N);
            y2_sr /= static_cast<double>(N);
            Dxx = (x2_sr - x_sr * x_sr) / (2*t);
            Dyy = (y2_sr - y_sr * y_sr) / (2*t);
            Dxy = (xy_sr - x_sr * y_sr) / (2*t);
            out_D << Dxx << "," << Dyy << "," << Dxy << endl;
            Dxx_sr += Dxx;
            Dyy_sr += Dyy;
            Dxy_sr += Dxy;
            Dxx_sr2 += Dxx * Dxx;
            Dyy_sr2 += Dyy * Dyy;
            Dxy_sr2 += Dxy * Dxy;
        }
        Dxx_sr /= static_cast<double>(Nt);
        Dyy_sr /= static_cast<double>(Nt);
        Dxy_sr /= static_cast<double>(Nt);
        Dxx_sr2 /= static_cast<double>(Nt);
        Dyy_sr2 /= static_cast<double>(Nt);
        Dxy_sr2 /= static_cast<double>(Nt);

        var_Dxx = sqrt((Dxx_sr2 - Dxx_sr * Dxx_sr) / static_cast<double>(Nt));
        var_Dyy = sqrt((Dyy_sr2 - Dyy_sr * Dyy_sr) / static_cast<double>(Nt));
        var_Dxy = sqrt((Dxy_sr2 - Dxy_sr * Dxy_sr) / static_cast<double>(Nt));
        cout << "N: " << N << " Dxx: " << Dxx_sr << " +/- " << var_Dxx << endl;
        cout << "N: " << N << " Dyy: " << Dyy_sr << " +/- " << var_Dyy << endl;
        cout << "N: " << N << " Dxy: " << Dxy_sr << " +/- " << var_Dxy << endl;
        out_xy.close();
        out_D.close();
    }
}

bool point_out(double x, double y, double cx, double cy, double R) {
    return sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy)) > R;
}

void dyf_abs(double Ra, double w) {
    double D = 1.;
    double dt = 0.01, t_max = 1000.;
    double sigma = sqrt(2. * D * dt);

    string fname_xy;
    string fname_N;

    if (Ra == 0.1) {
        fname_N = "dyf_abs_N_Ra0.1_w" + to_string(static_cast<int>(w)) + ".csv";
    }
    else {
        fname_N = "dyf_abs_N_Ra0.5_w" + to_string(static_cast<int>(w)) + ".csv";
    }
    ofstream out_N (fname_N);

    //generator liczb z rozkladu normalnego
    mt19937_64 rd(42);
    normal_distribution<double> gauss(0.0, sigma);
    //srodek i promien obszaru
    double xr = 0., yr = 0., Rr = 5.;
    //srodek obszaru absorpcyjnego
    double xa = 3., ya = 0.;
    //polozenie zrodla
    double xs = -4.5, ys = 0.;
    int N_max = 10000;
    int N = static_cast<int>(t_max/dt);
    double dn {0};
    //obiekt czastki
    struct particle {
        double x;
        double y;
        int theta;
    };
    //tablica czastek
    vector<particle> particles(N_max, {0.,0.,0});

    //punkty startowy i koncowy
    double x1 = 0, y1 = 0;
    double x2 = 0, y2 = 0;


    //petla czasowa
    for (int it = 0; it < N; it++) {
        dn += static_cast<double>(w) * dt;
        double t = it * dt;
        int n_new = 0; //dodawane ze zrodla
        int n = 0; //aktywne

        //sprawdzenie i dodanie nowych czastek
        for (int i =1; i <= N_max; i++) {
            if (particles[i].theta == 0 && n_new < dn) {
                particles[i].theta = 1;
                particles[i].x = xs;
                particles[i].y = ys;
                n_new++;
            }
        }
        dn -= n_new;

        //przemieszczanie
        for(auto &p: particles) if(p.theta==1) {
                x1 = p.x;
                y1 = p.y;
                x2 = p.x + gauss(rd);
                y2 = p.y + gauss(rd);
                double length=0;
                if (!point_out(x2, y2, xr, yr, Rr)) {
                    if (!point_out(x2, y2, xa, ya, Ra)) {
                        p.theta = 0;
                    }
                    else {
                        p.x = x2;
                        p.y = y2;
                    }
                }
                else {
                    int theta = 1;
                    do {
                        particle_translation(x1, y1, x2, y2, xr, yr, Rr, xa, ya, Ra, theta, length);
                    }while (length > 1e-6 && theta == 1);
                    p.theta = theta;
                    if (theta == 1) {
                        p.x = x1;
                        p.y = y1;
                    }
                }
            if (p.theta == 1) n++;
        }
        out_N << t << "," << n << endl;
        if (fabs(t - 0.1) < 1e-9 || fabs(t - 1) < 1e-9 || fabs(t - 100) < 1e-9 || fabs(t - 400) < 1e-9) {
            ostringstream oss;
            oss <<"dyf_abs_xy_Ra" << Ra << "_w" << w << "_t" << t << ".csv";
            ofstream out_xy(oss.str());
            for ( auto &p: particles) {
                if (p.theta == 1) {
                    out_xy << p.x << "," << p.y << endl;
                }
            }
            out_xy.close();
        }
    }
    out_N.close();
}


int main(){
    wiener();
    // vector<double> Ra_val = {0.1, 0.5};
    // vector<int> w_val = {10, 50, 100};
    // for (double Ra: Ra_val) {
    //     for (double w: w_val) {
    //         dyf_abs(Ra, w);
    //     }
    // }
}
