#include<iostream>
#include<cmath>
#include<cstdlib>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

const int N {100000};

double uniform(){
    return rand() / (double) RAND_MAX;
}

double gaussian(){
    double u1 = uniform();
    double u2 = uniform();
    return sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

void wiener(double dt, double t_max, ofstream& out_1, ofstream& out_2) {
    double D = 1., x_sr = 0., y_sr = 0., xy_sr = 0., x2_sr = 0., y2_sr = 0.;
    vector<double> x(N), y(N), t(N);
    x[0] = 0., y[0] = 0.;
    int Nt = static_cast<int>(t_max/dt);
    double sigma = sqrt(2. * D * dt);
    //wsp dyfuzji
    double Dxx = 0., Dyy = 0., Dxy = 0., Dxx_sr = 0., Dyy_sr = 0., Dxy_sr = 0.;

    for (int k = 1; k <= Nt; k++) {
        double t = k * dt;
        for (int i = 0; i < N; i++) {
            x[i] += sigma * gaussian();
            y[i] += sigma * gaussian();
            out_1 << x[i] << "," << y[i] << endl;
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
        Dxx_sr += Dxx;
        Dyy_sr += Dyy;
        Dxy_sr += Dxy;
    }
    Dxx_sr /= static_cast<double>(Nt);
    Dyy_sr /= static_cast<double>(Nt);
    Dxy_sr /= static_cast<double>(Nt);
    out_2 << Dxx_sr << "," << Dyy_sr << "," << Dxy_sr << endl;
}

int main(){
    double dt = 0.01, t_max = 10.;
    ofstream wiener_xy("wiener_xy.csv");
    ofstream wiener_D("wiener_D.csv");
    wiener(dt, t_max, wiener_xy, wiener_D);
    wiener_xy.close();
    wiener_D.close();
}
