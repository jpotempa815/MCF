#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;


double uniform() {
    return rand()/static_cast<double>(RAND_MAX);
}

double psi(double r, double a, double c) {
    return (1+c*r)*exp(-a*r);
}

double eps_loc(double r, double a, double c) {
    return (-a*a*c*r*r + (-a*a+4*a*c-2*c)*r + 2*a - 2*c - 2)/(2*c*r*r + 2*r);
}

double p_acc(double ri, double r_new, double a, double c) {
    return pow(r_new/ri, 2) * pow(psi(r_new, a, c), 2)/(pow(psi(ri, a, c), 2));
}

vector<double> energy(int N, double r, double a, double c) {
    vector<double> stat_param(2);
    stat_param[0] = 0.0;
    stat_param[1] = 0.0;
    double dr = 0.1;
    double ri = r;
    double eps_mean = 0.0, eps_mean2 = 0.0;
    for (int i = 0; i < N; i++) {
        double r_new = ri + dr*(2*uniform() - 1);
        double U2 = uniform();
        double p = min(p_acc(ri, r_new, a, c), 1.0);
        if (p >= U2 && r_new > 0.0) ri = r_new;
        eps_mean += eps_loc(ri, a, c);
        eps_mean2 += eps_loc(ri, a, c) * eps_loc(ri, a, c);
    }
    eps_mean /= static_cast<double>(N);
    eps_mean2 /= static_cast<double>(N);
    stat_param[0] = eps_mean;
    stat_param[1] = eps_mean2 - eps_mean*eps_mean;
    return stat_param;
}

vector<double> histogram(int N, double r, double a, double c) {
    double dr = 0.1;
    double r_max = 8.0;
    double ri = r;
    int M = 200;
    double delt_r = r_max / static_cast<double>(M);
    vector<double> histogram(M, 0.0);
    for (int i = 0; i < N; i++) {
        double r_new = ri + dr*(2*uniform() - 1);
        double p = min(p_acc(ri, r_new, a, c), 1.0);
        double U2 = uniform();
        if (p >= U2 && r_new > 0) ri = r_new;
        int k = floor(ri/delt_r);
        if(k < 200 && k >= 0) histogram[k] += 1;
    }
    for(int m = 0; m < M; m++){
        histogram[m] = histogram[m]/(static_cast<double>(N)*delt_r);
    }
    return histogram;
}

int main() {
    int N = 1000000;
    double a = 0.3;
    double c = -0.7;
    double r = 0.1;
    ofstream file("energy_stat.txt", ios::app);
    file << "a c średnia_energia wariancja log10" << endl;
    do {
        c = -0.7;
        do {
            file << a << " " << c << " " << energy(N, r, a, c)[0] << " " << energy(N, r, a, c)[1] << " " << log10(energy(N, r, a, c)[1] + 10e-20) << endl;
            c += 0.02;
        }while (c <= 0.3);
        a += 0.02;
    }while (a <= 1.2);
    file.close();
    //histogram
    ofstream hist_file("histogram.txt");
    a = 1.0;
    c = 0.0;
    r = 0.0;
    vector<double> hist = histogram(N, r, a, c);
    for (int i = 0; i < hist.size(); i++) {
        hist_file << hist[i] << endl;
    }
    hist_file.close();
    ofstream hist_file_2("histogram_200.txt");
    a = 0.5;
    c = -0.5;
    vector<double> hist_200 = histogram(N, r, a, c);
    for(int i = 0; i < hist_200.size(); i++) {
        hist_file_2 << hist_200[i] << endl;
    }
    hist_file_2.close();
    ofstream hist_file_3("histogram_200_2.txt");
    vector<double> hist_200_2 = histogram(10e7, r, a, c);
    for(int i = 0; i < hist_200_2.size(); i++) {
        hist_file_3 << hist_200_2[i] << endl;
    }
    hist_file_3.close();
    ofstream hist_file_4("histogram_200_3.txt");
    vector<double> hist_200_3 = histogram(10e8, r, a, c);
    for(int i = 0; i < hist_200_3.size(); i++) {
        hist_file_4 << hist_200_3[i] << endl;
    }
    hist_file_4.close();
}