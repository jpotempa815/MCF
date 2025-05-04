#include<iostream>
#include<cmath>
#include<cstdlib>
#include<string>
#include<vector>
#include<fstream>

using namespace std;

double uniform() {
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

void gillespie(int P_max, ofstream& out, ofstream& out2) {
    double x1_0 = 120, x2_0 = 80, x3_0 = 1;
    double k1 = 1, k2 = 1, k3 = 0.001, k4 = 0.01;
    int t_max = 200;

    double Gamma1 = 0, Gamma2 = 0, Gamma3 = 0, Gamma4 = 0, Gamma_max = 0;

    int N = 50;
    double dt = t_max/static_cast<double>(N);

    vector<double> h1(N, 0.0);
    vector<double> h2(N, 0.0);
    int m = 0;
    double t = 0;
    double x3_sr = 0, x3_sr2 = 0, sigma_x3t = 0;

    for (int p = 1; p <= P_max; p++) {
        t = 0;
        double x1 = x1_0, x2 = x2_0, x3 = x3_0;
        vector<double> h0(N, 0.0);
        vector<int> ncount(N, 0);
        while (t < t_max) {
            int l = 0;
            Gamma1 = k1;
            Gamma2 = k2;
            Gamma3 = k3 * x1 * x2;
            Gamma4 = k4 * x3;
            Gamma_max = Gamma1 + Gamma2 + Gamma3 + Gamma4;
            double U1 = uniform();
            double Dt;
            if (U1 <= 1e-7 || Gamma_max <= 0) Dt = 0;
            else Dt = -log(U1) / Gamma_max;
            double U2 = uniform();
            if ( U2 <= Gamma1/Gamma_max) m = 1;
            else if (U2 <= (Gamma1 + Gamma2)/Gamma_max) m = 2;
            else if (U2 <= (Gamma1 + Gamma2 + Gamma3)/Gamma_max) m = 3;
            else m = 4;
            if (m == 1) x1 += 1;
            else if (m == 2) x2 += 1;
            else if (m == 3) {
                x1 -= 1;
                x2 -= 1;
                x3 += 1;
            }
            else if (m == 4) {
                x3 -= 1;
            }
            t += Dt;
            if ( P_max == 1) out << t << "," << x1 << "," << x2 << "," << x3 << endl;
            else out << p << "," << t << "," << x1 << "," << x2 << "," << x3 << endl;
            //wklad sciezki
            l = floor(t/dt);
            if (l >= N) {
                l = N - 1;
            }
            h0[l] += x3;
            ncount[l]++;
        }
        //usrednianie po sciezce
        double x3_t = 0;
        for (int j = 0; j < N; ++j) {
            if (ncount[j] > 0) {
                x3_t = h0[j] / static_cast<double>(ncount[j]);
            }
            h1[j] += x3_t;
            h2[j] += x3_t * x3_t;
        }
    }

    //wszystkie sciezki
    for (int l = 0; l < N; ++l){
        double t_sr;
        x3_sr = h1[l]/static_cast<double>(P_max);
        x3_sr2 = h2[l]/static_cast<double>(P_max);
        // cout << x3_sr*x3_sr << " " << x3_sr2 << endl;
        sigma_x3t = sqrt((x3_sr2 - x3_sr*x3_sr)/static_cast<double>(P_max));
        t_sr = (l + 0.5) * dt;
        out2 << t_sr << "," << x3_sr << "," << sigma_x3t << endl;
    }
}

int main() {
    ofstream test;
    ofstream test2;
    test.open("test.txt");
    test2.open("test2.txt");
    gillespie(1, test, test2);
    test.close();
    test2.close();
    ofstream out1;
    out1.open("Pmax=5.txt");
    ofstream out1_stat;
    out1_stat.open("Pmax=5_stat.txt");
    gillespie(5, out1, out1_stat);
    out1.close();
    out1_stat.close();
    ofstream out2;
    out2.open("Pmax=100.txt");
    ofstream out2_stat;
    out2_stat.open("Pmax=100_stat.txt");
    gillespie(100, out2, out2_stat);
    out2.close();
    out2_stat.close();
    return 0;
}