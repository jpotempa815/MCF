#include<iostream>
#include<vector>
#include<cstdlib>
#include<cmath>
#include<fstream>

using namespace std;

double uniform() {
    return (double) rand() / (double) RAND_MAX;
}

void zapis_do_pliku(ostream& plik, vector<double>& x, vector<double>& y, int N) {
    for (int i=0; i<N; i++) {
        plik << x[i] << " " << y[i] << endl;
    }
}

void kolo(int N, vector<double>& x, vector<double>& y, double R, double xa) {
    for (int i = 0; i < N; i++) {
        double U_1 = uniform();
        double U_2 = uniform();
        x[i] = sqrt(-2.* log(1-U_1))* cos(2*M_PI*U_2);
        y[i] = sqrt(-2.* log(1-U_1)) * sin(2*M_PI*U_2);
        double r = sqrt(x[i]*x[i] + y[i]*y[i]);
        x[i] = x[i]/r;
        y[i] = y[i]/r;
        double q = sqrt(uniform());
        x[i] = q * x[i] * R + xa;
        y[i] = q * y[i] * R;
    }
}

void powierzchnia(int N, ofstream& plik, vector<double>& x, vector<double>& y, double R, double Ra, double Rb, double xa, double ya, double xb, double yb) {
    double mi_1 = 0, mi_2 = 0, sigma = 0, sigma_2 = 0;
    for (int i = 1; i <= N; i++) {
        double theta = 0;
        double odlegloscA = sqrt(pow(x[i-1] - xa, 2) + pow(y[i-1] - ya, 2));
        double odlegloscB = sqrt(pow(x[i-1] - xb, 2) + pow(y[i-1] - yb, 2));

        if (odlegloscA <= Ra && odlegloscB <= Rb) {
            theta = 1;
        }
        mi_1 += M_PI * R * R * theta;
        if (i == 100 || i == 1000 || i == 10000 || i == 100000 || i == 1000000) {
            double mi_s = mi_1/ (double)i;
            mi_2 = M_PI * R * R * mi_s;
            sigma_2 = (mi_2 - mi_s * mi_s) / (double)i;
            sigma = sqrt(sigma_2);
            plik << i << " " << mi_s << " " << sigma << endl;
        }
    }
}

int main(){
    double Ra {2}, Rb {0}, xa {0};
    Rb = sqrt(2.) * Ra;
    int N = 10000;
    vector<double> x(N), y(N);
    //test generatora
    ofstream test, test2;
    test.open("test_generatora_1.txt");
    test2.open("test_generatora_2.txt");
    xa = Ra + Rb;
    kolo(N, x, y, Ra, xa);
    zapis_do_pliku(test, x, y, N);
    kolo(N, x, y, Rb, 0);
    zapis_do_pliku(test2, x, y, N);
    test.close();
    test2.close();

    //powierzchnia wspolna
    N = 1000000;

    ofstream plik1;
    plik1.open("pow_A_xa.txt");
    xa = Rb + 0.5*Ra;
    vector<double> x1(N), y1(N);
    kolo(N, x1, y1, Ra, xa);
    powierzchnia(N, plik1, x1, y1, Ra, Ra, Rb, xa, 0, 0, 0);
    plik1.close();

    ofstream plik2;
    plik2.open("pow_B_xa.txt");
    kolo(N, x1, y1, Rb, xa);
    powierzchnia(N, plik2, x1, y1, Rb, Ra, Rb, xa, 0, 0, 0);
    plik2.close();

    ofstream plik3;
    plik3.open("pow_A_0.txt");
    xa = 0;
    kolo(N, x1, y1, Ra, xa);
    powierzchnia(N, plik3, x1, y1, Ra, Ra, Rb, xa, 0, 0, 0);
    plik3.close();

    ofstream plik4;
    plik4.open("pow_B_0.txt");
    kolo(N, x1, y1, Rb, xa);
    powierzchnia(N, plik4, x1, y1, Rb, Ra, Rb, xa, 0, 0, 0);
    plik4.close();
}
