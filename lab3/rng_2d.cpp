#include<iostream>
#include<vector>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include<math.h>

using namespace std;

// const int N = 10000;
const int N = 10000;

double uniform() {
    return (double) rand() / (double) RAND_MAX;
}

void zapis_do_pliku(ostream& plik, vector<double>& x, vector<double>& y) {
    for (int i=0; i<N; i++) {
            plik << x[i] << " " << y[i] << endl;
    }
}

void Box_Muller(ofstream& plik, vector<double>& x, vector<double>& y) {
    for (int i = 0; i < N; i++) {
        double U_1 = uniform();
        double U_2 = uniform();
        x[i] = sqrt(-2.* log(1-U_1))* cos(2*M_PI*U_2);
        y[i] = sqrt(-2.* log(1-U_1)) * sin(2*M_PI*U_2);
    }
    zapis_do_pliku(plik, x, y);
}

void rozklad_kolo(ofstream& plik, vector<double>& x, vector<double>& y) {
    double R {0};
    for (int i = 0; i < N; i++) {
        double r = sqrt(x[i]*x[i] + y[i]*y[i]);
        x[i] = x[i]/r;
        y[i] = y[i]/r;
        R = sqrt(uniform());
        x[i] = R * x[i];
        y[i] = R * y[i];
    }
    zapis_do_pliku(plik, x, y);
}

void afiniczna(ofstream& plik, double alpha, double b1, double b2, vector<double>& x, vector<double>& y) {
    vector<double> x_elip(N), y_elip(N);

    for (int i = 0; i < N; i++) {
        x_elip[i] = (b1 * cos(alpha) * x[i] - b2 * sin(alpha) * y[i]);
        y_elip[i] = (b1 * sin(alpha) * x[i] + b2 * cos(alpha) * y[i]);
    }
    x = x_elip;
    y = y_elip;
    zapis_do_pliku(plik, x_elip, y_elip);
}

void kowariancja(ofstream& plik, vector<double> x, vector<double> y) {
    vector<vector<double>> kowariancja(2, vector<double>(2));
    double x_sr {0}, x_sr_2 {0}, y_sr {0}, y_sr_2 {0}, xy_sr {0};
    double sigma_x2 {0}, sigma_y2 {0}, sigma_xy {0};
    for (int i = 0; i < N; i++) {
        x_sr += x[i];
        x_sr_2 += pow(x[i], 2);
        y_sr += y[i];
        y_sr_2 += pow(y[i],2);
        xy_sr += x[i]*y[i];
    }
    sigma_x2 = x_sr_2/(double)N - (x_sr*x_sr)/((double)N*N);
    sigma_y2 = y_sr_2/(double)N - (y_sr*y_sr)/((double)N*N);
    sigma_xy = xy_sr/(double)N - (x_sr*y_sr)/((double)N*N);
    kowariancja[0][0] = sigma_x2;
    kowariancja[1][1] = sigma_y2;
    kowariancja[0][1] = kowariancja[1][0] = sigma_xy;

    double r_xy = sigma_xy/sqrt(sigma_x2*sigma_y2);

    plik << r_xy << endl;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            plik << kowariancja[i][j] << " ";
        }
        plik << endl;
    }
}


int main(){
    vector<double> x_jedn(N), y_jedn(N);
    vector<double> x_norm(N), y_norm(N);
    //normalny w dwoch wymiarach
    ofstream box_muller;
    box_muller.open("box_muller.txt");
    Box_Muller(box_muller, x_norm, y_norm);
    box_muller.close();
    //kolo
    ofstream kolo;
    kolo.open("kolo.txt");
    Box_Muller(box_muller, x_jedn, y_jedn);
    rozklad_kolo(kolo, x_jedn, y_jedn);
    kolo.close();
    //elipsa jednorodna
    ofstream elip;
    elip.open("elip.txt");
    afiniczna(elip, M_PI / 4., 1, 0.2, x_jedn, y_jedn);
    elip.close();
    //elipsa normalna
    ofstream elip_norm;
    elip_norm.open("elip_norm.txt");
    afiniczna(elip_norm, M_PI / 4., 1, 0.2, x_norm, y_norm);
    elip_norm.close();
    //kowariancja dla elipsy normalnej
    ofstream kowar1;
    kowar1.open("kowar1.txt");
    kowariancja(kowar1, x_norm, y_norm);
    //kowariancja dla elipsy jednorodnej
    ofstream kowar2;
    kowar2.open("kowar2.txt");
    kowariancja(kowar2, x_jedn, y_jedn);
    kowar2.close();
}