#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cstdlib>
#include<cmath>

using namespace std;

int nx = 30, ny = 30;
double d = 0.1, Vl = 1., Vt = -1., Vb = -1., eps = 1, rho_max = 1.;
double x_max = nx * d, y_max = ny * d;
double sigma_p = x_max/10.;

double uniform() {
    return rand()/static_cast<double>(RAND_MAX);
}

//relaksacja
vector<vector<double>> relaxation_method() {
    int it_max = 1e4;
    double omega = 1.8, tol = 1e-6;
    double F_old = 0.0, F_new = 0.0;
    vector<vector<double>> V (nx+1, vector<double>(ny+1, 0.0));
    vector<vector<double>> rho (nx+1, vector<double>(ny+1, 0.0));

    //inicjalizacja rho
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            rho[i][j] = rho_max * exp(-((i * d - x_max / 2) * (i * d - x_max / 2) + (j * d - y_max / 2) * (j * d - y_max / 2)) / (2 * sigma_p * sigma_p));
        }
    }

    //WB Dirichleta
    for (int j = 0; j <= ny; j++) V[0][j] = Vl * sin(M_PI*j*d/(y_max));
    for (int i = 0; i <= nx; ++i) {
        V[i][0] = Vb * sin(M_PI* i*d/(x_max));
        V[i][ny] = Vt * sin(M_PI* i*d/(x_max));
    }

    for (int it = 1; it < it_max; it++) {
        for (int i = 1; i < nx; i++) {
            for (int j = 1; j < ny; j++) {
                V[i][j] = (1- omega)*V[i][j] + omega/4.0 * (V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + (d * d * rho[i][j]) / eps);
            }
        }
        //WB von Neumanna
        for (int j = 1; j < ny; j++) V[nx][j] = V[nx-1][j];

        //warunek zbieznosci
        F_old = F_new;
        F_new = 0.0;
        for (int i = 1; i < nx; i++) {
            for (int j = 1; j < ny; j++) {
                double Ex = (V[i+1][j] - V[i-1][j]) / (2 * d);
                double Ey = (V[i][j+1] - V[i][j-1]) / (2 * d);
                F_new += (Ex*Ex + Ey*Ey)/2 - rho[i][j]*V[i][j];
            }
        }
        if (fabs((F_new - F_old)/F_new) < tol) break;
    }
    return V;
}

vector<vector<double>> montecarlo(int N_chains, int n_length, bool block, int typ) {
    vector<vector<double>> V(nx + 1, vector<double>(ny + 1, 0.0));
    //odchylenie standardowe
    vector<vector<double>> sigma_V(nx + 1, vector<double>(ny + 1, 0.0));
    vector<vector<double>> rho(nx + 1, vector<double>(ny + 1, 0.0));
    //wskaznik WB Dirichleta
    vector<vector<double>> B(nx + 1, vector<double>(ny + 1, 0.0));
    //tablica lancuchow zakonczonych absorpcja
    vector<vector<double>> S(nx + 1, vector<double>(ny + 1, 0.0));

    //pliki do zapisu odchylenia standardowego i ulamka zaaabsorbowanych lancuchow
    ofstream sigma("sigma_" + to_string(typ) + ".dat");
    ofstream S_file("S_" + to_string(typ) + ".dat");
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            rho[i][j] = rho_max * exp(-((i * d - x_max / 2) * (i * d - x_max / 2) + (j * d - y_max / 2) * (j * d - y_max / 2)) / (2 * sigma_p * sigma_p));
        }
    }
    //WB Dirichleta
    for (int j = 0; j <= ny; j++) {
        V[0][j] = Vl * sin(M_PI * j * d / (y_max));
        B[0][j] = 1.0;
    }
    for (int i = 0; i <= nx; i++) {
        V[i][0] = Vb * sin(M_PI * i * d / (x_max));
        B[i][0] = 1.0;
        V[i][ny] = Vt * sin(M_PI * i * d / (x_max));
        B[i][ny] = 1.0;
    }
    for (int i_0 = 1; i_0 < nx; i_0++) {
        for (int j_0 = 1; j_0 < ny; j_0++) {
            double sum_V1 = 0.0, sum_V2 = 0.0;
            int k_chains = 0;
            for (int N = 1; N <= N_chains; N++ ) {
                int i = i_0;
                int j = j_0;
                double g = 0.0; // gestosc ladunkow
                for (int n = 1; n <= n_length; n++) {
                    int m = floor(4*uniform());
                    if (m == 0) i--;
                    else if (m == 1) i++;
                    else if (m == 2) j--;
                    else if (m == 3) j++;
                    //odbicie
                    if (i == nx+1) i = nx-1;
                    //absorpcja
                    if (B[i][j] == 1) {
                        double dV = V[i][j] + g;
                        sum_V1 += dV;
                        sum_V2 += dV * dV;
                        k_chains++;
                        break;
                    }
                    //brak absorpcji
                    g += rho[i][j] * d * d/(4*eps);
                }
            }
            double V1 = sum_V1 / static_cast<double>(k_chains);
            double V2 = sum_V2 / static_cast<double>(k_chains);
            V[i_0][j_0] = V1;
            sigma_V[i_0][j_0] = sqrt((V2 - V1 * V1)/static_cast<double>(k_chains));
            //aktualizacja wskaznika WB Dirichleta
            if (block == false) B[i_0][j_0] = 0.0;
            else B[i_0][j_0] = 1.0;
            S[i_0][j_0] = static_cast<double>(k_chains)/static_cast<double>(N_chains);
        }
    }
    //zapis do pliku
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            sigma << sigma_V[i][j] << " ";
            S_file << S[i][j] << " ";
        }
        sigma << endl;
        S_file << endl;
    }
    sigma.close();
    S_file.close();
    return V;
}

int main() {
    vector<vector<double>> V_rel = relaxation_method();
    vector<vector<double>> V_mc1 = montecarlo(100, 100, false, 1);
    vector<vector<double>> V_mc2 = montecarlo(100, 100, true, 2);
    vector<vector<double>> V_mc3 = montecarlo(300, 300, true, 3);
    ofstream pot_mc1("pot_mc1.dat");
    ofstream pot_mc2("pot_mc2.dat");
    ofstream pot_mc3("pot_mc3.dat");
    ofstream pot_rel1("pot_rel1.dat");
    ofstream pot_rel2("pot_rel2.dat");
    ofstream pot_rel3("pot_rel3.dat");
    for (int i = 0; i < V_rel.size(); i++) {
        for (int j = 0; j < V_rel[i].size(); j++) {
            pot_mc1 << V_mc1[i][j] << " ";
            pot_mc2 << V_mc2[i][j] << " ";
            pot_mc3 << V_mc3[i][j] << " ";
            pot_rel1 << fabs(V_rel[i][j]-V_mc1[i][j]) << " ";
            pot_rel2 << fabs(V_rel[i][j]-V_mc2[i][j]) << " ";
            pot_rel3 << fabs(V_rel[i][j]-V_mc3[i][j]) << " ";
        }
        pot_mc1 << endl;
        pot_mc2 << endl;
        pot_mc3 << endl;
        pot_rel1 << endl;
        pot_rel2 << endl;
        pot_rel3 << endl;
    }
    pot_mc1.close();
    pot_mc2.close();
    pot_mc3.close();
    pot_rel1.close();
    pot_rel2.close();
    pot_rel3.close();
}