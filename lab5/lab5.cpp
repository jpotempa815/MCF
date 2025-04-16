#include<iostream>
#include<cmath>
#include<cstdlib>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

double uniform(){
    return rand() / (double) RAND_MAX;
}

double gx(double x, string& calka){
  if(calka == "tanh"){
      return 1. + tanh(x);
  }
  else if(calka == "1/x"){
      return 1.0 / (1 + x*x);
  }
  else if(calka == "cos"){
      return pow(cos(M_PI * x), 10);
  }
    return 0.0;
}

void podstawowa(int N, double a, double b, string& calka, ofstream& nazwa_pliku, ofstream& stat){
    double g_sr = 0, g_x = 0, g_sr2 = 0, var = 0, x = 0, err = 0;
    for(int i = 1; i <= N; i++){
        x = a + (b-a) * uniform();
        g_x = gx(x, calka);
        nazwa_pliku << x << endl;
        g_sr += (b-a) * g_x;
        g_sr2 += (b-a) * (b-a) * g_x * g_x;
        if (i == 100 || i == 1000 || i == 10000 || i == 100000) {
            double g_sr2_i = g_sr2/static_cast<double>(i);
            double g_sr_i = g_sr/static_cast<double>(i);
            var = (g_sr2_i - g_sr_i*g_sr_i)/static_cast<double>(i);
            err = sqrt(var)/(g_sr_i) * 100;
            cout << i << " " << g_sr_i << " " << sqrt(var) << " " << err << "%" << endl;
            stat << i << "," << g_sr_i << "," << sqrt(var) << "," << err << endl;
        }
    }
}

void los_system(int N, double a, double b, string& calka, ofstream& nazwa_pliku, ofstream& stat, bool opt, vector<double>& var_est_1, vector<double>& var_est_2){
    const int M = 10;
    int Nm = 0;
    double dx = (b - a) / M;
    double pm = 1.0 / static_cast<double>(M);
    double C = 0.0, var = 0.0, err = 0.0;

    vector<int> Nm_vec(M, 0);

    if (opt) {
        vector<double> var_est = (N < 1000) ? var_est_1 : var_est_2;
        double var_sum = 0.0;
        for (int i = 0; i < M; ++i) {
            var_sum += pm * var_est[i];
        }

        for (int i = 0; i < M; ++i) {
            Nm_vec[i] = ceil((pm * var_est[i] / var_sum) * N);
        }
    }

    vector<double> g_sr_vec(M, 0.0);
    vector<double> g_sr2_vec(M, 0.0);
    vector<double> var_m_vec(M, 0.0);
    vector<int> actual_Nm(M, 0);

    for (int i = 0; i < M; ++i) {
        double x_l = a + i * dx;
        double x_p = x_l + dx;
        if (opt) Nm = Nm_vec[i];
        else Nm = pm * N;

        actual_Nm[i] = Nm;
        double g_sr = 0.0, g_sr2 = 0.0;
        for (int j = 0; j < Nm; j++) {
            double x = x_l + (x_p - x_l) * uniform();
            nazwa_pliku << x << endl;
            double g_x = gx(x, calka);
            g_sr_vec[i] += g_x;
            g_sr2_vec[i] += g_x * g_x;
        }
        g_sr_vec[i] = ((b-a) * g_sr_vec[i]) / static_cast<double>(Nm);
        g_sr2_vec[i] = ((b-a) * (b-a) * g_sr2_vec[i]) / static_cast<double>(Nm);
        var_m_vec[i] = g_sr2_vec[i] - g_sr_vec[i] * g_sr_vec[i];
        // C += pm * g_sr;
        C += pm * g_sr_vec[i];
        var += (pm * pm / static_cast<double>(Nm)) * var_m_vec[i];

        if (!opt && N == 100) {
            if (i < var_est_1.size()) {
                var_est_1[i] = sqrt(var_m_vec[i]);
            } else {
                var_est_1.push_back(sqrt(var_m_vec[i]));
            }
        }
        else if (!opt && N == 1000){
            if (i < var_est_2.size()) {
                var_est_2[i] = sqrt(var_m_vec[i]);
            } else {
                var_est_2.push_back(sqrt(var_m_vec[i]));
            }
        }
    }

    err = sqrt(var) / C * 100;

    cout << N << " " << C << " " << sqrt(var) << " " << err << "%" << endl;
    stat << N << "," << C << "," << sqrt(var) << "," << err << endl;

}



int main(){
    ofstream all_stat;
    all_stat.open("stat.csv");

    // tanh
    string calka {"tanh"};
    ofstream tanh;
    // tanh.open("tanh.txt");
    tanh.open("tanh.txt");
    podstawowa(100000, -3., 3., calka, tanh, all_stat);
    tanh.close();
    vector<double> var_est_1(10,0.0), var_est_2(10, 0.0);
    int N {100000};
    ofstream los_tanh;
    los_tanh.open("tanh_los.txt");
    for(int i = 100; i <= N; i*=10){
        los_system(i, -3., 3., calka, los_tanh, all_stat, false, var_est_1, var_est_2);
    }
    los_tanh.close();
    ofstream opt_tanh;
    opt_tanh.open("tanh_opt.txt");
    for (int i = 100; i <= N; i*=10) {
        los_system(i, -3., 3., calka, opt_tanh, all_stat, true, var_est_1, var_est_2);
    }
    opt_tanh.close();

    // 1/x
    calka = "1/x";
    ofstream x_2;
    x_2.open("x_2.txt");
    podstawowa(100000, 0, 10., calka, x_2, all_stat);
    x_2.close();
    var_est_1.clear();
    var_est_2.clear();
    ofstream los_x_2;
    los_x_2.open("x_2_los.txt");
    for(int i = 100; i <= N; i*=10){
        los_system(i, 0, 10., calka, los_x_2, all_stat, false, var_est_1, var_est_2);
    }
    los_x_2.close();
    ofstream opt_x_2;
    opt_x_2.open("x_2_opt.txt");
    for (int i = 100; i <= N; i*=10) {
        los_system(i, 0, 10., calka, opt_x_2, all_stat, true, var_est_1, var_est_2);
    }
    opt_x_2.close();

    // cos
    calka = "cos";
    ofstream cos;
    cos.open("cos.txt");
    podstawowa(100000, 0, 1., calka, cos, all_stat);
    cos.close();
    var_est_1.clear();
    var_est_2.clear();
    ofstream los_cos;
    los_cos.open("cos_los.txt");
    for(int i = 100; i <= N; i*=10){
        los_system(i, 0., 1., calka, los_cos, all_stat, false, var_est_1, var_est_2);
    }
    los_cos.close();
    ofstream opt_cos;
    opt_cos.open("cos_opt.txt");
    for (int i = 100; i <= N; i*=10) {
        los_system(i, 0., 1., calka, opt_cos, all_stat, true, var_est_1, var_est_2);
    }
    opt_cos.close();

    all_stat.close();
}