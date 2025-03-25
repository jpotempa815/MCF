#include<iostream>
#include<cmath>
#include<fstream>
#include <vector>

using namespace std;

int N {1000000};
//liczba przedzialow
int k {10};

double uniform(){
    return rand()/(double)RAND_MAX;
}

vector<double> zlozony(){
    double g_1 = 4/5.;
    double U_1 {0};
    double U_2 {0};
    double x {0};
    vector<double> X(N);

    for(int i = 0; i < N; i++)
    {
        U_1 = uniform();
        U_2 = uniform();
        if(U_1 <= g_1) x = U_2;
        else x = sqrt(1-sqrt(1-U_2));
        X[i] = x;
    }
    return X;
}

double f(double x)
{
    return 4/5. * (1 + x - x*x*x);
}

double F(double x)
{
    return 4/5. * (x + x*x/2. - x*x*x*x/4.);
}

vector<double> lancuch(double delta){
    // wartosci lancucha
    vector<double> X(N);
    X[0] = uniform();
    //inicjalizacja wartosci startowych
    double U_1 {0};
    double U_2 {0};
    // prawdopodobienstwo akceptacji, w którym skraca się czynnik T(X_i|X_i+1)/T(X_i+1|X_i) dla delta \in [0,1]
    double p_acc {0};
    //inicjalizacja nowego elementu lancucha
    double x_new {0};
    for(int i = 1; i < N; i++)
    {
        //losujemy wartosci startowe
        U_1 = uniform();
        U_2 = uniform();
        x_new = X[i-1] + (2* U_1 - 1)*delta;
        if(x_new <= 1 && x_new >= 0){
            p_acc = min((double)f(x_new)/f(X[i-1]), (double)1.);
            if(U_2 <= p_acc) X[i] = x_new;
            else X[i] = X[i-1];
        }
        else X[i] = X[i-1];
    }
    return X;
}

double eliminacja(double G){
    double U1, U2, G2;
    vector<double> X(N);
    do{
        U1 = uniform();
        U2 = uniform();
        G2 = G*U2;
    }while ((G2 > ( (double) 4/5*(1 + U1 - pow(U1,3)) ) ));

    return U1;
}

void zapis_do_pliku(vector<double> X, ostream& plik)
{
    for(int i = 0; i < N; i++) plik << X[i] << endl;
    for (int i = 0; i < N; i++) cout << i << endl;
}


int main(){
    vector<double> rozkl_zlozony, rozkl_lancuch_1, rozkl_lancuch_2, rozkl_eliminacja(N);
    rozkl_zlozony = zlozony();
    rozkl_lancuch_1 = lancuch(0.5);
    rozkl_lancuch_2 = lancuch(0.05);
    double G {1.15};
    for (int i=0; i < N; i++) {
        cout<< i <<endl;
        rozkl_eliminacja[i] = eliminacja(G);
    }

    //rozklad zlozony
    ofstream plik_1;
    plik_1.open("zlozony.txt");

    zapis_do_pliku(rozkl_zlozony, plik_1);

    plik_1.close();

    //lancuch Markowa
    ofstream lancuch_1;
    ofstream lancuch_2;
    lancuch_1.open("lancuch_1.txt");
    lancuch_2.open("lancuch_2.txt");

    zapis_do_pliku(rozkl_lancuch_1, lancuch_1);
    zapis_do_pliku(rozkl_lancuch_2, lancuch_2);

    lancuch_1.close();
    lancuch_2.close();

    //metoda eliminacji
    ofstream plik_4;
    plik_4.open("eliminacja.txt");

    zapis_do_pliku(rozkl_eliminacja, plik_4);

    plik_4.close();

    return 0;
}
