#include<iostream>
#include<cstdlib>
#include<fstream>
#include<cmath>

using namespace std;

double uniform()
{
    return (rand()/(double)RAND_MAX);
}

void bernuli(int N, float p, int k, ofstream& bernouli)
{
    double U1 {0};
    float sum_X {0}, sum_X2 {0};
    //statystyka
    float X_sr {0}, X2_sr {0}, err_X {0}, var_num {0}, var_teo {0}, err_var{0};
    int save = 100;

    for(int i=1; i<=N; i++)
    {
        int X;
        U1 = uniform();
        if (U1 <= p) X = 1;
        else X = 0;
        sum_X = sum_X + X;
        sum_X2 = sum_X2 + X*X;

        if(i == save)
        {
            // k++;
            save *= 10;
            X_sr = sum_X/(float)i;
            X2_sr = sum_X2/(float)i;
            err_X = abs((X_sr-p)/p);
            var_num = (X2_sr - X_sr*X_sr)/(float)i;
            var_teo = (p - p*p)/(float)i;
            err_var = abs((var_num-var_teo)/var_teo);
            bernouli << i << " " << X_sr << " " << err_X << " " << err_var << endl;
        }

    }
}

int main(){
    int N = 1e7;
    int k = 2;
    
    float p1 {0.1}, p2 {0.5}, p3 {0.9};
    
    ofstream bernouli;
    bernouli.open("lab1.txt");
    
    bernuli(N, p1, k, bernouli);
    bernuli(N, p2, k, bernouli);
    bernuli(N, p3, k, bernouli);
    
    bernouli.close();
    return 0;
}