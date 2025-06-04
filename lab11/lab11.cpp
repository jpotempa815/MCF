#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include"photon_diffusion.cpp"

using namespace std;

void save_data(const string& filename, const vector<vector<double>>& data) {
    ofstream file(filename);
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[i].size(); j++) {
            file << data[i][j] << " ";
        }
        file << endl;
    }
    file.close();
}

void save_data(const string& filename, const vector<double>& data) {
    ofstream file(filename);
    for (double val : data) {
        file << val << " ";
    }
    file << endl;
    file.close();
}

double normalization_check(const vector<vector<double>>& absorbance, const vector<double>& transmitance, const vector<double>& reflectance) {
    double norm = 0;
    for (int i = 0; i < absorbance.size(); i++) {
        for (int j = 0; j < absorbance[i].size(); j++) {
            norm += absorbance[i][j];
        }
    }
    for (int i = 0; i < transmitance.size(); i++) {
        norm += transmitance[i];
    }
    for (int i = 0; i < reflectance.size(); i++) {
        norm += reflectance[i];
    }
    return norm;
}


int main(){  
    for (int i = 1; i <= 8; i++) {
        PHOTON_DIFFUSION_2D ob;

        ob.xmax = 0.2;
        ob.x_source = 0.1;
        ob.dx_source = 0.0;
        ob.x_detect = 0.15;
        ob.dx_detect = 0.01;
        ob.nx = 100;
        ob.ny = 100;

        ob.rx0 = 0.0;
        ob.ry0 = 1.0;

        ob.nlayers = 3;

        ob.layers_data[1][0] = 1.;
        ob.layers_data[1][1] = 10;
        ob.layers_data[1][2] = 0.02;
        ob.layers_data[1][3] = 0.75;
        ob.layers_data[1][4] = 1.3;

        ob.layers_data[2][0] = 1.0;
        ob.layers_data[2][1] = 190;
        ob.layers_data[2][2] = 0.02;
        ob.layers_data[2][3] = 0.075;
        ob.layers_data[2][4] = 1.0;

        ob.layers_data[3][0] = 10.0;
        ob.layers_data[3][1] = 90;
        ob.layers_data[3][2] = 0.02;
        ob.layers_data[3][3] = 0.95;
        ob.layers_data[3][4] = 1.0;

        switch (i) {
            case 1:
                //ppkt 2: 1
                ob.rx0 = 0.8;
                ob.ry0 = 0.6;
                ob.layers_data[2][4] = 1.5;
                break;
            case 2:
                //ppkt 2: 2
                ob.layers_data[2][4] = 2.5;
                ob.rx0 = 0.8;
                ob.ry0 = 0.6;
                break;
            case 3:
                //ppkt 2: 3
                ob.layers_data[1][4] = 1.0;
                ob.rx0 = 0.8;
                ob.ry0 = 0.6;
                ob.layers_data[2][4] = 1.5;
                break;
            case 4:
                //ppkt 2: 4
                ob.layers_data[2][1] = 10.;
                ob.layers_data[2][4] = 1.5;
                ob.rx0 = 0.8;
                ob.ry0 = 0.6;
                ob.layers_data[1][4] = 1.0;
                break;
            case 5:
                break;
            case 6:
                //ppkt 3: 2
                ob.layers_data[1][4] = 1.0;
                ob.layers_data[2][0] = 10.;
                ob.layers_data[2][1] = 210.;
                ob.layers_data[2][4] = 1.5;
                break;
            case 7:
                //ppkt 3: 3
                ob.layers_data[1][4] = 1.0;
                ob.layers_data[2][0] = 1.;
                ob.layers_data[2][1] = 210.;
                ob.layers_data[2][4] = 1.5;
                break;
            case 8:
                //ppkt 3: 4
                ob.layers_data[1][4] = 1.0;
                ob.layers_data[2][0] = 10.;
                ob.layers_data[2][1] = 210.;
                ob.layers_data[2][4] = 1.5;
                ob.layers_data[2][3] = 0.75;
                break;
            default:
                break;
        }
        ob.init();

        int N = 200000;
        ob.write_all_paths = 0;
        ob.write_source_detection_paths = 0;

        for (int k = 0; k < N; k++) {
            ob.single_path();
        }
        save_data("abs_" + to_string(i) + ".dat", ob.absorption);
        save_data("trans_" + to_string(i) + ".dat", ob.transmittance);
        save_data("refl_" + to_string(i) + ".dat", ob.reflectance);
        cout << "Normalization check for " << i << ": " << normalization_check(ob.absorption, ob.transmittance, ob.reflectance)/static_cast<double>(N) << endl;
    }
}