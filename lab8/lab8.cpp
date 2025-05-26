#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <algorithm>

using namespace std;

const double R0 = 1.315;
const double R1 = 1.70;
const double R2 = 2.00;
const double De = 6.325;
const double S = 1.29;
const double lambda = 1.5;
const double delta = 0.80469;
const double a0 = 0.011304;
const double c0 = 19.0;
const double d0 = 2.5;
const double PI = 3.14159265359;


struct Atom {
    double r, phi, theta;
    double x, y, z;
};

vector<Atom> atoms;
int n_atoms;

double beta_min = 1.0;
double beta_max = 100.0;
double p_exponent = 2.0;
int it_max = 100000;
double w_r = 1e-4;
double w_phi = 0.05;
double w_theta = 0.05;
double W_all = 1e-4;

double uniform() {
    return rand() / static_cast<double>(RAND_MAX);
}
double f_cut(double r) {
    if (r <= R1) return 1.0;
    if (r > R2) return 0.0;
    return 0.5 * (1.0 + cos(PI * (r - R1) / (R2 - R1)));
}


double V_R(double r) {
    return (De / (S - 1.0)) * exp(-sqrt(2.0 * S) * lambda * (r - R0));
}

double V_A(double r) {
    return (De * S / (S - 1.0)) * exp(-sqrt(2.0 / S) * lambda * (r - R0));
}

double g_theta(double cos_theta) {
    double c0_2 = c0 * c0;
    double d0_2 = d0 * d0;
    double denominator = d0_2 + (1.0 + cos_theta) * (1.0 + cos_theta);
    return a0 * (1.0 + c0_2/d0_2 - c0_2/denominator);
}

void spherical_to_cartesian(int i) {
    atoms[i].x = atoms[i].r * sin(atoms[i].theta) * cos(atoms[i].phi);
    atoms[i].y = atoms[i].r * sin(atoms[i].theta) * sin(atoms[i].phi);
    atoms[i].z = atoms[i].r * cos(atoms[i].theta);
}

double distance(int i, int j) {
    double dx = atoms[i].x - atoms[j].x;
    double dy = atoms[i].y - atoms[j].y;
    double dz = atoms[i].z - atoms[j].z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

double calculate_cos_angle(int i, int j, int k) {
    double r_ij_x = atoms[j].x - atoms[i].x;
    double r_ij_y = atoms[j].y - atoms[i].y;
    double r_ij_z = atoms[j].z - atoms[i].z;

    double r_ik_x = atoms[k].x - atoms[i].x;
    double r_ik_y = atoms[k].y - atoms[i].y;
    double r_ik_z = atoms[k].z - atoms[i].z;

    double dot_product = r_ij_x * r_ik_x + r_ij_y * r_ik_y + r_ij_z * r_ik_z;
    double mag_ij = sqrt(r_ij_x*r_ij_x + r_ij_y*r_ij_y + r_ij_z*r_ij_z);
    double mag_ik = sqrt(r_ik_x*r_ik_x + r_ik_y*r_ik_y + r_ik_z*r_ik_z);

    if (mag_ij < 1e-10 || mag_ik < 1e-10) return 0.0;

    double cos_val = dot_product / (mag_ij * mag_ik);
    // Ograniczenie do zakresu [-1, 1]
    if (cos_val > 1.0) cos_val = 1.0;
    if (cos_val < -1.0) cos_val = -1.0;

    return cos_val;
}

double calculate_zeta(int i, int j) {
    double zeta = 0.0;

    for (int k = 0; k < n_atoms; k++) {
        if (k == i || k == j) continue;

        double r_ik = distance(i, k);
        if (r_ik > R2) continue;

        double cos_theta_ijk = calculate_cos_angle(i, j, k);
        double g_val = g_theta(cos_theta_ijk);
        double fc_ik = f_cut(r_ik);

        zeta += fc_ik * g_val;
    }

    return zeta;
}

//modyfikacje dla punktu 4
double calculate_zeta_modified(int i, int j) {
    double zeta = 0.0;

    for (int k = 0; k < n_atoms; k++) {
        if (k == i || k == j) continue;

        double r_ik = distance(i, k);
        if (r_ik > R2) continue;

        double cos_theta_ijk = calculate_cos_angle(i, j, k);

        if (cos_theta_ijk > 0) {
            zeta = 10.0;
            break;
        }

        double g_val = g_theta(cos_theta_ijk);
        double fc_ik = f_cut(r_ik);

        zeta += fc_ik * g_val;
    }

    return zeta;
}

double calculate_atom_energy(int i, bool modified) {
    double energy = 0.0;

    for (int j = 0; j < n_atoms; j++) {
        if (i == j) continue;

        double r_ij = distance(i, j);
        if (r_ij > R2) continue;

        double fc_ij = f_cut(r_ij);
        if (fc_ij < 1e-10) continue;
        double zeta_ij, zeta_ji;
        //modyfikacja dla pkt 4
        if (modified == true) {
            zeta_ij = calculate_zeta_modified(i, j);
            zeta_ji = calculate_zeta_modified(j, i);
        }
        else {
            zeta_ij = calculate_zeta(i, j);
            zeta_ji = calculate_zeta(j, i);
        }


        double B_ij = pow(1.0 + zeta_ij, -delta);
        double B_ji = pow(1.0 + zeta_ji, -delta);
        double B = 0.5 * (B_ij + B_ji);

        double V_rep = V_R(r_ij);
        double V_att = V_A(r_ij);

        energy += fc_ij * (V_rep - B * V_att);
    }

    return 0.5 * energy;
}

double calculate_total_energy_test() {
    double total = 0.0;

    for (int i = 0; i < n_atoms; i++) {
        for (int j = i + 1; j < n_atoms; j++) {
            double r_ij = distance(i, j);
            if (r_ij > R2) continue;

            double fc_ij = f_cut(r_ij);
            if (fc_ij < 1e-10) continue;

            double zeta_ij = calculate_zeta(i, j);
            double zeta_ji = calculate_zeta(j, i);

            double B_ij = pow(1.0 + zeta_ij, -delta);
            double B_ji = pow(1.0 + zeta_ji, -delta);
            double B = 0.5 * (B_ij + B_ji);

            double V_rep = V_R(r_ij);
            double V_att = V_A(r_ij);

            total += fc_ij * (V_rep - B * V_att);
        }
    }

    return total;
}

double calculate_total_energy(bool modified) {
    double total = 0.0;

    for (int i = 0; i < n_atoms; i++) {
        for (int j = i + 1; j < n_atoms; j++) {
            double r_ij = distance(i, j);
            if (r_ij > R2) continue;

            double fc_ij = f_cut(r_ij);
            if (fc_ij < 1e-10) continue;
            double zeta_ij, zeta_ji;
            //modyfikacja dla pkt 4
            if (modified == true) {
                zeta_ij = calculate_zeta_modified(i, j);
                zeta_ji = calculate_zeta_modified(j, i);
            }
            else {
                zeta_ij = calculate_zeta(i, j);
                zeta_ji = calculate_zeta(j, i);
            }

            double B_ij = pow(1.0 + zeta_ij, -delta);
            double B_ji = pow(1.0 + zeta_ji, -delta);
            double B = 0.5 * (B_ij + B_ji);

            double V_rep = V_R(r_ij);
            double V_att = V_A(r_ij);

            total += fc_ij * (V_rep - B * V_att);
        }
    }

    return total;
}

void initialize_atoms(int n, double r_init) {
    n_atoms = n;
    atoms.resize(n_atoms);

    for (int i = 0; i < n_atoms; i++) {
        atoms[i].r = r_init;
        atoms[i].phi = uniform() * 2.0 * PI;
        atoms[i].theta = uniform() * PI;
        spherical_to_cartesian(i);
    }
}


bool load_positions_from_file(const string& filename, int n) {
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "Nie można otworzyć pliku: " << filename << endl;
        return false;
    }

    n_atoms = n;
    atoms.resize(n_atoms);

    for (int i = 0; i < n_atoms; i++) {
        if (!(file >> atoms[i].x >> atoms[i].y >> atoms[i].z)) {
            cout << "Błąd podczas wczytywania pozycji atomu " << i << endl;
            return false;
        }

        atoms[i].r = sqrt(atoms[i].x*atoms[i].x + atoms[i].y*atoms[i].y + atoms[i].z*atoms[i].z);
        atoms[i].phi = atan2(atoms[i].y, atoms[i].x);
        if (atoms[i].phi < 0) atoms[i].phi += 2.0 * PI;

        if (atoms[i].r > 1e-10) {
            atoms[i].theta = acos(atoms[i].z / atoms[i].r);
        } else {
            atoms[i].theta = 0.0;
        }
    }

    file.close();
    return true;
}

double calculate_average_radius() {
    double sum = 0.0;
    for (int i = 0; i < n_atoms; i++) {
        sum += atoms[i].r;
    }
    return sum / n_atoms;
}

void calculate_pcf(vector<double>& pcf, int M) {
    double r_sr = calculate_average_radius();
    double r_max = 2.5 * r_sr;
    double dr = r_max / M;
    double Omega = 4.0 * PI * r_sr * r_sr;

    fill(pcf.begin(), pcf.end(), 0.0);

    for (int i = 0; i < n_atoms; i++) {
        for (int j = i + 1; j < n_atoms; j++) {
            double r = distance(i, j);
            int m = static_cast<int>(floor(r / dr));
            if (m < M) {
                double r_m = (m + 0.5) * dr;
                if (r_m > 1e-10) {
                    pcf[m] += (2.0 * Omega) / (n_atoms * n_atoms) / (2.0 * PI * r_m * dr);
                }
            }
        }
    }
}

bool move_atom(int i, double beta, bool modified) {
    Atom old_atom = atoms[i];
    double old_energy = calculate_atom_energy(i, modified);

    double U1 = uniform();
    double U2 = uniform();
    double U3 = uniform();

    double dr = atoms[i].r * (2.0 * U1 - 1.0) * w_r;
    double dphi = (2.0 * U2 - 1.0) * w_phi;
    double dtheta = (2.0 * U3 - 1.0) * w_theta;

    atoms[i].r += dr;
    atoms[i].phi += dphi;
    atoms[i].theta += dtheta;

    if (atoms[i].phi < 0) atoms[i].phi += 2.0 * PI;
    if (atoms[i].phi > 2.0 * PI) atoms[i].phi -= 2.0 * PI;
    if (atoms[i].theta < 0) atoms[i].theta = old_atom.theta;
    if (atoms[i].theta > PI) atoms[i].theta = old_atom.theta;
    if (atoms[i].r < 0.5) atoms[i].r = old_atom.r; // Minimalna odległość od środka

    spherical_to_cartesian(i);

    double new_energy = calculate_atom_energy(i, modified);
    double delta_E = new_energy - old_energy;

    double p_acc = min(1.0, exp(-beta * delta_E));
    double U4 = uniform();

    if (U4 <= p_acc) {
        return true;
    } else {
        atoms[i] = old_atom;
        return false;
    }
}

bool global_radius_change(double beta, bool modified) {
    vector<double> old_radii(n_atoms);
    for (int i = 0; i < n_atoms; i++) {
        old_radii[i] = atoms[i].r;
    }

    double old_energy = calculate_total_energy(modified);

    double U1 = uniform();
    double scale_factor = 1.0 + W_all * (2.0 * U1 - 1.0);

    for (int i = 0; i < n_atoms; i++) {
        atoms[i].r *= scale_factor;
        spherical_to_cartesian(i);
    }

    double new_energy = calculate_total_energy(modified);
    double delta_E = new_energy - old_energy;

    double p_acc = min(1.0, exp(-beta * delta_E));
    double U2 = uniform();

    if (U2 <= p_acc) {
        return true;
    } else {
        for (int i = 0; i < n_atoms; i++) {
            atoms[i].r = old_radii[i];
            spherical_to_cartesian(i);
        }
        return false;
    }
}


void save_structure(const string& filename, int n_atoms) {
    ofstream file(filename);
    file << n_atoms << endl;
    for (int i = 0; i < n_atoms; i++) {
        file << fixed << setprecision(6)
             << "C " << atoms[i].x << " " << atoms[i].y << " " << atoms[i].z << "\n";
    }
    file.close();
}

void save_pcf(const string& filename, int M = 100) {
    vector<double> pcf(M, 0.0);
    calculate_pcf(pcf, M);

    double r_sr = calculate_average_radius();
    double r_max = 2.5 * r_sr;
    double dr = r_max / M;

    ofstream file(filename);
    for (int m = 0; m < M; m++) {
        double r_m = (m + 0.5) * dr;
        file << fixed << setprecision(6) << r_m << " " << pcf[m] << "\n";
    }
    file.close();
}

void test_c60_energy(bool modified) {
    cout << "=== Test poprawności dla C60 ===" << endl;
    load_positions_from_file("atoms_positions_c60.dat", 60);

    cout << "Obliczanie energii dla " << n_atoms << " atomów..." << endl;

    double total_energy = calculate_total_energy(modified);
    double binding_energy = total_energy / n_atoms;
    double avg_radius = calculate_average_radius();

    cout << "\n=== Wyniki testu ===" << endl;
    cout << "Całkowita energia: " << fixed << setprecision(3) << total_energy << " eV" << endl;
    cout << "Energia wiązania na atom: " << fixed << setprecision(3) << binding_energy << " eV" << endl;
    cout << "Średni promień: " << fixed << setprecision(3) << avg_radius << " A" << endl;
    cout << "\n=== Wartości oczekiwane ===" << endl;
    cout << "Oczekiwana energia całkowita: -421.6 eV" << endl;
    cout << "Oczekiwana energia wiązania: -7.027 eV" << endl;
    cout << "Oczekiwany średni promień: 3.52 A" << endl;

    save_structure("c60_test_structure.dat", n_atoms);
    cout << "\nZapisano strukturę testową do pliku: c60_test_structure.dat" << endl;
}


void run_simulation(int zad, bool modified) {
    cout << "\n=== Zadanie " << zad << ": Symulacja C60 ===" << endl;
    if ( zad == 5 ) initialize_atoms(60, 2.5);
    else initialize_atoms(60, 3.5);

    cout << "Parametry symulacji:" << endl;
    cout << "n = " << n_atoms << endl;
    cout << "beta_min = " << beta_min << ", beta_max = " << beta_max << endl;
    cout << "p = " << p_exponent << endl;
    cout << "it_max = " << it_max << endl;
    cout << "w_r = " << w_r << ", w_phi = " << w_phi << ", w_theta = " << w_theta << endl;
    cout << "W_all = " << W_all << endl;

    cout << "\nRozpoczynam symulację..." << endl;

    ofstream energy_file("energy" + to_string(zad) + ".dat");
    ofstream radius_file("radius" + to_string(zad) + ".dat");


    for (int it = 1; it <= it_max; it++) {
        double beta = beta_min + pow(static_cast<double>(it) / it_max, p_exponent) * (beta_max - beta_min);

        for (int i = 0; i < n_atoms; i++) {
            move_atom(i, beta, modified);
        }

        global_radius_change(beta, modified);

        if (it % 100 == 0) {
            double total_energy = calculate_total_energy(modified);
            double avg_radius = calculate_average_radius();

            energy_file << it << " " << fixed << setprecision(6) << total_energy << "\n";
            radius_file << it << " " << fixed << setprecision(6) << avg_radius << "\n";

            if (it % 1000 == 0) {
                cout << "Iteracja " << it << ": E = " << total_energy
                     << " eV, r_sr = " << avg_radius << " A" << endl;
            }
        }
    }

    energy_file.close();
    radius_file.close();

    double final_energy = calculate_total_energy(modified);
    double final_radius = calculate_average_radius();

    cout << "\n=== Wyniki końcowe ===" << endl;
    cout << "Końcowa energia całkowita: " << fixed << setprecision(3) << final_energy << " eV" << endl;
    cout << "Końcowa energia wiązania na atom: " << fixed << setprecision(3) << final_energy/n_atoms << " eV" << endl;
    cout << "Końcowy średni promień: " << fixed << setprecision(3) << final_radius << " A" << endl;

    save_structure("final_structure" + to_string(zad) + ".dat", n_atoms);
    save_pcf("pcf_histogram" + to_string(zad) + ".dat", 100);

}

struct SimulationParams {
    double beta_min, beta_max, p_exp, w_r_val, w_phi_val, w_theta_val;
    string description;
};

void test_parameters(bool modified) {
    cout << "\n=== Zadanie 6: Test różnych parametrów ===" << endl;

    vector<SimulationParams> param_sets = {
        {1.0, 100.0, 2.0, 1e-4, 0.05, 0.05, "Parametry bazowe"},
        {0.5, 50.0, 2.0, 1e-4, 0.05, 0.05, "Niższe beta"},
        {2.0, 200.0, 2.0, 1e-4, 0.05, 0.05, "Wyższe beta"},
        {1.0, 100.0, 1.0, 1e-4, 0.05, 0.05, "p=1"},
        {1.0, 100.0, 3.0, 1e-4, 0.05, 0.05, "p=3"},
        {1.0, 100.0, 2.0, 1e-3, 0.05, 0.05, "Większe w_r"},
        {1.0, 100.0, 2.0, 1e-5, 0.05, 0.05, "Mniejsze w_r"},
        {1.0, 100.0, 2.0, 1e-4, 0.1, 0.1, "Większe w_phi, w_theta"},
        {1.0, 100.0, 2.0, 1e-4, 0.01, 0.01, "Mniejsze w_phi, w_theta"}
    };

    ofstream results_file("parameter_test_results.dat");

    for (size_t test_nr = 0; test_nr < param_sets.size(); test_nr++) {
        cout << "\nTest " << test_nr + 1 << ": " << param_sets[test_nr].description << endl;

        double orig_beta_min = beta_min;
        double orig_beta_max = beta_max;
        double orig_p_exponent = p_exponent;
        double orig_w_r = w_r;
        double orig_w_phi = w_phi;
        double orig_w_theta = w_theta;

        beta_min = param_sets[test_nr].beta_min;
        beta_max = param_sets[test_nr].beta_max;
        p_exponent = param_sets[test_nr].p_exp;
        w_r = param_sets[test_nr].w_r_val;
        w_phi = param_sets[test_nr].w_phi_val;
        w_theta = param_sets[test_nr].w_theta_val;

        initialize_atoms(60, 2.5);

        int orig_it_max = it_max;
        it_max = 50000; 

        for (int it = 1; it <= it_max; it++) {
            double beta = beta_min + pow(static_cast<double>(it) / it_max, p_exponent) * (beta_max - beta_min);

            for (int i = 0; i < n_atoms; i++) {
                move_atom(i, beta, modified);
            }
            global_radius_change(beta, modified);
            if (it%1000==0) cout << it << endl;
        }

        double final_energy = calculate_total_energy(modified);
        double binding_energy = final_energy / n_atoms;
        double final_radius = calculate_average_radius();

        cout << "Energia wiązania: " << fixed << setprecision(3) << binding_energy << " eV" << endl;
        cout << "Średni promień: " << fixed << setprecision(3) << final_radius << " A" << endl;

        results_file << test_nr + 1 << " \"" << param_sets[test_nr].description << "\" "
                     << fixed << setprecision(6) << final_energy << " "
                     << binding_energy << " " << final_radius << "\n";

        beta_min = orig_beta_min;
        beta_max = orig_beta_max;
        p_exponent = orig_p_exponent;
        w_r = orig_w_r;
        w_phi = orig_w_phi;
        w_theta = orig_w_theta;
        it_max = orig_it_max;
    }

    results_file.close();
    cout << "\nWyniki testów parametrów zapisano do: parameter_test_results.dat" << endl;
}

void run_multiple_n(bool modified) {
    cout << "\n=== Zadanie 7: Seria symulacji dla n=30-40 ===" << endl;

    ofstream stability_file("stability_series.dat");

    for (int n = 20; n <= 100; n+=10) {
        cout << "\nSymulacja dla n = " << n << endl;

        initialize_atoms(n, 2.5);

        for (int it = 1; it <= it_max; it++) {
            double beta = beta_min + pow(static_cast<double>(it) / it_max, p_exponent) * (beta_max - beta_min);

            for (int i = 0; i < n_atoms; i++) {
                move_atom(i, beta, modified);
            }
            global_radius_change(beta, modified);

            if (it % 10000 == 0) {
                double total_energy = calculate_total_energy(modified);
                double avg_radius = calculate_average_radius();
                cout << "  Iteracja " << it << ": E = " << total_energy
                     << " eV, r_sr = " << avg_radius << " A" << endl;
            }
        }

        double final_energy = calculate_total_energy(modified);
        double binding_energy = final_energy / n_atoms;
        double final_radius = calculate_average_radius();

        cout << "n=" << n << ": E_b = " << fixed << setprecision(3) << binding_energy
             << " eV, r_sr = " << final_radius << " A" << endl;

        stability_file << n << " " << fixed << setprecision(6) << final_energy << " "
                       << binding_energy << " " << final_radius << "\n";

        // Zapisz strukturę dla każdego n
        string filename = "structure_n" + to_string(n) + ".dat";
        save_structure(filename, n_atoms);
    }

    stability_file.close();
    cout << "\nWyniki serii stabilności zapisano do: stability_series.dat" << endl;
}

void run_all_tasks() {
    cout << "Symulacja Monte Carlo dla fullerenow" << endl;
    cout << "====================================" << endl;

    // Zadanie 2: Test poprawności dla C60
    test_c60_energy(false);

    // Zadanie 3: Symulacja SA
    run_simulation(3, false);

    // Zadanie 4: Zmodyfikowany potencjał
    run_simulation(4, true);

    // Zadanie 5: Start od r=2.5A
    run_simulation(5, true);

    // Zadanie 6: Test parametrów
    test_parameters(true);

    // Zadanie 7: Seria stabilności
    run_multiple_n(true);
}

int main() {

    run_all_tasks();

    return 0;
}