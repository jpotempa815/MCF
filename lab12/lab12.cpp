#include "dsmc_2d.cpp"

using namespace std;

int main() {
    DSMC_2D ob;

    // ob.read("zad1_r1.dat");
    // ob.init();
    // ob.write_position_velocity("rv_ini_1_r1.dat");
    // ob.nthreads = 7;
    // ob.icol=1;
    // ob.evolution(0.0, 10000, "1r1");
    // ob.hist_velocity_all("hist1.dat",5.0, 50);
    // ob.write_position_velocity("rv_end_1_r1.dat");

    // DSMC_2D ob2;
    // ob2.read("zad1_r2.dat");
    // ob2.init();
    // ob2.write_position_velocity("rv_ini_1_r2.dat");
    // ob2.nthreads = 7;
    // ob2.icol=1;
    // ob2.evolution(0.0, 10000, "1r2");
    // ob2.hist_velocity_all("hist2.dat",5.0, 50);
    // ob2.write_position_velocity("rv_end_1_r2.dat");

    DSMC_2D ob3;
    ob3.read("zad2_r1.dat");
    ob3.init();
    ob3.write_position_velocity("rv_ini_2_r1.dat");
    ob3.nthreads = 7;
    ob3.icol=1;
    ob3.evolution(0.0, 10000, "2r1");
    ob3.hist_velocity_all("hist3.dat",5.0, 50);
    ob3.write_position_velocity("rv_end_2_r1.dat");

    // DSMC_2D ob4;
    // ob4.read("zad3_r1.dat");
    // ob4.init();
    // ob4.write_position_velocity("rv_ini_3_r1.dat");
    // ob4.nthreads = 7;
    // ob4.icol=1;
    // ob4.evolution(0.0, 10000, "3r1");
    // ob4.hist_velocity_all("hist5.dat",5.0, 50);
    // ob4.write_position_velocity("rv_end_3_r1.dat");
    //
    //
    // DSMC_2D ob5;
    // ob5.read("zad4_r1.dat");
    // ob5.init();
    // ob5.write_position_velocity("rv_ini_4_r1.dat");
    // ob5.nthreads = 7;
    // ob5.icol=1;
    // ob5.evolution(0.0, 10000, "4r1");
    // ob5.hist_velocity_all("hist7.dat",5.0, 50);
    // ob5.write_position_velocity("rv_end_4_r1.dat");

    return 0;
}