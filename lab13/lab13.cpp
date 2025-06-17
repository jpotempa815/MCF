#include "dsmc_2d.cpp"

void mergeFiles(const string& file1, const string& file2, const string& outputFile) {
    ifstream input1(file1);
    ifstream input2(file2);
    ofstream output(outputFile);

    string line;
    while (getline(input1, line)) {
        output << line << endl;
    }

    while (getline(input2, line)) {
        output << line << endl;
    }

    input1.close();
    input2.close();
    output.close();
}

int main() {
    DSMC_2D left;
    left.read("left.dat");
    left.init();
    left.write_position_velocity ("rv_left.dat");


    DSMC_2D right;
    right.read("right.dat");
    right.init();
    right.write_position_velocity ("rv_right.dat");

    mergeFiles("rv_left.dat", "rv_right.dat", "pos_vel_start.dat");

    DSMC_2D all;
    all.read("all.dat");
    all.init();
    all.nthreads = 7;
    all.icol=1;
    all.evolution(0.0, 2000);
    all.write_position_velocity("pos_vel_end.dat");

    return 0;
}