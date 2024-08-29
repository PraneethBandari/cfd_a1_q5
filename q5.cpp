#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int main()
{
    double L = 0.5;          // Wall thickness (m)
    int N = 50;              // Number of nodes
    double dx = L / (N - 1); // Step size (m)
    double T_left = 100;
    double T_right = 500;
    double accuracy = 1e-6;

    vector<double> x(N);
    for (int i = 0; i < N; i++) {
        x[i] = dx * i;
    }

    // Case 1
    double B = 0.001;
    cout << "Case 1\n";
    
    double C2 = 0.1 * ((T_left + 273.15) + B * (T_left + 273.15) * (T_left + 273.15) * 0.5); // Derived from B.c at x=0
    double C1 = (0.1 * ((T_right + 273.15) + B * (T_right + 273.15) * (T_right + 273.15) * 0.5) - C2) / L;
    vector<double> T_analyticalcase1(N);
    for (int i = 0; i < N; i++) {
        double T = (-1 + sqrt(1 + 2 * B * (C1 * x[i] + C2) * 10)) / B;
        T_analyticalcase1[i] = T - 273.15;
        cout << "Analytical solution: " << T_analyticalcase1[i] << endl;
    }

    vector<double> Tempcase1(N, (T_left + T_right) * 0.5 + 273.15);
    vector<double> K(N, 0.1 * (1 + B * Tempcase1[0]));
    vector<double> T_newcase1(N, 0);
    T_newcase1[0] = T_left + 273.15;
    T_newcase1[N - 1] = T_right + 273.15;

    double max_dt = 1;
    int n = 0; // Reset iteration count for case 1
    while (max_dt > accuracy) {
        max_dt = 0;
        for (int i = 1; i < N - 1; i++) {
            T_newcase1[i] = (K[i - 1] * Tempcase1[i - 1] + K[i + 1] * Tempcase1[i + 1]) / (K[i + 1] + K[i - 1]);
        }
        // Check for convergence and update temperatures
        for (int i = 0; i < N; i++) {
            double dt = fabs(Tempcase1[i] - T_newcase1[i]);
            if (dt > max_dt) {
                max_dt = dt;
            }
            Tempcase1[i] = T_newcase1[i];
            K[i] = 0.1 * (1 + B * Tempcase1[i]);
        }
        n++;
    }

    cout << "Number of iterations: " << n << "\n";
    for (int i = 0; i < N; i++) {
        cout << "Node " << i + 1 << ": Analytical temp = " << T_analyticalcase1[i] << " FDM temp = " << Tempcase1[i] - 273.15 << " \n";
    }

    ofstream myfile("fdmq5case1sol.txt");
    if (myfile.is_open()) {
        for (int i = 0; i < N; i++) {
            myfile << i + 1 << ":" << Tempcase1[i] - 273.15 << endl;
        }
        myfile.close();
        cout << "Temperature values have been written to fdmq5case1sol.txt" << endl;
    } else {
        cout << "Error opening file for FDM solution!" << endl;
    }

    ofstream outfile("fdmq5case1analyticalsol.txt");
    if (outfile.is_open()) {
        for (int i = 0; i < N; i++) {
            outfile << i + 1 << ":" << T_analyticalcase1[i] << endl;
        }
        outfile.close();
        cout << "Temperature values have been written to fdmq5case1analyticalsol.txt" << endl;
    } else {
        cout << "Error opening file for Analytical solution!" << endl;
    }

    // Case 2
    B = -0.001;
    cout << "Case 2\n";
    
    C2 = 0.1 * ((T_left + 273.15) + B * (T_left + 273.15) * (T_left + 273.15) * 0.5); // Derived from B.c at x=0
    C1 = (0.1 * ((T_right + 273.15) + B * (T_right + 273.15) * (T_right + 273.15) * 0.5) - C2) / L;
    vector<double> T_analyticalcase2(N);
    for (int i = 0; i < N; i++) {
        double T = (-1 + sqrt(1 + 2 * B * (C1 * x[i] + C2) * 10)) / B;
        T_analyticalcase2[i] = T - 273.15;
        cout << "Analytical solution: " << T_analyticalcase2[i] << endl;
    }

    vector<double> Tempcase2(N, (T_left + T_right) * 0.5 + 273.15);
    vector<double> Kcase2(N, 0.1 * (1 + B * Tempcase2[0]));
    vector<double> T_newcase2(N, 0);
    T_newcase2[0] = T_left + 273.15;
    T_newcase2[N - 1] = T_right + 273.15;

    max_dt = 1;
    n = 0; // Reset iteration count for case 2
    while (max_dt > accuracy) {
        max_dt = 0;
        for (int i = 1; i < N - 1; i++) {
            T_newcase2[i] = (Kcase2[i - 1] * Tempcase2[i - 1] + Kcase2[i + 1] * Tempcase2[i + 1]) / (Kcase2[i + 1] + Kcase2[i - 1]);
        }
        // Check for convergence and update temperatures
        for (int i = 0; i < N; i++) {
            double dt = fabs(Tempcase2[i] - T_newcase2[i]);
            if (dt > max_dt) {
                max_dt = dt;
            }
            Tempcase2[i] = T_newcase2[i];
            Kcase2[i] = 0.1 * (1 + B * Tempcase2[i]);
        }
        n++;
    }

    cout << "Number of iterations: " << n << "\n";
    for (int i = 0; i < N; i++) {
        cout << "Node " << i + 1 << ": Analytical temp = " << T_analyticalcase2[i] << " FDM temp = " << Tempcase2[i] - 273.15 << " \n";
    }

    ofstream myfile2("fdmq5case2sol.txt");
    if (myfile2.is_open()) {
        for (int i = 0; i < N; i++) {
            myfile2 << i + 1 << ":" << Tempcase2[i] - 273.15 << endl;
        }
        myfile2.close();
        cout << "Temperature values have been written to fdmq5case2sol.txt" << endl;
    } else {
        cout << "Error opening file for FDM solution!" << endl;
    }

    ofstream outfile2("fdmq5case2analyticalsol.txt");
    if (outfile2.is_open()) {
        for (int i = 0; i < N; i++) {
            outfile2 << i + 1 << ":" << T_analyticalcase2[i] << endl;
        }
        outfile2.close();
        cout << "Temperature values have been written to fdmq5case2analyticalsol.txt" << endl;
    } else {
        cout << "Error opening file for Analytical solution!" << endl;
    }

    return 0;
}
