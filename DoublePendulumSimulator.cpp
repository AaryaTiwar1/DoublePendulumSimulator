
Copyright (C) [1/3/2024] [Aarya Tiwari]

DOUBLE PENDULUM SIMULATOR IN C++ [Written and signed using JetBrains CLion Education]
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * 
 *
 
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>

using namespace std;

// Define constants
const double g = 9.81;  // Acceleration due to gravity (m/s^2)

// Function to calculate derivatives (angular velocities)
void derivatives(const complex<double>& t, const complex<double>& y,
                 complex<double>& dydt, double l1, double l2, double m1, double m2) {
    complex<double> delta_theta = y[2] - y[0];
    complex<double> sin_delta_theta = sin(delta_theta);
    complex<double> cos_delta_theta = cos(delta_theta);

    complex<double> denominator1 = (m1 + m2) * l1 - m2 * l1 * cos_delta_theta * 
    cos_delta_theta;
    complex<double> denominator2 = (m1 + m2) * l2 - m2 * l2 * cos_delta_theta * 
    cos_delta_theta;

    dydt[0] = y[1];
    dydt[1] = (m2 * l2 * y[3] * y[3] * sin_delta_theta * cos_delta_theta +
               (m1 + m2) * g * sin(y[2]) - m2 * l2 * y[3] * y[3] * sin(y[0])) / 
               denominator1;

    dydt[2] = y[3];
    dydt[3] = ((m1 + m2) * l1 * y[1] * y[1] * sin_delta_theta * cos_delta_theta +
               (m1 + m2) * g * sin(y[0]) + m2 * l1 * y[1] * y[1] * sin(y[2])) / 
               denominator2;
}

// Function to perform Euler integration
void euler(const complex<double>& t, complex<double>& y, double h,
           double l1, double l2, double m1, double m2) {
    complex<double> dydt[4];
    derivatives(t, y, dydt, l1, l2, m1, m2);

    for (int i = 0; i < 4; ++i) {
        y[i] += h * dydt[i];
    }
}

// Function to perform Runge-Kutta integration
void runge_kutta(const complex<double>& t, complex<double>& y, double h,
                 double l1, double l2, double m1, double m2) {
    complex<double> k1[4], k2[4], k3[4], k4[4];

    derivatives(t, y, k1, l1, l2, m1, m2);
    for (int i = 0; i < 4; ++i) k2[i] = y[i] + 0.5 * h * k1[i];

    derivatives(t + 0.5 * h, k2, k2, l1, l2, m1, m2);
    for (int i = 0; i < 4; ++i) k3[i] = y[i] + 0.5 * h * k2[i];

    derivatives(t + 0.5 * h, k3, k3, l1, l2, m1, m2);
    for (int i = 0; i < 4; ++i) k4[i] = y[i] + h * k3[i];

    for (int i = 0; i < 4; ++i) {
        y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

int main() {
    // Get initial conditions from user
    complex<double> y[4];
    cout << "Enter initial conditions (theta1 omega1 theta2 omega2): ";
    cin >> y[0] >> y[1] >> y[2] >> y[3];

    // Pendulum parameters
    double l1 = 1.0, l2 = 1.0, m1 = 1.0, m2 = 1.0;

    // Time step and total simulation time
    double h = 0.01, total_time = 10.0;

    // Output to a file for visualization
    ofstream output_file("double_pendulum_simulation_complex.dat");
    output_file << "# time theta1 omega1 theta2 omega2\n";

    complex<double> t(0.0, 0.0);  // Complex time for future expansion

    for (double current_time = 0.0; current_time <= total_time; current_time += h) {
        // Output current state to file
        output_file << current_time << " " << y[0].real() << " " << y[1].real() << " "
                     << y[2].real() << " " << y[3].real() << "\n";

        // Perform Euler integration
        euler(t, y, h, l1, l2, m1, m2);

        // Update complex time
        t += complex<double>(h, 0.0);
    }

    output_file.close();

    cout << "Simulation complete. Data saved to double_pendulum_simulation_complex.dat\n";

    return 0;
}


\end{verbatim}

\subsection*{Python Simulation}
Python is a more widestream language, however it is less scientifically accurate. However if only simple Euler methods are to be used to simulate a double pendulum in Python, you can use the following script. The script utilizes the \texttt{matplotlib} library for plotting and the \texttt{scipy.integrate} module for solving the differential equations governing the double pendulum's motion.

\begin{verbatim}
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def double_pendulum(t, y, l1, l2, m1, m2, g):
    # ... (equations of motion)

# Example usage
initial_conditions = [np.pi/4, 0, np.pi/2, 0]
length1 = 1.0
length2 = 1.0
mass1 = 1.0
mass2 = 1.0
gravity = 9.8

# ... (function to simulate and plot)

# Example usage
simulate_and_plot(initial_conditions, length1, length2, mass1, mass2, gravity, 

'double_pendulum_simulation.png')
