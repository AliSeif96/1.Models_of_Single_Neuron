/************************************************************************************************/
/*** Topic: Wang-Buzsaki model with Runge-Kutta 4th Order Method for one neuron    Ali-Seif   ***/
/*** Version Release 17.12 rev 11256                                                          ***/
/*** Date: 2/27/2021                                                                          ***/
/*** Code implemented in Microsoft Visual Studio Enterprise 2019 C++ compiler                 ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;
//_________________________________Calculate alpha and betas_________________________________//

double alpha_n(double v) {
    return   0.01 * (v + 60.0) / (1.0-exp(-0.1 * (v + 60.0)));
}

double beta_n(double v) {
    return   0.125 * exp(-(v + 70.0) / 80.0);
}

double alpha_m(double v) {
    double ab = abs(v + 45.0);
    double alpha = 0.0;

    if (ab > 0.00000001) {
        alpha = 0.1 * (v + 45.0) / (1.0 - exp(-0.1 * (v + 45.0)));
    }
    else {
        alpha = 1;
    }


    return   alpha;
}

double beta_m(double v) {
    return   4.0 * exp(-(v + 70.0) / 18.0);
}

double alpha_h(double v) {
    return   0.07 * exp(-(v + 70.0) / 20.0);
}

double beta_h(double v) {
    return   1.0 / (exp(-(v + 40.0)/10) + 1.0);
}
//__________________________Calculate infinite activation variables__________________________//

double n_inf(double v) {
    return   alpha_n(v) / (alpha_n(v) + beta_n(v));
}

double h_inf(double v) {
    return   alpha_h(v) / (alpha_h(v) + beta_h(v));
}

double m_inf(double v) {
    return   alpha_m(v) / (alpha_m(v) + beta_m(v));
}
//__________________________________Calculation of currents__________________________________//

double INa1(double v, double h) {
    int   g_Na = 120;
    int   E_Na = 45;
    return    g_Na * h * pow(m_inf(v), 3) * (v - E_Na);
}

double IK1(double v, double n) {
    int   g_K = 36;
    int   E_K = -82;
    return    g_K * pow(n, 4) * (v - E_K);
}

double Il1(double v) {
    float g_l = 0.3;
    int   E_l = -59;
    return        g_l * (v - E_l);
}
//___________________________________Differential Equations__________________________________//

double dvdt(double t, double v, double n, double h) {
    float Iapp = 7;
    float C_m = 1.0;
    return   (1 / C_m) * (Iapp - (INa1(v, h) + IK1(v, n) + Il1(v)));
}

double dmdt(double t, double m, double v) {
    int   phi = 1;
    return   phi * ((alpha_m(v) * (1 - m)) - beta_m(v) * m);
}

double dndt(double t, double n, double v) {
    int   phi = 1;
    return   phi * ((alpha_n(v) * (1 - n)) - beta_n(v) * n);
}

double dhdt(double t, double h, double v) {
    int   phi = 1;
    return   phi * ((alpha_h(v) * (1 - h)) - beta_h(v) * h);
}
//__________________________________Runge-Kutta calculations_________________________________//

double rk4thOrder_v(double t0, double v, double dt, double n, double h) {
    double  k1, k2, k3, k4;
    k1 = dt * dvdt(t0, v, n, h);
    k2 = dt * dvdt((t0 + dt / 2), (v + k1 / 2), n, h);
    k3 = dt * dvdt((t0 + dt / 2), (v + k2 / 2), n, h);
    k4 = dt * dvdt((t0 + dt), (v + k3), n, h);
    v = v + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));
    return   v;
}

double rk4thOrder_m(double t0, double v, double dt, double m) {

    double  k1, k2, k3, k4;
    k1 = dt * dmdt(t0, m, v);
    k2 = dt * dmdt((t0 + dt / 2), (m + k1 / 2), v);
    k3 = dt * dmdt((t0 + dt / 2), (m + k2 / 2), v);
    k4 = dt * dmdt((t0 + dt), (m + k3), v);
    m = m + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));
    return   m;
}

double rk4thOrder_n(double t0, double v, double dt, double n) {

    double  k1, k2, k3, k4;
    k1 = dt * dndt(t0, n, v);
    k2 = dt * dndt((t0 + dt / 2), (n + k1 / 2), v);
    k3 = dt * dndt((t0 + dt / 2), (n + k2 / 2), v);
    k4 = dt * dndt((t0 + dt), (n + k3), v);
    n = n + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));
    return   n;
}

double rk4thOrder_h(double t0, double v, double dt, double h) {

    double  k1, k2, k3, k4;
    k1 = dt * dhdt(t0, h, v);
    k2 = dt * dhdt((t0 + dt / 2), (h + k1 / 2), v);
    k3 = dt * dhdt((t0 + dt / 2), (h + k2 / 2), v);
    k4 = dt * dhdt((t0 + dt), (h + k3), v);
    h = h + double((1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4));
    return   h;
}
//*******************************************************************************************//
//                                                                                           //
//________________________________The principle of the program_______________________________//
//                                                                                           //
//*******************************************************************************************//
int main() {

    double  t0 = 0, t_final = 200, v = -70.0, dt = 0.01, n, h, m;
    m = m_inf(v);
    n = 0.6;
    h = 0.7;
    ofstream temp("temp.txt", ios::out | ios::trunc);

    for (t0 = 0; t0 <= t_final; t0 = t0 + dt) {

        v = rk4thOrder_v(t0, v, dt, n, h);
        m = rk4thOrder_m(t0, v, dt, m);
        n = rk4thOrder_n(t0, v, dt, n);
        h = rk4thOrder_h(t0, v, dt, h);
        temp << t0 << '\t' << v << endl;
    }
    temp.close();
    cout << "\nFinish" << endl;
    return 0;
}