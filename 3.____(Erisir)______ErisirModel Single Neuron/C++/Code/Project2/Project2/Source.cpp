/************************************************************************************************/
/*** Topic: modified Erisir model with Runge-Kutta 4th Order Method for one neuron            ***/
/*** Version Release 17.12 rev 11256                                                Ali-Seif  ***/
/*** Date: 3/1/2021                                                                           ***/
/*** Code implemented in CodeBlocks C++ compiler (v. 17.12),                                  ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;
//_________________________________Calculate alpha and betas_________________________________//


double alpha_n(double v) {
    return   (95.0 - v) / ((exp((95.0 - v) / 11.8)) - 1);
}

double beta_n(double v) {
    return   0.025 * exp(-(v) / 22.222);
}

double alpha_m(double v) {
    return  40.0 * (75.5 - v) / ((exp((75.5 - v) / 13.5)) - 1.0);
}

double beta_m(double v) {
    return  1.2262 * exp(-(v) / 42.248);
}

double alpha_h(double v) {
    return   0.0035 * exp(-(v) / 24.186);
}

double beta_h(double v) {
    return   -0.017 * (v + 51.25) / ((exp(-(v + 51.25) / 5.2)) - 1.0);
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

double INa1(double v, double m, double h) {
    int   g_Na = 112;
    int   E_Na = 60;
    return    g_Na * h * pow(m, 3) * (v - E_Na);
}

double IK1(double v, double n) {
    int   g_K = 224;
    int   E_K = -90;
    return    g_K * pow(n, 4) * (v - E_K);
}

double Il1(double v) {
    float g_l = 0.5;
    int   E_l = -70;
    return        g_l * (v - E_l);
}
//___________________________________Differential Equations__________________________________//

double dvdt(double t, double v, double m, double n, double h) {
    float Iapp = 7;
    float C_m = 1.0;
    return   (1 / C_m) * (Iapp - (INa1(v, m, h) + IK1(v, n) + Il1(v)));
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

//*******************************************************************************************//
//                                                                                           //
//________________________________The principle of the program_______________________________//
//                                                                                           //
//*******************************************************************************************//
int main() {

    double  t0 = 0, t_final = 101, v = -70.0, dt = 0.01, dt2 = dt / 2, n, h, m;
    m = m_inf(v);
    n = n_inf(v);
    h = h_inf(v);
    ofstream temp("temp.txt", ios::out | ios::trunc);
    double v_inc = 0.0, v_temp = 0.0, m_inc = 0.0, m_temp = 0.0, h_inc = 0.0, h_temp = 0.0, n_inc = 0.0, n_temp = 0.0;

    for (t0 = 0; t0 <= t_final; t0 = t0 + dt) {

        v_inc = dvdt(t0, v, m, n, h);
        m_inc = dmdt(t0, m, v);
        h_inc = dhdt(t0, h, v);
        n_inc = dndt(t0, n, v);

        v_temp = v + dt2 * v_inc;
        m_temp = m_inf(v);
        h_temp = h + dt2 * h_inc;
        n_temp = n + dt2 * n_inc;


        v_inc = dvdt(t0, v_temp, m_temp, n_temp, h_temp);
        m_inc = dmdt(t0, m_temp, v_temp);
        h_inc = dhdt(t0, h_temp, v_temp);
        n_inc = dndt(t0, n_temp, v_temp);

        v = v + dt * v_inc;
        m = m_inf(v);
        h = h + dt * h_inc;
        n = n + dt * n_inc;
        temp << t0 << '\t' << v << endl;


    }
    temp.close();
    cout << "\nFinish" << endl;
    return 0;
}