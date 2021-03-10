/************************************************************************************************/
/*** Topic: Voltage trace of a WB neuron with an inhibitory autapse.                          ***/
/*** Version Release 17.12 rev 11256                                                          ***/
/*** Date: 3/10/2021                                                               Ali-Seif   ***/
/*** Code implemented in Microsoft Visual Studio Enterprise 2019 C++ compiler                 ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;
//##############################################################
//####                                                      ####
//####               Calculate alpha and betas              ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//______________________________________I-cell____________________________________________________
//________________________________________________________________________________________________
double alpha_n_i(double v_i) {
    return   -0.01 * (v_i + 34.0) / (exp(-0.1 * (v_i + 34.0)) - 1.0);
}
double beta_n_i(double v_i) {
    return   0.125 * exp(-(v_i + 44.0) / 80.0);
}
double alpha_m_i(double v_i) {
    return  0.1 * (v_i + 35.0) / (1.0 - exp(-(v_i + 35.0) / 10.0));
}
double beta_m_i(double v_i) {
    return  4.0 * exp(-(v_i + 60.0) / 18.0);
}
double alpha_h_i(double v_i) {
    return   0.07 * exp(-(v_i + 58.0) / 20.0);
}
double beta_h_i(double v_i) {
    return   1.0 / (exp(-0.1 * (v_i + 28.0)) + 1.0);
}
//##############################################################
//####                                                      ####
//####      Calculate infinite activation variables         ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//______________________________________I-cell____________________________________________________
//________________________________________________________________________________________________
double n_inf_i(double v_i) {
    return   alpha_n_i(v_i) / (alpha_n_i(v_i) + beta_n_i(v_i));
}
double h_inf_i(double v_i) {
    return   alpha_h_i(v_i) / (alpha_h_i(v_i) + beta_h_i(v_i));
}
double m_inf_i(double v_i) {
    return   alpha_m_i(v_i) / (alpha_m_i(v_i) + beta_m_i(v_i));
}
double tau_n_i(double v_i) {
    double tau_n = 1.0 / (alpha_n_i(v_i) + beta_n_i(v_i));
    int phi = 5;
    return   (tau_n / phi);
}
double tau_h_i(double v_i) {
    double tau_h = 1.0 / (alpha_h_i(v_i) + beta_h_i(v_i));
    int phi = 5;
    return    (tau_h / phi);
}
//##############################################################
//####                                                      ####
//####             Calculation of currents                  ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//______________________________________I-cell____________________________________________________
//________________________________________________________________________________________________
double INa_i(double v_i, double h_i) {
    int   g_Na_i = 35;
    int   E_Na_i = 55;
    return    g_Na_i * h_i * pow(m_inf_i(v_i), 3) * (v_i - E_Na_i);
}
double IK_i(double v_i, double n_i) {
    int   g_K = 9;
    int   E_K = -90;
    return    g_K * pow(n_i, 4) * (v_i - E_K);
}
double Il_i(double v_i) {
    float g_l = 0.1;
    int   E_l = -65;
    return        g_l * (v_i - E_l);
}
double Isyn_i(double v_i, double g_ii, double s_i,double v_rev_i) {
    int   E_rev_i = v_rev_i;
    return        g_ii * s_i * (E_rev_i - v_i);
}

//##############################################################
//####                                                      ####
//####            Differential Equations                    ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//______________________________________I-cell____________________________________________________
//________________________________________________________________________________________________
double dvdt_i(double t, double v_i, double n_i, double h_i, double g_ii, double s_i,double I_i,double v_rev_i) {
    float I_app_i = I_i;
    float C_m = 1.0;
    return   (1 / C_m) * (I_app_i - (INa_i(v_i, h_i) + IK_i(v_i, n_i) + Il_i(v_i)) + Isyn_i(v_i, g_ii, s_i, v_rev_i));
}
double dndt_i(double t, double n_i, double v_i) {
    return   ((n_inf_i(v_i) - n_i) / tau_n_i(v_i));
}
double dhdt_i(double t, double h_i, double v_i) {
    return   ((h_inf_i(v_i) - h_i) / tau_h_i(v_i));
}
double dqdt_i(double t, double q_i, double v_i, double tau_dq_i) {
    return   0.5 * (1 + tanh(0.1 * v_i)) * (1.0 - q_i) * 10.0 - q_i / tau_dq_i;
}
double dsdt_i(double t, double s_i, double q_i, double v_i, double tau_r_i, double tau_d_i) {
    return   q_i * (1.0 - s_i) / tau_r_i - s_i / tau_d_i;
}
//##############################################################
//####                                                      ####
//####               tau_dq_function                        ####
//####                                                      ####
//##############################################################
//________________________________________________________________________________________________
//_____________________________________tau_peak___________________________________________________
//________________________________________________________________________________________________
double tau_peak_function(double tau_d, double tau_r, double tau_d_q) {

    double dt = 0.01;
    double dt05 = dt / 2;
    double s = 0;
    double t = 0;
    double s_inc = exp(-t / tau_d_q) * (1 - s) / tau_r - s * tau_d;
    double t_old;
    double s_inc_old;
    double s_tmp;
    double s_inc_tmp;
    while (s_inc > 0) {
        t_old = t;
        s_inc_old = s_inc;
        s_tmp = s + dt05 * s_inc;
        s_inc_tmp = exp(-(t + dt05) / tau_d_q) * (1 - s_tmp) / tau_r - s_tmp / tau_d;
        s = s + dt * s_inc_tmp;
        t = t + dt;
        s_inc = exp(-t / tau_d_q) * (1 - s) / tau_r - s / tau_d;
    }
    return   (t_old * (-s_inc) + t * s_inc_old) / (s_inc_old - s_inc);
}
//________________________________________________________________________________________________
//______________________________________tau_dq____________________________________________________
//________________________________________________________________________________________________
double tau_dq_function(double tau_d, double tau_r, double tau_hat) {
    double tau_d_q_left = 1.0;
    while (tau_peak_function(tau_d, tau_r, tau_d_q_left) > tau_hat) {
        tau_d_q_left = tau_d_q_left / 2;
    }

    double tau_d_q_right = tau_r;
    while (tau_peak_function(tau_d, tau_r, tau_d_q_right) < tau_hat) {
        tau_d_q_right = tau_d_q_right * 2;
    }
    double tau_d_q_mid;

    while (tau_d_q_right - tau_d_q_left > pow(10, -12)) {
        tau_d_q_mid = (tau_d_q_left + tau_d_q_right) / 2;
        if (tau_peak_function(tau_d, tau_r, tau_d_q_mid) <= tau_hat) {
            tau_d_q_left = tau_d_q_mid;

        }
        else {
            tau_d_q_right = tau_d_q_mid;

        }
    }
    return (tau_d_q_left + tau_d_q_right) / 2;
}
//_______________________________________________________________________________________\\
//_____________              The principle of the program                   _____________\\
//_____________                                      @                      _____________\\
//_____________           @@       @@       @            @@     @           _____________\\
//_____________           @ @     @ @      @ @       @   @ @    @           _____________\\
//_____________           @  @   @  @     @   @      @   @  @   @           _____________\\
//_____________           @   @@@   @    @@@@@@@     @   @   @  @           _____________\\
//_____________           @    @    @   @       @    @   @    @ @           _____________\\
//_____________           @         @  @         @   @   @     @@           _____________\\
//_______________________________________________________________________________________
int main() {
    //time variables
    double  
        t0 = 0, 
        t_final = 201, 
        dt = 0.01, 
        dt2 = dt / 2;
    //Define network parameters
    double  
        I_i = 1.5,
        g_ii = 0.5,
        v_rev_i = -75.0,
        tau_r_i = 0.5,
        tau_peak_i = 0.5,
        tau_d_i = 9;
    //function of tau dqi
    double  tau_dq_i = tau_dq_function(tau_d_i, tau_r_i, tau_peak_i); //0.1163;
    //initial conditions
    double 
        v_i = -75.0,
        h_i = 0.1,
        n_i = 0.1,
        q_i = 0,
        s_i = 0;

    ofstream temp("temp.txt", ios::out | ios::trunc);
    double v_inc_i = 0.0, v_temp_i = 0.0, h_inc_i = 0.0, h_temp_i = 0.0, n_inc_i = 0.0, n_temp_i = 0.0, q_inc_i = 0.0, q_temp_i = 0.0, s_inc_i = 0.0, s_temp_i = 0.0;
    //Solve the system using the midpoint method
    for (t0 = 0; t0 <= t_final; t0 = t0 + dt) {

        v_inc_i = dvdt_i(t0, v_i, n_i, h_i, g_ii, s_i, I_i, v_rev_i);
        n_inc_i = dndt_i(t0, n_i, v_i);
        h_inc_i = dhdt_i(t0, h_i, v_i);
        q_inc_i = dqdt_i(t0, q_i, v_i, tau_dq_i);
        s_inc_i = dsdt_i(t0, s_i, q_i, v_i, tau_r_i, tau_d_i);

        v_temp_i = v_i + dt2 * v_inc_i;
        h_temp_i = h_i + dt2 * h_inc_i;
        n_temp_i = n_i + dt2 * n_inc_i;
        q_temp_i = q_i + dt2 * q_inc_i;
        s_temp_i = s_i + dt2 * s_inc_i;

        v_inc_i = dvdt_i(t0, v_temp_i, n_temp_i, h_temp_i, g_ii, s_temp_i, I_i, v_rev_i);
        n_inc_i = dndt_i(t0, n_temp_i, v_temp_i);
        h_inc_i = dhdt_i(t0, h_temp_i, v_temp_i);
        q_inc_i = dqdt_i(t0, q_temp_i, v_temp_i, tau_dq_i);
        s_inc_i = dsdt_i(t0, s_temp_i, q_temp_i, v_temp_i, tau_r_i, tau_d_i);

        v_i = v_i + dt * v_inc_i;
        h_i = h_i + dt * h_inc_i;
        n_i = n_i + dt * n_inc_i;
        q_i = q_i + dt * q_inc_i;
        s_i = s_i + dt * s_inc_i;

        temp << t0  << '\t' << v_i << endl;
    }
    temp.close();
    cout << "\nFinish" << endl;
    return 0;
}