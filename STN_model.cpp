#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define CHECK 4
#define EPS6  1.e-6

//model parameters
double g_L = 2.25;
double g_K = 45.0;
double g_Na = 0.0;
double g_T = 0.5;
double g_Ca = 0.5;
double g_AHP = 9.0;
double v_L = -60.0;
double v_K = -80.0;
double v_Na = 55.0;
double v_Ca = 140.0;
double tau_1_h = 500.0;
double tau_1_n = 100.0;
double tau_1_r = 17.5;
double tau_0_h = 1.0;
double tau_0_n = 1.0;
double tau_0_r = 40.0;
double phi_h = 0.75;
double phi_n = 0.75;
double phi_r = 0.2;
double k_1 = 15.0;
double k_Ca = 22.5;
double epsilon = 3.75e-5;
double theta_m = -30.0;
double theta_h = -39.0;
double theta_n = -32.0;
double theta_r = -67.0;
double theta_a = -63.0;
double theta_b = 0.4;
double theta_s = -39.0;
double theta_tau_h = -57.0;
double theta_tau_n = -80.0;
double theta_tau_r = 68.0;
double theta_H_g = -39.0;
double theta_g = 30.0;
double alpha = 5.0;
double v_GS = -85.0;
double sigma_m = 15.0;
double sigma_h = -3.1;
double sigma_n = 8.0;
double sigma_r = -2.0;
double sigma_a = 7.8;
double sigma_b = -0.1;
double sigma_s = 8.0;
double sigma_tau_h = -3.0;
double sigma_tau_n = -26.0;
double sigma_tau_r = -2.2;
double sigma_H_g = 8.0;
double beta = 1.0;

void STN_neuron();
double f_v(double v, long check);
double Current_L(double v);
double differential_n(double v, long check);
double Rungekkuta_n(double tau, double n_inf, long check);
double f_n(double n_i, double tau, double n_inf);
double Function_tau_n(double v);
double Infinity_n(double v);
double differential_h(double v, long check);
double Rungekkuta_h(double tau, double h_inf, long check);
double f_h(double h_i, double tau, double h_inf);
double Function_tau_h(double v);
double Infinity_h(double v);
double differential_r(double v, long check);
double Rungekkuta_r(double tau, double r_inf, long check);
double f_r(double r_i, double tau, double r_inf);
double Function_tau_r(double v);
double Infinity_r(double v);
double Current_K(double v, long check);
double Current_Na(double v, long check);
double Infinity_m(double v);
double Current_T(double v, long check);
double Infinity_a(double v);
double Infinity_b(double r_b);
double Current_Ca(double v);
double Infinity_s(double v);
double Current_AHP(double v, double I_Ca, double I_T, long check);
void Calucium_concentration(double I_Ca, double I_T, long check);
double f_ca(double ca, double I_Ca, double I_T);
double Current_GS(double v, long check);
double Sumtation_GPe(double v_gj, long check);
double f_GPe(double s_j, double H_infinity);
double Infinity_H(double v);


double H = 0.05e0;
double v_STN = -71.4;

//voltage
double n[CHECK];
double h[CHECK];
double r[CHECK];
double Ca[CHECK];
double GPe_s_j[CHECK];

//current
double I_L;
double I_K;
double I_Na;
double I_T;
double I_Ca;
double I_AHP;
double I_GS;

FILE *fp;

int main(void){
  	long i;
  	double t=0;
	long ns;
	
  	for(i = 0; i < 4; i++){
      	Ca[i] = 1.0e-7; h[i] = 0;
      	n[i] = 0; r[i] = 0;
      	GPe_s_j[i] = 1.0;
    }
  
  	if((fp = fopen("data.dat", "w"))==NULL){
      	printf("error in fopen\n");
      	exit(1);
    }
  
    printf("-------- t = %.1lf -------\n", t);
    printf("v_STN = %.7lf\n\n", v_STN);
	fprintf(fp, "t, v_STN\n");
    fprintf(fp, "%.7lf, %.7lf\n", t, v_STN);
  
  	for(;;){
      	t += H;

      	STN_neuron();

		fprintf(fp, "%.7lf, %.7lf\n", t, v_STN);

      	if(t > 100)break;
    }

  	return 0;
}

//voltage
void STN_neuron(){
  	double k1, k2, k3, k4;

  	k1 = H * f_v(v_STN, 0);
  	k2 = H * f_v(v_STN + k1/2.0, 1);
  	k3 = H * f_v(v_STN + k2/2.0, 2);
  	k4 = H * f_v(v_STN + k3, 3);

  	v_STN = v_STN + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
}

double f_v(double v, long check){
  	I_L = Current_L(v);
  	I_K = Current_K(v, check);
  	I_Na = Current_Na(v, check);
  	I_T = Current_T(v, check);
  	I_Ca = Current_Ca(v);
 	I_AHP = Current_AHP(v, I_Ca, I_T, check);
 	I_GS = Current_GS(v, check);

	return (-I_L -I_K -I_Na -I_T -I_Ca -I_AHP -I_GS);
}

double Current_L(double v){
	return g_L * (v - v_L);
}

//calculation of factor n
double differential_n(double v, long check){
	double Rungekkuta, tau, n_inf;

	tau = Function_tau_n(v);
 
	n_inf = Infinity_n(v);
 
 	Rungekkuta = Rungekkuta_n(tau, n_inf, check);

 	return Rungekkuta;
}

double Rungekkuta_n(double tau, double n_inf, long check){
	double k1, k2, k3, k4;

	k1 = H * f_n(n[check], tau, n_inf);
	k2 = H * f_n(n[check] + k1/2.0, tau, n_inf);
	k3 = H * f_n(n[check] + k2/2.0, tau, n_inf);
	k4 = H * f_n(n[check] + k3, tau, n_inf);

	return n[check] + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
}

double f_n(double n_i, double tau, double n_inf){
	return phi_n * (n_inf - n_i) / tau;
}

//tau_n(v)
double Function_tau_n(double v){
	double exponential;

	exponential = exp(-(v - theta_tau_n) / sigma_tau_n);
 
	return tau_0_n + (tau_1_n / (1 + exponential));
}

//X = n
double Infinity_n(double v){
	double exponential;

	exponential = exp(-(v - theta_n) / sigma_n);
 
	return 1 / (1 + exponential);
}

//calculation of factor h
double differential_h(double v, long check){
	double Rungekkuta, tau, h_inf;

	tau = Function_tau_h(v);
 
	h_inf = Infinity_h(v);

 	Rungekkuta = Rungekkuta_h(tau, h_inf, check);

	return Rungekkuta;
}

double Rungekkuta_h(double tau, double h_inf, long check){
	double k1, k2, k3, k4;

	k1 = H * f_h(h[check], tau, h_inf);
	k2 = H * f_h(h[check] + k1/2.0, tau, h_inf);
	k3 = H * f_h(h[check] + k2/2.0, tau, h_inf);
	k4 = H * f_h(h[check] + k3, tau, h_inf);

	return h[check] + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
}

double f_h(double h_i, double tau, double h_inf){
	return phi_h * (h_inf - h_i) / tau;
}

//tau_h(v)
double Function_tau_h(double v){
	double exponential;

	exponential = exp(-(v - theta_tau_h) / sigma_tau_h);
 
	return tau_0_h + (tau_1_h / (1 + exponential));
}

//Infinity_h
double Infinity_h(double v){
	double exponential;

	exponential = exp(-(v - theta_h) / sigma_h);
 
	return 1 / (1 + exponential);
}

//calculation of factor r
double differential_r(double v, long check){
	double Rungekkuta, tau, r_inf;

	tau = Function_tau_r(v_STN);
 
	r_inf = Infinity_r(v_STN);

 	Rungekkuta = Rungekkuta_r(tau, r_inf, check);

	return Rungekkuta;
}

double Rungekkuta_r(double tau, double r_inf, long check){
	double k1, k2, k3, k4;

	k1 = H * f_r(r[check], tau, r_inf);
	k2 = H * f_r(r[check] + k1/2.0, tau, r_inf);
	k3 = H * f_r(r[check] + k2/2.0, tau, r_inf);
	k4 = H * f_r(r[check] + k3, tau, r_inf);

	return r[check] + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
}

double f_r(double r_i, double tau, double r_inf){
	return phi_r * (r_inf - r_i) / tau;
}

//tau_r(v)
double Function_tau_r(double v){
	double exponential;

	exponential = exp(-(v - theta_tau_r) / sigma_tau_r);
 
	return tau_0_r + (tau_1_r / (1 + exponential));
}

//Infinity_r
double Infinity_r(double v){
	double exponential;

	exponential = exp(-(v - theta_r) / sigma_r);
 
	return 1 / (1 + exponential);
}

//I_K
double Current_K(double v, long check){
	double n4;

 	n[check] = differential_n(v_STN - v_K, check);

	n4 = n[check] * n[check] * n[check] * n[check];

	return g_K * n4;
}

//I_Na
double Current_Na(double v, long check){
	double m_inf, m3, h_Na;

	m_inf = Infinity_m(v);

	m3 = m_inf * m_inf * m_inf;

 	h[check] = differential_h(v_STN - v_Na, check);
 
	return g_Na * m3 * h[check];
}

//m
double Infinity_m(double v){
	double exponential;

	exponential = exp(-(v - theta_m) / sigma_m);

	return 1 / (1 + exponential);
}

//I_T
double Current_T(double v, long check){
	double a, b, a3, b2;

 	r[check] = differential_r(v_STN, check);

	a = Infinity_a(v);

	b = Infinity_b(r[check]);

	a3 = a * a * a;

	b2 = b * b;

	return g_T * a3 * b2 * (v - v_Ca);
}

//X = a
double Infinity_a(double v){
	double exponential;

	exponential = exp(-(v - theta_a) / sigma_a);

	return 1 / (1 + exponential);
}

//b_Åá(r)
double Infinity_b(double r_b){
	double exponential_1, exponential_2;
	double x, y;

	exponential_1 = exp((r_b - theta_b) / sigma_b);
	exponential_2 = exp(- theta_b / sigma_b);

	x = 1 / (1 + exponential_1);
	y = 1 / (1 + exponential_2);
 
	return x - y;
}

//I_Ca
double Current_Ca(double v){
	double s, s2;

	s = Infinity_s(v);

	s2 = s * s;

	return g_Ca * s2 * (v - v_Ca);
}

//X = s
double Infinity_s(double v){
	double exponential;

	exponential = exp(-(v - theta_s) / sigma_s);

	return 1 / (1 + exponential);
}

//I_AHP
double Current_AHP(double v, double I_Ca, double I_T, long check){
	double calcium;

 	Calucium_concentration(I_Ca, I_T, check);

	calcium = Ca[check] / (Ca[check] + k_1);

	return g_AHP * (v - v_K) * calcium;
}

void Calucium_concentration(double I_Ca, double I_T, long check){
	double k1, k2, k3, k4;

	k1 = H * f_ca(Ca[check], I_Ca, I_T);
	k2 = H * f_ca(Ca[check] + k1/2.0, I_Ca, I_T);
	k3 = H * f_ca(Ca[check] + k2/2.0, I_Ca, I_T);
	k4 = H * f_ca(Ca[check] + k3, I_Ca, I_T);

	Ca[check] = Ca[check] + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
}

double f_ca(double ca, double I_Ca, double I_T){
	return epsilon * (-I_Ca -I_T -k_Ca * ca);
}

//I_GS
double Current_GS(double v, long check){
	double GPe_synaps;

	/*GPe_s_j[check] = Sumtation_GPe(30, check);*/
 
	return 1.0 * (v - v_GS) * 1.0;
}

double Sumtation_GPe(double v_gj, long check){
	double H_inf;
	double k1, k2, k3, k4;
 
	H_inf = Infinity_H(v_gj - theta_g);

 	k1 = H * f_GPe(GPe_s_j[check], H_inf);
 	k2 = H * f_GPe(GPe_s_j[check] + k1/2.0, H_inf);
 	k3 = H * f_GPe(GPe_s_j[check] + k2/2.0, H_inf);
 	k4 = H * f_GPe(GPe_s_j[check] + k3, H_inf);

 	return GPe_s_j[check] + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
}

double f_GPe(double s_j, double H_infinity){
  	return alpha * H_infinity * (1 - s_j) - beta * s_j;
}

//H_Åá(v)
double Infinity_H(double v){
	double exponential;

	exponential = exp(-(v - theta_H_g) / sigma_H_g);

	return 1 / (1 + exponential);
}