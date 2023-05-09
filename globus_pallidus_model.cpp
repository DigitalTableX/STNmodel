#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <time.h>

//mode
#define CAL_MODE 0		//0:V_G¨S

#define EPS6	1.e-6
#define EPS9	1.e-9

#define	STATE_MAX	8

struct STATE_stn;

class STATE_stn{
public:
	double v;			//membrane potential
	double h;			//I_Na
	double n;			//I_K
	double r;			//I_T
	double ca;			//Ca2+
	double hn;			//NMDA
	double hd;			//I_DIC
	double ld;			//Ca2+
};

//variable
double Gl, Vl;										//I_leak

double Gna, Vna;									//I_Na
double Phi_h;
double Tau_0h, Tau_1h, Theta_tau_h, Sigma_tau_h;
double Theta_h, Sigma_h;
double Theta_m, Sigma_m;

double Gk, Vk;										//I_K
double Phi_n;
double Tau_0n, Tau_1n;
double Theta_tau_n, Sigma_tau_n;
double Theta_n, Sigma_n;

double Gt, Vca;										//I_T
double Phi_r;
double Tau_0r, Tau_1r, Theta_tau_r, Sigma_tau_r;
double Theta_r, Sigma_r;
double Theta_a, Sigma_a;
double Theta_b, Sigma_b;

double Gca;											//I_Ca
double Theta_s, Sigma_s;

double Gahp;										//I_AHP
double Epsilon, K1, Kca;

double Gn, Vn;										//I_NMDA
double Mn, Eta, Mg, Gam;
double Theta_nh, Sigma_nh, Tau_hn;
double Tau_ca_n;
double Kn1;

double Gd, Vd;										//I_DIC
double Kd;
double Theta_d, Sigma_d;
double Tau_0d, Tau_1d, Theta_tau_d, Sigma_tau_d;

double Ggs, Vgs;									//I_G¨S
double Alpha, Theta_g, Beta;
double Theta_h_g, Sigma_h_g;

//globus pallidus variable
double Gl_g, Vl_g;									//I_leak

double Gk_g, Vk_g;									//I_K
double Phi_n_g;
double Tau_0n_g, Tau_1n_g, Theta_t_n_g, Sigma_t_n_g;
double Theta_n_g, Sigma_n_g;

double Gna_g, Vna_g;								//I_Na
double Phi_h_g;
double Tau_0h_g, Tau_1h_g, Theta_t_h_g, Sigma_t_h_g;
double Theta_h_g_g, Sigma_h_g_g;
double Theta_m_g, Sigma_g;

double Gca_g, Vca_g;								//I_ca
double Theta_s_g, Sigma_s_g;

double Gahp_g, K1_g;								//I_AHP
double Epsilon_g, K_ca_g;

double Gt_g;										//I_T
double Phi_r_g;
double Theta_r_g, Sigma_r_g;
double Theta_a_g, Sigma_a_g;

double Ggg, Vgg;									//I_G¨G
double Alpha_gg, Beta_gg;
double Theta_hg_gg, Sigma_hg_gg;

double Gsg, Vsg;									//I_S¨G
double Alpha_sg, Beta_sg;
double Theta_hg_sg, Sigma_hg_sg;

double V_init;
double H_init;
double N_init;
double R_init;
double Ca_init;
double Hn_init;
double Hd_init;
double Ld_init;

double T_init;										//start time
double T_step;
double T_fin;										//final time
long Step_gamen;

double I_inj_stn;									//current to STN neurons

FILE *Fp_state;
FILE *Fp_current;

//function
void STN_model();
void init_QSM(STATE_stn &xa);
void str_to_arr(STATE_stn xa, double x_all[]);
void arr_to_str(double x_all[], STATE_stn &xa);
void ruku(double t, double x_all[], double h);
void dx_dt(double t, double state[], double K[]);
void calculate(STATE_stn x, STATE_stn &dxdt);
double hinf(double v);
double tauh(double v);
double ninf(double v);
double taun(double v);
double rinf(double v);
double taur(double v);
double hninf(double ld);
double hdinf(double v);
double taud(double v);
double Gnmda(double v, double mn, double hn);
double I_Ca_NMDA(double v, double mn, double hn);
double I_leak(double v);
double minf(double v);
double  I_Na(double v, double h);
double I_K(double v, double n);
double ainf(double v);
double binf(double r);
double I_T(double v, double r);
double sinf(double v);
double I_Ca(double v);
double I_AHP(double v, double ca);
double I_NMDA(double v, double mn, double hn);
double m_Dinf(double ld);
double I_DIC(double v, double ld, double hd);
double I_GPe(void);
double all_current(STATE_stn x);
void output(long ns, double t, STATE_stn x);
void output_state(long ns, double t, STATE_stn x);
void output_current(long ns, double t, STATE_stn x);

int main(void){
	Gl = 2.25, Vl = -60.;							//I_leak
	
	Gna = 37.5, Vna = 55., Phi_h = 0.75;			//I_Na
	Tau_0h = 1.0, Tau_1h = 500.;
	Theta_tau_h = -85., Sigma_tau_h = 3.;
	Theta_h = -39., Sigma_h = 3.1;
	Theta_m = -30., Sigma_m = -14.5;
	
	Gk = 45., Vk = -80., Phi_n = 0.75;				//I_K
	Tau_0n = 1., Tau_1n = 100.;
	Theta_tau_n = -80., Sigma_tau_n = 26.;
	Theta_n = -32., Sigma_n = -8.;
	
	Gt = 0.1, Vca = 120., Phi_r = 0.2;				//I_T
	Tau_0r = 40., Tau_1r = 17.5;
	Theta_tau_r = 68., Sigma_tau_r = 2.2;
	Theta_r = -67., Sigma_r = 2.;
	Theta_a = -63., Sigma_a = -7.8;
	Theta_b = 0.4, Sigma_b = -0.1;
	
	Gca = 0.5;										//I_Ca
	Theta_s = -39., Sigma_s = -8.;
	
	Gahp = 20.;										//I_AHP
	Epsilon = 3.75e-5, K1 = 7.5, Kca = 22.5;
	
	Gn = 20., Vn = -20.;							//I_NMDA
	Mn = 1., Eta = 0.33, Mg = 1.3, Gam = 0.05;
	Theta_nh = 0.7, Sigma_nh = 0.05, Tau_hn = 3000.;
	Tau_ca_n = 80., Kn1 = 0.005;
	
	Gd = 20., Vd = -18.;							//I_DIC
	Kd = 1.15;
	Theta_d = -95., Sigma_d = 14.;
	Tau_0d = 300., Tau_1d = 350.;
	Theta_tau_d = -60., Sigma_tau_d = 3.;
	
	Ggs = 1.0, Vgs = -90.;							//I_G¨S
	Alpha = 5., Theta_g = 30., Beta = 1.;
	Theta_h_g = -39., Sigma_h_g = 8.;
	I_inj_stn = -90.;								//current
	
	T_init = -1.e3;
	T_step = 0.05e0;
	T_fin = 10.e3;
	Step_gamen = 200;
	
	V_init = -71.4;
	H_init = 0.54;
	N_init = 0.0;
	R_init = 0.664;
	Ca_init = 0.446;
	Hn_init = 0.388;
	Hd_init = 0.115;
	Ld_init = 0.613;
	
	if(CAL_MODE == 0)STN_model();
}

//STN model
void STN_model(){
	double t, x_all[STATE_MAX];
	STATE_stn xa;
	long ns;
	
	init_QSM(xa);
	for(ns = 0; ;ns++){
		t = T_init + ns * T_step;
		if(ns % Step_gamen == 0)printf("t = %.7lf\n", t/1000.e0);
		output(ns, t, xa);							//file output
		str_to_arr(xa, x_all);
		ruku(t, x_all, T_step);						//Runge-Kutta
		arr_to_str(x_all, xa);
		if(t > T_fin - EPS6)break;
	}
}

void init_QSM(STATE_stn &xa){
	xa.v = V_init;
	xa.h = H_init;
	xa.n = N_init;
	xa.r = R_init;
	xa.ca = Ca_init;
	xa.hn = Hn_init;
	xa.hd = Hd_init;
	xa.ld = Ld_init;
}

void str_to_arr(STATE_stn xa, double x_all[]){
	x_all[0] = xa.v;
	x_all[1] = xa.h;
	x_all[2] = xa.n;
	x_all[3] = xa.r;
	x_all[4] = xa.ca;
	x_all[5] = xa.hn;
	x_all[6] = xa.hd;
	x_all[7] = xa.ld;
}

void arr_to_str(double x_all[], STATE_stn &xa){
	xa.v = x_all[0];
	xa.h = x_all[1];
	xa.n = x_all[2];
	xa.r = x_all[3];
	xa.ca = x_all[4];
	xa.hn = x_all[5];
	xa.hd = x_all[6];
	xa.ld = x_all[7];
}

void ruku(double t, double x_all[], double h){
	long i;
	double K1[STATE_MAX], K2[STATE_MAX], K3[STATE_MAX], K4[STATE_MAX];
	double state[STATE_MAX];
	
	//1
	dx_dt(t, x_all, K1);
	
	//2
	for(i = 0; i < STATE_MAX; i++)state[i] = x_all[i] + h * K1[i] / 2.e0;
	dx_dt(t + h/2.e0, state, K2);
	
	//3
	for(i = 0; i < STATE_MAX; i++)state[i] = x_all[i] + h * K2[i] / 2.e0;
	dx_dt(t + h/2.e0, state, K3);
	
	//4
	for(i = 0; i < STATE_MAX; i++)state[i] = x_all[i] + h * K3[i];
	dx_dt(t + h, state, K4);
	
	//update
	for(i = 0; i < STATE_MAX; i++){
		x_all[i] = x_all[i] + h * (K1[i] + K2[i]*2.e0 + K3[i]*2.e0 + K4[i]) / 6.e0;
	}
}

void dx_dt(double t, double state[], double K[]){
	STATE_stn x, dxdt;
	arr_to_str(state, x);
	calculate(x, dxdt);
	str_to_arr(dxdt, K);
}

//calculate
void calculate(STATE_stn x, STATE_stn &dxdt){
	dxdt.v = -all_current(x) + I_GPe();
	dxdt.h = Phi_h * (hinf(x.v) - x.h)/tauh(x.v);
	dxdt.n = Phi_n * (ninf(x.v) - x.n)/taun(x.v);
	dxdt.r = Phi_r * (rinf(x.v) - x.r)/taur(x.v);
	dxdt.ca = Epsilon * (-I_Ca(x.v) - I_T(x.v, x.r) - Kca * x.ca);
	dxdt.hn = (-x.hn + hninf(x.ld))/Tau_hn;
	dxdt.hd = (hdinf(x.v) - x.hd) / taud(x.v);
	dxdt.ld = (-I_Ca_NMDA(x.v, Mn, x.hn) - x.ld) / Tau_ca_n;
}

//total current sum
double all_current(STATE_stn x){
	return I_leak(x.v)+I_Na(x.v, x.h)+I_K(x.v, x.n)+I_T(x.v, x.r)+I_Ca(x.v)+I_AHP(x.v, x.ca)+I_NMDA(x.v, Mn, x.hn)+I_DIC(x.v, x.ld, x.hd);
}

//hyperpolarizing current
double I_GPe(void){
	if(CAL_MODE == 0) return I_inj_stn;
}
//h
double hinf(double v){
	return 1. / (1. + exp((v - Theta_h) / Sigma_h));
}

double tauh(double v){
	return Tau_0h + Tau_1h / (1. + exp((v - Theta_tau_h) / Sigma_tau_h));
}

//n
double ninf(double v){
	return 1. / (1. + exp((v - Theta_n) / Sigma_n));
}

double taun(double v){
	return Tau_0n + Tau_1n / (1. + exp((v - Theta_tau_n) / Sigma_tau_n));
}

//r
double rinf(double v){
	return 1. / (1. + exp((v - Theta_r) / Sigma_r));
}

double taur(double v){
	return Tau_0r + Tau_1r / (1. + exp((v - Theta_tau_r) / Sigma_tau_r));
}

//hn
double hninf(double ld){
	return 1. / (1. + exp((ld - Theta_nh) / Sigma_nh));
}

double hdinf(double v){
	return 1. / (1. + exp((v - Theta_d) / Sigma_d));
}

double taud(double v){
	return Tau_0d + Tau_1d / (1. + exp((v-Theta_tau_d) / Sigma_tau_d));
}

//Gnmda
double Gnmda(double v, double mn, double hn){
	return Gn * mn * hn / (1. + Eta * Mg * exp(-Gam * v));
}	

//I_Ca_NMDA
double I_Ca_NMDA(double v, double mn, double hn){
	return Kn1 * Gnmda(v, mn, hn) * (v - Vca);
}

//leak current
double I_leak(double v){
	return Gl * (v - Vl);
}

//I_Na
double  I_Na(double v, double h){
	return Gna * pow(minf(v), 3.0) * h * (v - Vna);
}

//minf
double minf(double v){
	return 1. / (1. + exp((v - Theta_m) / Sigma_m));
}

//I_K
double I_K(double v, double n){
	return Gk * pow(n, 4.0) * (v - Vk);
}
 
//I_T
double I_T(double v, double r){
	return Gt * pow(ainf(v), 3.0) * pow(binf(r), 2.0) * (v - Vca);
}

//ainf
double ainf(double v){
	return 1. / (1. + exp((v - Theta_a) / Sigma_a));
}

//binf
double binf(double r){
	return 1. / (1. + exp((r - Theta_b) / Sigma_b)) - 1. / (1. + exp(-Theta_b / Sigma_b));
}

//I_Ca
double I_Ca(double v){
	return Gca * pow(sinf(v), 2.0) * (v - Vca);
}

//sinf
double sinf(double v){
	return 1. / (1. + exp((v - Theta_s) / Sigma_s));
}

//I_AHP
double I_AHP(double v, double ca){
	return Gahp * (v - Vk) * (ca / (ca + K1));
}

//I_NMDA
double I_NMDA(double v, double mn, double hn){
	return Gnmda(v, mn, hn) * (v - Vn);
}

//I_DIC
double I_DIC(double v, double ld, double hd){
	return Gd * m_Dinf(ld) * hd * (v - Vd);
}

//m_Dinf
double m_Dinf(double ld){
	return pow(ld, 3.0) / (pow(ld, 3.0) + pow(Kd, 3.0));
}

//file output
void output(long ns, double t, STATE_stn x){
	output_state(ns, t, x);
	output_current(ns, t, x);
}

void output_state(long ns, double t, STATE_stn x){
	if(ns == 0){
		Fp_state = fopen("EDCS_out\\EDCS404_1\\state.dat", "w");
		fprintf(Fp_state, "t, t(sec), v, h, n, r, ca, hn, hd, ld\n");
	}
	
	fprintf(Fp_state, "%.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf\n", t, t/1000.0, x.v, x.h, x.n, x.r, x.ca, x.hn, x.hd, x.ld);
}

void output_current(long ns, double t, STATE_stn x){
	double ileak, iNa, iK, iT, iCa, iAHP, iNMDA, iDIC, iGPe;
	
	if(ns == 0){
		Fp_current = fopen("EDCS_out\\EDCS404_1\\current.dat", "w");
		fprintf(Fp_current, "t, t(sec), I_leak, I_Na, I_K, I_T, I_Ca, I_AHP, I_NMDA, I_DIC, I_GPe");
	}
	
	ileak = I_leak(x.v);
	iNa = I_Na(x.v, x.h);
	iK = I_K(x.v, x.n);
	iT = I_T(x.v, x.r);
	iCa = I_Ca(x.v);
	iAHP = I_AHP(x.v, x.ca);
	iNMDA = I_NMDA(x.v, Mn, x.hn);
	iDIC = I_DIC(x.v, x.ld, x.hd);
	iGPe = I_GPe();
	
	fprintf(Fp_current, "%.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf\n", ileak, iNa, iK, iT, iCa, iAHP, iNMDA, iDIC, iGPe);
}