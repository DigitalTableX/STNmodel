#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <time.h>

#include "math.c"		//random number (0~1) generation routine
#include "time.c"		//return date and time

//output folder
#define		OUT_FOLDER	"out"

//data label
#define		DLB_MAIN_BASE		"PRG"

#define 	PAI					(atan(1.e0) * 4.e0)	//Pi
#define		EPS3				1.e-3
#define		EPS6				1.e-6
#define		EPS9				1.e-9
#define		EPS12				1.e-12
#define		BIG30				1.e30
#define		NBIG3				1000
#define		NBIG9				1000000000
#define		VERR				-2.e6
#define		NERR				-20000

//neuron model
#define		NUM_X_ALL			8					//total number of state quantities

#define		MAX_FILE_OPEN		500					//maximum number of files that can be open

//random number
#define		RANSUU_METHOD		1					//0: ran_fast_RCP
													//1: ran1_RCP
													//2: ran2_RCP
										
//mode
#define		CALMODE			0						//0:basic simulation(one time only)
									
#define		MAX(a,b)		((a) >= (b) ? (a) : (b))
#define		MIN(a,b)		((a) <= (b) ? (a) : (b))
#define		EQ(a,b)			(fabs((a)-(b))<EPS6 ? 1 : 0)//a=b:1
#define		NE(a,b)			(fabs((a)-(b))>EPS6 ? 1 : 0)//a!=b:1
#define		LE(a,b)			((a) < (b)+EPS6 ? 1 : 0)	//a<=b:1
#define		LT(a,b)			((a) < (b)-EPS6 ? 1 : 0)	//a<b:1
#define		GE(a,b)			((a) > (b)-EPS6 ? 1 : 0)	//a>=b:1
#define		GT(a,b)			((a) > (b)+EPS6 ? 1 : 0)	//a>b:1
#define		EQS(s1,s2)		(!strcmp((s1),(s2)))


struct STATE_str;

class STATE_str{
public:
	double v;		//KR(1):V
	double h;		//KR Appendix:h
	double n;		//KR Appendix:n
	double r;		//KR Appendix:r
	double ca;		//KR Appendix:[Ca]AHP
	double hn;		//KR(4):h_N
	double hd;		//KR(6):h_D
	double ld;		//KR(5):[Ca]N
};

//function
void check_fun();
void check_para_main();
void disp_start_main();
void disp_end_main();
void main_1kai();
void qs_main(void);
void ruku_for_xa(double t, STATE_str &xa);
void system_variable(long ns, double t, long n_post);
void c_fr_post(long ns, double t, long n_post);
void firing_process(long ns, double t, STATE_str xa, long *n_post);
long ch_hakka(long ns, double v);
void xa_to_x_all(STATE_str xa, double x_all[]);
void x_all_to_xa(double x_all[], STATE_str &xa);
void voltage_dyn(double v, double h, double n, double r, double ca, double hn, double hd, double ld, 
 		double *v_dt, double *h_dt, double *n_dt, double *r_dt, double *ca_dt, double *hn_dt, 
 		double *hd_dt, double *ld_dt);
double curr_all(double v, double h, double n, double r, double ca, double hn, double hd, double ld, double mn);
double curr_org(double v, double h, double n, double r, double ca);
double minf(double v);
double hinf(double v);
double ninf(double v);
double rinf(double v);
double ainf(double v);
double binf(double r);
double sinf(double v);
double tauh(double v);
double taun(double v);
double taur(double v);
double il(double v);
double ina(double v, double h);
double ik(double v, double n);
double iahp(double v, double ca);
double ica(double v);
double it(double v, double r);
double id(double v, double hd, double ld);
double mdinf(double ld);
double hdinf(double v);
double taud(double v);
double inmda(double v, double mn, double hn);
double gnmda(double v, double mn, double hn);
double hninf(double ld);
double inmdaca(double v, double mn, double hn);
void output_QSM(long ns, double t, STATE_str xa);
void output_state_jikei(long ns, double t, STATE_str xa);
void output_current_jikei(long ns, double t, STATE_str xa);
void init_QSM(STATE_str &xa);
void ruku(double t, double x[], double h);
void c_x_dt(double t, double x_all[], double x_all_dt[]);
void my_fopen(FILE **fp, char fname[], long n_label);
void set_dlb_main(int argc, char *argv[]);
float ransuu(void);
void ransuu_init(long i_seed);
void disp_start_end_time(char md[], char c[]);
void disp_err(char c[]);

//variable
double T_step, T_init, T_fin;

double Gl,Gna,Gk,Gahp,Gca,Gt;
double Vl,Vna,Vk,Vca;
double Thetam,Sigmam,Thetah,Sigmah;
double Thetan,Sigman,Thetar,Sigmar;
double Thetaa,Sigmaa,Thetab,Sigmab;
double Thetas,Sigmas;
	
double Tauh0,Tauh1,Thetaht,Sigmaht;
double Taun0,Taun1,Thetant,Sigmant;	
double Taur0,Taur1,Thetart,Sigmart;

double Phi,Phir;
double K1;
double Eps,Kca;

double Gn,Gd;
double Vn,Vd;
	
double Tauhn,Thetahn,Sigmahn;
double Taud0,Taud1,Thetadt,Sigmadt,Phid;
double Thetad,Sigmad;
	
double Tauld;
	
double Eta,Mg,Gam;
double K1n;
double Kd,Nd;

double V_init, H_init, N_init, R_init, Ca_init, Hn_init, Hd_init, Ld_init;

double V_th;
double I_inj0, Mn0;

long Nstep_gamen;

long Fout_jikei;
double Tget_jikei_min, Tget_jikei_max;
long Nstep_jikei;

long Fout_fr_post;										
double Tave_fr_post;

long Ransuu_seed;

//--QSM--
double Fr_post;			//firing rate[1/ms]
double Cv_post;			//interspike interval(ISI):CV

char Dlb_main[200];
char Time_text[200];

long I_ransuu;

//file pointer
FILE *Fp_state_jikei, *Fp_current_jikei;

int main(int argc, char *argv[]){
	T_step 		= 0.05e0;	//[ms]
	T_init		= -1.e3;	//[ms]
	T_fin		= 10.e3;	//[ms]

	//parameters for Ileak,INa,IK,IAHP,ICA,&IT
	Gl=2.25,Gna=37.5,Gk=45.,Gahp=20.,Gca=.5,Gt=0.1;	//peak conductance[mS/cm2]
	
	Vl=-60.,Vna=55.,Vk=-80.,Vca=120.;				//reversal potential[mV]
	
	Thetam=-30.,Sigmam=-14.5,Thetah=-39.,Sigmah=3.1;//X_inf(V)(X=m,h,n,r,a,b,s)
	Thetan=-32.,Sigman=-8.,Thetar=-67.,Sigmar=2.;
	Thetaa=-63.,Sigmaa=-7.8,Thetab=0.4,Sigmab=-0.1;
	Thetas=-39.,Sigmas=-8.;
	
	Tauh0=1.,Tauh1=500.,Thetaht=-85.,Sigmaht=3.;	//tau_X(V)(X=h,n,r)
	Taun0=1.,Taun1=100.,Thetant=-80.,Sigmant=26.;
	Taur0=40.,Taur1=17.5,Thetart=68.,Sigmart=2.2;
	
	Phi=0.75,Phir=0.2;								//phi_X(X=n,h,r)
	
	K1=7.5;											//I_AHP
	
	Eps=3.75e-5,Kca=22.5;							//Ca2+ dynamics

	//parameters for INMDA & IDIC
	Gn=20.; Gd=20.;				//I_NMDA,I_DIC:peak conductance[mS/cm2]
	Vn=-20.,Vd=-18.;			//I_NMDA,I_DIC:reversal potential[mV]
	
	Tauhn=3000.,Thetahn=0.7,Sigmahn=0.05;					//h_N
	Taud0=300.,Taud1=350.,Thetadt=-60.,Sigmadt=3.,Phid=1.;	//h_DIC
	Thetad=-95.,Sigmad=14.;
	
	Tauld=80.;						//[Ld]
	
	Eta=0.33,Mg=1.3,Gam=0.05;		//NMDA
	
	K1n=0.005;						//I_NMDA^CA
	
	Kd=1.15,Nd=3.0;					//m_inf^DIC

	//initial value
	V_init	=-71.4;			//v
	H_init	= 0.54;			//h
	N_init	= 0.0;			//n
	R_init 	= 0.664;		//r
	Ca_init	= 0.446;		//ca
	Hn_init	= 0.388;		//hn
	Hd_init	= 0.115;		//hd
	Ld_init	= 0.613;		//ld
	
	V_th	= -20.e0;		//threshold
	
	I_inj0	= -90.;
	Mn0		= 1.e0;			//KR(2):m_N
	
	Nstep_gamen			= 200;				//display
	
	//output
	//time series data
	Fout_jikei			= 2;

	Tget_jikei_min 		= T_init;					//T_init <= Tget_jikei_min, Tget_jikei_max <= T_fin
	Tget_jikei_max		= T_fin;
	
	Nstep_jikei			= 20;

	//post firing
	Fout_fr_post 		= 1;						//1:output
	Tave_fr_post   		= T_fin;					//ISI CV
	
	Ransuu_seed  = -1;
	
	ransuu_init(Ransuu_seed);
    set_dlb_main(argc, argv);
	check_para_main();
	
	disp_start_main();
    disp_start_end_time("s", Dlb_main);				//start time
	//check_fun();
	
	if (CALMODE == 0) main_1kai();
	
	disp_start_end_time("e", Dlb_main);				//end time
	disp_end_main();
}

//check function
void check_fun(){
	ransuu_init(Ransuu_seed);
	printf("ransuu_init\n");
	printf("ransuu() = %.7lf\n", ransuu());
	printf("ransuu() = %.7lf\n", ransuu());
	printf("ransuu() = %.7lf\n", ransuu());
	printf("ransuu() = %.7lf\n", ransuu());
	printf("ransuu() = %.7lf\n", ransuu());
	ransuu_init(Ransuu_seed);
	printf("ransuu_init\n");
	printf("ransuu() = %.7lf\n", ransuu());
	printf("ransuu() = %.7lf\n", ransuu());
	exit(1);
}

//check parameters
void check_para_main(){
	//Tget_jikei_min(max)
	if (Fout_jikei>= 1){
		if (!(GE(Tget_jikei_min, T_init) && LE(Tget_jikei_max, T_fin))) disp_err("err-CPM-TJKa");
	}
}

//display start
void disp_start_main(){
	printf("--start of %s--\n", Dlb_main);
}

//display end
void disp_end_main(){
	printf("--end of %s--\n", Dlb_main);
}

void main_1kai(){
	qs_main();
	printf("Fr_post = %.7lf[Hz]\n", Fr_post * 1000.e0);
	printf("Cv_post = %.7lf\n", Cv_post);
}

//integral
void qs_main(void){
	double t, x_all[NUM_X_ALL];
	STATE_str xa;
	long n_post;
	long ns;
	
	init_QSM(xa);
	for (ns = 0; ;ns++){
		t = T_init + ns * T_step; 			//time
		if (ns % Nstep_gamen == 0) printf("t = %.3lf(sec)\n", t/1000.e0);
		
		firing_process(ns, t, xa, &n_post);	//firing process
		system_variable(ns, t, n_post);		//system variable
		output_QSM(ns, t, xa);				//file output
		ruku_for_xa(t, xa);
		
		if (t > T_fin - EPS6) return;		//final process
	}
}

//Runge-Kutta
void ruku_for_xa(double t, STATE_str &xa){
	double x_all[NUM_X_ALL];
	xa_to_x_all(xa, x_all);
	ruku(t, x_all, T_step);
	x_all_to_xa(x_all, xa);
}

void system_variable(long ns, double t, long n_post){
	if (Fout_fr_post == 1) c_fr_post(ns, t, n_post);
}

//Fr_post),ISI CV(Cv_post)
void c_fr_post(long ns, double t, long n_post){
	static long n_sum;
	static double m1_sum, m2_sum;
	static double t_post_bef;			//[ms]
	double isi;							//ISI[ms]
	double isi_ave, isi_var, isi_sd;	//average, variance, standard deviation
	if (ns == 0){ n_sum = 0; m1_sum = 0.e0; m2_sum = 0.e0; t_post_bef = VERR; return;}
	
	if (t < T_fin - Tave_fr_post - EPS6) return;
	
	//m1_sum, m2_sum
	if (n_post == 1){					//firing process
		if (t_post_bef < VERR + EPS6){ 	//first firing
			t_post_bef = t;
		}
		else{							//after first firing
			isi = t - t_post_bef;
			m1_sum += isi;	
			m2_sum += isi * isi;
			n_sum++;
			t_post_bef = t;	
		}
	}
	
	//final time
	if ( abs(t - T_fin) < EPS6 ){
		if (n_sum == 0){
			Fr_post = 0.e0;
			Cv_post = -1.e0;
		}
		else{
			isi_ave	= m1_sum / n_sum;
			isi_var = m2_sum / n_sum - isi_ave * isi_ave;
			isi_sd  = sqrt(isi_var);
			Fr_post	= 1.e0 / isi_ave;
			Cv_post	= isi_sd / isi_ave;
		}
	}
}

//postsynaptic cell
void firing_process(long ns, double t, STATE_str xa, long *n_post){
	*n_post = ch_hakka(ns, xa.v);
}

//firing
long ch_hakka(long ns, double v){
	long nwk;
	static double v_bef;
	if (ns == 0){ v_bef = VERR; }
	
	if (v_bef < V_th && v >= V_th) nwk = 1;
	else nwk = 0;
	
	v_bef = v;
	
	return nwk;
} 

//xa¨x_all
void xa_to_x_all(STATE_str xa, double x_all[]){
	x_all[0] = xa.v;
	x_all[1] = xa.h;
	x_all[2] = xa.n;
	x_all[3] = xa.r;
	x_all[4] = xa.ca;
	x_all[5] = xa.hn;
	x_all[6] = xa.hd;
	x_all[7] = xa.ld;
}

//x_all¨xa
void x_all_to_xa(double x_all[], STATE_str &xa){
	xa.v = x_all[0];
	xa.h = x_all[1];
	xa.n = x_all[2];
	xa.r = x_all[3];
	xa.ca = x_all[4];
	xa.hn = x_all[5];
	xa.hd = x_all[6];
	xa.ld = x_all[7];
}

//membrane potential
void voltage_dyn(double v, double h, double n, double r, double ca, double hn, double hd, double ld, 
 		double *v_dt, double *h_dt, double *n_dt, double *r_dt, double *ca_dt, double *hn_dt, 
 		double *hd_dt, double *ld_dt){
 	double i_inj;
 	double mn;			//KR(2):m_N
 	i_inj = I_inj0;
 	mn = Mn0;
 	
	*v_dt  = -curr_all(v,h,n,r,ca,hn,hd,ld,mn)+i_inj;	//KR(1)
	*h_dt  = Phi*(hinf(v)-h)/tauh(v);					//KR Appendix(dX/dt = ¥¥, (X=n, h, r))
	*n_dt  = Phi*(ninf(v)-n)/taun(v);
	*r_dt  = Phir*(rinf(v)-r)/taur(v);
	*ca_dt = Eps*(-ica(v)-it(v,r)-Kca*ca);				//KR Appendix(d[Ca]_AHP/dt = ¥¥)
	*hn_dt = (-hn+hninf(ld))/Tauhn;						//KR(4)
	*hd_dt = Phid*(hdinf(v)-hd)/taud(v);				//KR(7)
	*ld_dt = (-inmdaca(v,mn,hn)-ld)/Tauld;				//KR(5)
}

//KR(1)
double curr_all(double v, double h, double n, double r, double ca, double hn, double hd, double ld, double mn){
	return curr_org(v,h,n,r,ca) + id(v,hd,ld) + inmda(v, mn, hn); 
}

double curr_org(double v, double h, double n, double r, double ca){
	return il(v)+ik(v,n)+ina(v,h)+it(v,r)+ica(v)+iahp(v,ca);
}

//X_inf(v)(X=m,h,n,r,a,s)(KR Appendix)
double minf(double v){
	return 1. / (1. + exp((v-Thetam)/Sigmam));
}

double hinf(double v){
	return 1. / (1. + exp((v-Thetah)/Sigmah));
}

double ninf(double v){
	return 1. / (1. + exp((v-Thetan)/Sigman));
}

double rinf(double v){
	return 1. / (1. + exp((v-Thetar)/Sigmar));
}

double ainf(double v){
	return 1. / (1. + exp((v-Thetaa)/Sigmaa));
}

//b_inf(r)(KR Appendix)
double binf(double r){
	return 1. / (1. + exp((r-Thetab)/Sigmab)) - 1. / ( 1.+exp(-Thetab/Sigmab) );
}

double sinf(double v){
	return 1. / ( 1. + exp((v-Thetas)/Sigmas) );
}

//tau_X(v)(X=h,n,r)(KR Appendix)
double tauh(double v){
	return Tauh0 + Tauh1 / ( 1. + exp((v-Thetaht)/Sigmaht) );
}

double taun(double v){
	return Taun0 + Taun1 / ( 1. + exp((v-Thetant)/Sigmant));
}

double taur(double v){
	return Taur0 + Taur1 / (1. + exp((v-Thetart)/Sigmart));
}

//Ileak
double il(double v){
	return Gl * (v-Vl);
}

//I_Na(KR Appendix)
double ina(double v, double h){
	return Gna * pow(minf(v),3.0) * h * (v-Vna);
}

//I_K(KR Appendix)
double ik(double v, double n){
	return Gk * pow(n,4.0) * (v-Vk);
}

//I_AHP(KR Appendix)
double iahp(double v, double ca){
	return Gahp * ca / (ca+K1) * (v-Vk);
}

//I_Ca(KR Appendix)
double ica(double v){
	return Gca * pow(sinf(v),2.0) * (v-Vca);
}

//I_T(KR Appendix)
double it(double v, double r){
	return Gt * pow(ainf(v),3.0) * pow(binf(r),2.0) * (v-Vca);
}

//I_DIC(KR(6))
double id(double v, double hd, double ld){
	return Gd * mdinf(ld) * hd * (v - Vd);
}

//m_d^inf([Ca]N)
double mdinf(double ld){
	return pow(ld,Nd) / ( pow(ld,Nd) + pow(Kd,Nd) );
}

//h_D^inf(v)
double hdinf(double v){
	return 1. / ( 1. + exp((v-Thetad)/Sigmad) );
}

//tau_D(v)
double taud(double v){
	return Taud0 + Taud1 / ( 1. + exp((v-Thetadt)/Sigmadt) );
}

//I_NMDA(KR(2))
double inmda(double v, double mn, double hn){
	return gnmda(v,mn,hn) * (v-Vn);
}

//G_NMDA(KR(3))
double gnmda(double v, double mn, double hn){
	return Gn * mn * hn / ( 1. + Eta*Mg*exp(-Gam*v) );
}

//h_N^inf
double hninf(double ld){
	return 1. / (1. + exp((ld-Thetahn)/Sigmahn));
}

//I_NMDA^Ca
double inmdaca(double v, double mn, double hn){
	return K1n * gnmda(v,mn,hn) * (v - Vca);
}

//QSM
void output_QSM(long ns, double t, STATE_str xa){
	if (Fout_jikei >= 1){
		output_state_jikei(ns, t, xa); 
		
		if (Fout_jikei == 2) output_current_jikei(ns, t, xa);
	}
} 

//data output
void output_state_jikei(long ns, double t, STATE_str xa){
	if (ns == 0){
		my_fopen(&Fp_state_jikei, "out_state_jikei", 1);
		fprintf(Fp_state_jikei, "t, t(sec), v, h, n, r, ca, hn, hd, ld\n");
	}
	
	if (t > Tget_jikei_min-EPS6 && t < Tget_jikei_max+EPS6 && ns % Nstep_jikei == 0){
		fprintf(Fp_state_jikei, "%.7le, %.7le, %.7le, %.7le, %.7le, %.7le, %.7le, %.7le, %.7le, %.7le\n", 
			t, t/1000.e0, xa.v, xa.h, xa.n, xa.r, xa.ca, xa.hn, xa.hd, xa.ld);
	}
	
	if (abs(t - T_fin) < EPS6) fclose(Fp_state_jikei);
}

//current output
void output_current_jikei(long ns, double t, STATE_str xa){
		double il_, ik_, ina_, it_, ica_, iahp_, curr_org_, id_, inmda_, inmdaca_;
	double mn;
	mn = Mn0;
	
	if (ns == 0){
		my_fopen(&Fp_current_jikei, "out_current_jikei", 1);
		fprintf(Fp_current_jikei, "t, t(sec), il, ik, ina, it, ica, iahp, curr_org, id, inmda, inmdaca\n");
	}
	
	if (t > Tget_jikei_min-EPS6 && t < Tget_jikei_max+EPS6 && ns%Nstep_jikei == 0){
		il_ = il(xa.v); ik_ = ik(xa.v, xa.n); ina_ = ina(xa.v, xa.h); it_ = it(xa.v, xa.r); ica_ = ica(xa.v);
		iahp_ = iahp(xa.v, xa.ca);
		curr_org_ = curr_org(xa.v, xa.h, xa.n, xa.r, xa.ca);
		id_ = id(xa.v, xa.hd, xa.ld); inmda_ = inmda(xa.v, mn, xa.hn);
		inmdaca_ = inmdaca(xa.v, mn, xa.hn);
	
		fprintf(Fp_current_jikei, "%.7le, %.7le, %.7le, %.7le, %.7le, %.7le, %.7le, %.7le, %.7le, %.7le, %.7le, %.7le\n", 
			t, t/1000.e0, il_, ik_, ina_, it_, ica_, iahp_, curr_org_, id_, inmda_, inmdaca_);
	}
	
	if (abs(t - T_fin) < EPS6) fclose(Fp_current_jikei);
	
}

//initial
void init_QSM(STATE_str &xa){
	//xa
	xa.v = V_init;
	xa.h = H_init;
	xa.n = N_init;
	xa.r = R_init;
	xa.ca = Ca_init;
	xa.hn = Hn_init;
	xa.hd = Hd_init;
	xa.ld = Ld_init;
	
	//Fr_post, Cv_post
	Fr_post = VERR;
	Cv_post = VERR;
}

//Runge-Kutta
void ruku(double t, double x[], double h){
	long i;
	double k1[NUM_X_ALL], k2[NUM_X_ALL], k3[NUM_X_ALL], k4[NUM_X_ALL];
	double wk1[NUM_X_ALL];
	
	//1
	c_x_dt(t, x, k1);
	
	//2
	for (i = 0; i < NUM_X_ALL; i++) wk1[i] = x[i] + h * k1[i]/2.e0;
	c_x_dt(t + h/2.e0, wk1, k2);
	
	//3
	for (i = 0; i < NUM_X_ALL; i++) wk1[i] = x[i] + h * k2[i]/2.e0;
	c_x_dt(t + h/2.e0, wk1, k3);
	
	//4
	for (i = 0; i < NUM_X_ALL; i++) wk1[i] = x[i] + h * k3[i];
	c_x_dt(t + h, wk1, k4);

	//update
	for (i = 0; i < NUM_X_ALL; i++){
		x[i] = x[i] + h * (k1[i]/6.e0 + k2[i]/3.e0 + k3[i]/3.e0 + k4[i]/6.e0);
	}
}

void c_x_dt(double t, double x_all[], double x_all_dt[]){
	STATE_str xa, xa_dt;
	x_all_to_xa(x_all, xa);
	voltage_dyn(xa.v, xa.h, xa.n, xa.r, xa.ca, xa.hn, xa.hd, xa.ld, 
 		 &(xa_dt.v),  &(xa_dt.h),  &(xa_dt.n),  &(xa_dt.r),  &(xa_dt.ca),  &(xa_dt.hn),  
 		 &(xa_dt.hd),  &(xa_dt.ld));
 	xa_to_x_all(xa_dt, x_all_dt);
}

//file
void my_fopen(FILE **fp, char fname[], long n_label){
	static long n_fop = 0;
	char text_wk1[200], text_folder[200];
	long i;	
	
	if (n_fop >= MAX_FILE_OPEN){ disp_err("err in my-fopen for opening too many files\n"); }
	if (n_label >= 2) disp_err("err my_fopen n-label\n");
	
	//file name
	sprintf(text_wk1, "%s\\%s", OUT_FOLDER, fname);
	if (n_label == 1) strcat(text_wk1, Dlb_main);
	strcat(text_wk1, ".dat");
		
	if ( (*fp = fopen(text_wk1, "w") ) == NULL ){ 
		printf("err in fopcl for write mode open %s", text_wk1); 
		exit(1);
	}
	n_fop++;
}

//data label
void set_dlb_main(int argc, char *argv[]){
	if ( !strcmp(DLB_MAIN_BASE, "PRG") ){
		sprintf(Dlb_main, "(%s", argv[0]);
	}
	else{
		sprintf(Dlb_main, "(%s", DLB_MAIN_BASE); //DLB_MAIN_BASE
	}
	
	strcat(Dlb_main, ")");
}

//random number
float ransuu(void){
	if (RANSUU_METHOD == 2){
		return ran2_RCP(&I_ransuu);
	}
	else if (RANSUU_METHOD == 1){
		return ran1_RCP(&I_ransuu);
	}
	else if (RANSUU_METHOD == 0){
		return ran_fast_RCP(&I_ransuu);
	}
	else{
		printf("err in ran-suu\n"); 
		exit(1);
	}
}

void ransuu_init(long i_seed){
	if (i_seed >= 0){ printf("err i_seed >= 0\n"); exit(1); }
	I_ransuu = i_seed;
}

//time output
void disp_start_end_time(char md[], char c[]){
	FILE *fp;
	if ( (fp = fopen("time_now.dat", "a") ) == NULL){printf("err DPSET\n");exit(1); }
	
	c_time_now(Time_text);
	if (!strcmp(md, "s")){
		fprintf(fp, "start time of %s = %s\n", c, Time_text);
	}
	else{
		fprintf(fp, "end time of %s = %s\n", c, Time_text);
	}
	fclose(fp);
}

//error message
void disp_err(char c[]){
	printf("%s\n", c);
	exit(1);
}