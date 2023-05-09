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
#define		OUT_FOLDER		"out"

//data label
#define 	DLB_MAIN_BASE	"PRG"

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
#define		NUM_X_ALL			4					//total number of state quantities

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
	double hn;		//KR(4):h_N
	double hd;		//KR(6):h_D
	double ld;		//KR(5):[Ca]N
};

//variable
double T_step;
double T_init;
double T_fin;

double Gl;
double Gn;
double Gd;
double Vl;
double Vn;
double Vd;
double Vca;

double Tauhn;
double Thetahn;
double Sigmahn;
double Taud0, Taud1;
double Thetadt, Sigmadt;
double Phid;
double Thetad, Sigmad;

double Tauld;

double Eta, Mg, Gam;

double K1n;

double Kd, Nd;

double V_init;
double Hn_init;
double Hd_init;
double Ld_init;

double I_inj0;
double Mn0;

double Nstep_gamen;
double Fout_jikei;
double Tget_jikei_min;
double Tget_jikei_max;

double Nstep_jikei;

char Dlb_main[200];
char Time_text[200];

//function
void main_1kai();
void qs_main(void);
void set_dlb_main(int argc, char *argv[]);
void ruku_for_xa(double t, STATE_str &xa);
void ruku(double t, double x[], double h);
void c_x_dt(double t, double x_all[], double x_all_dt[]);
void voltage_dyn(double v, double hn, double hd, double ld, double *v_dt, double *hn_dt, 
 		double *hd_dt, double *ld_dt);
void xa_to_x_all(STATE_str xa, double x_all[]);
void x_all_to_xa(double x_all[], STATE_str &xa);
void init_QSM(STATE_str &xa);
void output_QSM(long ns, double t, STATE_str xa);
void output_state_jikei(long ns, double t, STATE_str xa);
void output_current_jikei(long ns, double t, STATE_str xa);
void my_fopen(FILE **fp, char fname[], long n_label);
double il(double v);
double inmda(double v, double mn, double hn);
double gnmda(double n, double mn, double hn);
double hninf(double ld);
double inmdaca(double v, double mn, double hn, double ld2);
double id(double v, double hd, double ld2);
double mdinf(double ld2);
double hdinf(double v);
double taud(double v);
void disp_start_main();
void disp_end_main();
void disp_start_end_time(char md[], char c[]);
double curr_all(double v, double hn, double hd, double ld, double mn);

//file pointer
FILE *Fp_state_jikei, *Fp_current_jikei;

int main(int argc, char *argv[]){
	T_step		= 0.05e0;						//[ms]
	T_init		= -1.e3;						//start time[ms]
	T_fin		= 10.e3;						//final time[ms]
	
	Gl=2.25; Gn=20.e0; Gd=20;					//peak conductance
	Vl=-60.e0; Vn=-20.e0; 						//reversal potential
	Vd=-18.e0; Vca=120.e0;
	
	Tauhn=3000.e0; Thetahn=0.7; Sigmahn=0.05;	//h_N
	Taud0=300.e0; Taud1=350.e0;					//h_DIC
	Thetadt=-60.e0; Sigmadt=3.e0;
	Phid=1.e0;
	Thetad=-95.e0; Sigmad=14.e0;
	
	Tauld=80.e0;								//[Ld]
	
	Eta=0.33; Mg=1.3; Gam=0.05;					//NMDA
	
	K1n=0.005;									//I_NMDA^CA
	
	Kd=1.15; Nd=3.0;							//m_inf^DIC
	
	//initial value
	V_init	= -71.4;	//v
	Hn_init = 0.388;	//hn
	Hd_init = 0.115;	//hd
	Ld_init = 0.613;	//ld
	
	I_inj0	= -90;)
	Mn0		= 1.e0;
	
	//display
	Nstep_gamen	= 200;
	
	//output
	Fout_jikei	= 2;
						
	Tget_jikei_min	= T_init;	//T_init <= Tget_jikei_min, Tget_jikei_max <= T_fin
	Tget_jikei_max	= T_fin;
	
	Nstep_jikei	= 20;
	
	set_dlb_main(argc, argv);
	
	disp_start_main();
    disp_start_end_time("s", Dlb_main);	//start time
    
    if(CALMODE == 0)main_1kai();
    
    disp_start_end_time("e", Dlb_main);	//end time
	disp_end_main();
}

//display
void disp_start_main(){
	printf("--start of %s--\n", Dlb_main);
}

void disp_end_main(){
	printf("--end of %s--\n", Dlb_main);
}

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

void main_1kai(){
	qs_main();
	printf("==================== I—¹ =====================\n");
}

//integral
void qs_main(void){
	double t, x_all[NUM_X_ALL];
	STATE_str xa;
	long ns;
	
	init_QSM(xa);
	for(ns = 0; ;ns++){
		t = T_init + ns * T_step;	//time
		if(ns % Nstep_gamen == 0) printf("t = %.3lf(sec)\n", t/1000.e0);
		
		output_QSM(ns, t, xa);		//file output
		ruku_for_xa(t, xa);			//Runge-Kutta
		
		if(t > T_fin - EPS6)return;	//final process
	}
}

//Runge-Kutta
void ruku_for_xa(double t, STATE_str &xa){
	double x_all[NUM_X_ALL];
	xa_to_x_all(xa, x_all);
	ruku(t, x_all, T_step);
	x_all_to_xa(x_all, xa);
}

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
	voltage_dyn(xa.v, xa.hn, xa.hd, xa.ld, &(xa_dt.v),  &(xa_dt.hn),  
 		 &(xa_dt.hd),  &(xa_dt.ld));
 	xa_to_x_all(xa_dt, x_all_dt);
}

//membrane potential dynamics
void voltage_dyn(double v, double hn, double hd, double ld, double *v_dt, double *hn_dt, 
 		double *hd_dt, double *ld_dt){
 	double i_inj;
 	double mn;											//KR(2):m_N)
 	i_inj = I_inj0;
 	mn = Mn0;
 	
	*v_dt  = -curr_all(v,hn,hd,ld,mn)+i_inj;			//KR(1)
	*hn_dt = (-hn+hninf(ld))/Tauhn;						//KR(4)
	*hd_dt = Phid*(hdinf(v)-hd)/taud(v);				//KR(7)
	*ld_dt = (-inmdaca(v,mn,hn)-ld)/Tauld;				//KR(5)
}

//xa¨x_all
void xa_to_x_all(STATE_str xa, double x_all[]){
	x_all[0] = xa.v;
	x_all[1] = xa.hn;
	x_all[2] = xa.hd;
	x_all[3] = xa.ld;
}

//x_all¨xa
void x_all_to_xa(double x_all[], STATE_str &xa){
	xa.v = x_all[0];
	xa.hn = x_all[1];
	xa.hd = x_all[2];
	xa.ld = x_all[3];
}

//initial
void init_QSM(STATE_str &xa){
	//xa
	xa.v	= V_init;
	xa.hn	= Hn_init;
	xa.hd	= Hd_init;
	xa.ld	= Ld_init;
}

//QSM
void output_QSM(long ns, double t, STATE_str xa){
	if(Fout_jikei >= 1){
		output_state_jikei(ns, t, xa);
		
		if(Fout_jikei == 2)output_current_jikei(ns, t, xa);
	}
}

//output
void output_state_jikei(long ns, double t, STATE_str xa){
	if(ns == 0){
		my_fopen(&Fp_state_jikei, "out_state_jikei", 1);
		fprintf(Fp_state_jikei, "Fp_state_jikei, "t, t(sec), v, hn, hd, ld\n");
	}
	
	if(t > Tget_jikei_min - EPS6 && t > Tget_jikei_max + EPS6 && ns % Nstep_jikei == 0){
		fprintf(Fp_state_jikei, "%.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf\n", t, t/1000.e0, xa.v, xa.hn, xa.hd, xa.ld);
	}
	
	if(abs(t - T_fin) < EPS6)fclose(Fp_state_jikei);
}

void output_current_jikei(long ns, double t, STATE_str xa){
	double il, id, inmda, inmdaca;
	double mn;
	double ld2;
	mn = Mn0;
	
	if(ns == 0){
		my_fopen(&Fp_current_jikei, "out_current_jikei", 1);
		fprintf(Fp_current_jikei, "t, t(sec), il, id, inmda, inmdaca\n");
	}
	
	if(t > Tget_jikei_min - EPS6 && t > Tget_jikei_max + EPS6 && ns % Nstep_jikei == 0){
		il = il(xa.v);
		inmda = inmda(xa.v, mn, xa.hn);
		inmdaca = inmdaca(xa.v, mn, xa.hn, ld2);
		id = id(xa.v, xa.hd, ld2);
		
		fprintf(Fp_current_jikei, "%.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf, %.7lf\n", t, t/1000.e0, v, il, id, inmda, inmdaca);
	}
	if(abs(t - T_fin) < EPS6)fclose(Fp_current_jikei);
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

//I_leak
double il(double v){
	return Gl * (v - Vl);
}

//I_NMDA
double inmda(double v, double mn, double hn){
	return gnmda(v, mn, hn) * (v - Vn);
}

//gnmda
double gnmda(double n, double, mn, double hn){
	return (Gn * mn * hn) / (1.e0 + Eta * Mg * exp(-Gam * v));
}

//h_N^inf
double hninf(double ld){
	return 1. / (1. + exp((ld-Thetahn)/Sigmahn));
}

//I_NMDA^Ca
double inmdaca(double v, double mn, double hn, double ld2){
	ld2 = K1n * gnmda(v,mn,hn) * (v - Vca);
	return ld2;
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

//I_DIC(KR(6))
double id(double v, double hd, double ld2){
	return Gd * mdinf(ld2) * hd * (v - Vd);
}

//m_d^inf
double mdinf(double ld2){
	return pow(ld2,Nd) / ( pow(ld2,Nd) + pow(Kd,Nd) );
}

//h_D^inf(v)
double hdinf(double v){
	return 1. / ( 1. + exp((v-Thetad)/Sigmad) );
}

//tau_D(v)
double taud(double v){
	return Taud0 + Taud1 / ( 1. + exp((v-Thetadt)/Sigmadt) );
}

//KR(1)
double curr_all(double v, double hn, double hd, double ld, double mn){
	return il(v) + id(v,hd,ld) + inmda(v, mn, hn); 
}