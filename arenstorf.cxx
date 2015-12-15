#include <iostream>
#include <math.h>
#include <cmath>

using namespace std;


void f(double* v, double mu){
	
	// save current values
	double x = v[0];
	double xd = v[1];
	double y = v[2];
	double yd = v[3];
	double z = v[4];
	
	double s = sqrt((x - 1 + mu)*(x - 1 + mu) + y*y + z*z);
	double r = sqrt((x + mu)*(x + mu) + y*y + z*z);

	v[0] = xd; // x dot
	v[1] = x + 2*yd - (1 - mu)*(x + mu)/(r*r*r) - mu*(x - 1 + mu)/(s*s*s); // x dot dot
	v[2] = yd; // y dot
	v[3] = y - 2*xd - (1 - mu)*y/(r*r*r) - mu*y/(s*s*s); // y dot dot
	v[4] = v[5]; // z dot
	v[5] = -(1 - mu)*z/(r*r*r) - mu*z/(s*s*s); // z dot dot
	}



void fill(double* Ki,double* Kv, double* yn,  double dt){
	for(int i = 0; i<6 ; i++){
		Ki[i] = yn[i] + dt*Kv[i];
	
	}
}



int main(){

	double dt = 1e-3;
	double Xi;
	double error = 1e-5;

	double yn[6]; // yn
	double yn14[6]; // yn+1 Rk4
	double yn15[6]; // yn+1 Rk5
	double K1[6];
	double K2[6];
	double K3[6];
	double K4[6];
	double K5[6];
	double Kzw[6];

	double b4[4];
	double b5[5];

	// coefficients
	b4[0] = 35.0/384.0;
	b4[1] = 0.0;
	b4[2] = 500.0/1113.0;
	b4[3] = 125.0/192.0;

	b5[0] = 5179.0/57600.0;
	b5[1] = 0;
	b5[2] = 7571.0/16695.0;
	b5[3] = 393.0/640.0;
	b5[4] = -92097.0/339200.0;
	
	double mu = 0.012277471;

	// initial values

	yn[0] = 0.994; // x0
	yn[1] = 0; // x dot
	yn[2] = 0; // y0
	yn[3] = -2.00158510637908; // y dot
	yn[4] = 0; // z0
	yn[5] = 0; // z dot

	// output for initial values
	cout << 0 << "\t"<< yn[0] <<"\t"<< yn[2] << endl;

  for(double t=dt;t<2;t+=dt){

	// Calculating K1
	K1[0] = yn[0];
	K1[1] = yn[1];
	K1[2] = yn[2];
	K1[3] = yn[3];
	K1[4] = yn[4];
	K1[5] = yn[5];

	f(K1, mu); 

	// Calculating K2
	fill(K2,K1,yn,0.2*dt);
	f(K2,mu);

	// Calculating K3
	Kzw[0] = (3.0/40.0)*K1[0] + (9.0/40.0)*K2[0];
	Kzw[1] = (3.0/40.0)*K1[1] + (9.0/40.0)*K2[1];
	Kzw[2] = (3.0/40.0)*K1[2] + (9.0/40.0)*K2[2];
	Kzw[3] = (3.0/40.0)*K1[3] + (9.0/40.0)*K2[3];
	Kzw[4] = (3.0/40.0)*K1[4] + (9.0/40.0)*K2[4];
	Kzw[5] = (3.0/40.0)*K1[5] + (9.0/40.0)*K2[5];
	
	fill(K3,Kzw,yn,dt);
	f(K3,mu);

	// Calculating K4
	Kzw[0] = (44.0/45.0)*K1[0] + (-56.0/15.0)*K2[0] + (32.0/9.0)*K3[0];
	Kzw[1] = (44.0/45.0)*K1[1] + (-56.0/15.0)*K2[1] + (32.0/9.0)*K3[1];
	Kzw[2] = (44.0/45.0)*K1[2] + (-56.0/15.0)*K2[2] + (32.0/9.0)*K3[2];
	Kzw[3] = (44.0/45.0)*K1[3] + (-56.0/15.0)*K2[3] + (32.0/9.0)*K3[3];
	Kzw[4] = (44.0/45.0)*K1[4] + (-56.0/15.0)*K2[4] + (32.0/9.0)*K3[4];
	Kzw[5] = (44.0/45.0)*K1[5] + (-56.0/15.0)*K2[5] + (32.0/9.0)*K3[5];
	
	fill(K4,Kzw,yn,dt);
	f(K4,mu);

	// Calculating K5
	Kzw[0] = (19372.0/6561.0)*K1[0] + (-25360.0/2187.0)*K2[0] + (64448.0/6561.0)*K3[0] + (-212.0/729.0)*K4[0];
	Kzw[1] = (19372.0/6561.0)*K1[1] + (-25360.0/2187.0)*K2[1] + (64448.0/6561.0)*K3[1] + (-212.0/729.0)*K4[1];
	Kzw[2] = (19372.0/6561.0)*K1[2] + (-25360.0/2187.0)*K2[2] + (64448.0/6561.0)*K3[2] + (-212.0/729.0)*K4[2];
	Kzw[3] = (19372.0/6561.0)*K1[3] + (-25360.0/2187.0)*K2[3] + (64448.0/6561.0)*K3[3] + (-212.0/729.0)*K4[3];
	Kzw[4] = (19372.0/6561.0)*K1[4] + (-25360.0/2187.0)*K2[4] + (64448.0/6561.0)*K3[4] + (-212.0/729.0)*K4[4];
	Kzw[5] = (19372.0/6561.0)*K1[5] + (-25360.0/2187.0)*K2[5] + (64448.0/6561.0)*K3[5] + (-212.0/729.0)*K4[5];
	
	fill(K5,Kzw,yn,dt);
	f(K5,mu);

	// yn+1 for RK4 method
	for(int i=0;i<6;i++){
		yn14[i] = yn[i] + dt*(b4[0]*K1[i] + b4[1]*K2[i] + b4[2]*K3[i] + b4[3]*K4[i]);
	}
	
	// yn+1 for RK5 method
	for(int i=0;i<6;i++){
		yn15[i] = yn[i] + dt*(b5[0]*K1[i] + b5[1]*K2[i] + b5[2]*K3[i] + b5[3]*K4[i] + b5[4]*K5[i]);
	}

	// output
	//cout << t << "\t"<< yn14[0] <<"\t"<< yn14[2] << endl;	
	cout << t << "\t" << dt << endl;

	// Step-size control
	Xi = abs(yn14[0] - yn15[0]) + abs(yn14[1] - yn15[1]) + abs(yn14[2] - yn15[2]) + abs(yn14[3] - yn15[3]) + abs(yn14[4] - yn15[4]) + abs(yn14[5] - yn15[5]);

	dt = 0.5*dt*pow(error/Xi,1.0/5.0);

	for(int i=0; i<6;i++) {
		yn[i] = yn14[i];
	}

}
}
