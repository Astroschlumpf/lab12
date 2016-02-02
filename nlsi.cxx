#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
#include <cmath>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double t);
void stepL(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int N);
void stepN(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int N);

//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double eta = 0.2;
        double t=0;
	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;

	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0", dx,Nx,xmin,t);
// 	int i = 1;

	for (int i = 1; i <= Na; i++) {

		for (int j = 1; j <= Nk-1; j++) {
		  stepL(psi1,psi0,dt,dx,Nx);
		  stepN(psi1,psi0,dt,dx,Nx);
		  h=psi0;
		  psi0=psi1;
		  psi1=h;
		  t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin,t);
	}

	return 0;
}

//-----------------------------------
void stepL(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int N)
{

  cmplx* d=new cmplx[N];
  cmplx* u=new cmplx[N];
  cmplx* l=new cmplx[N];
 
  cmplx alpha = - cmplx(0.,dt/(dx*dx));
  
  for(int i=0;i<N;i++) d[i] = cmplx(1.0, -2.0*dt/(dx*dx));
  for(int i=0;i<N;i++) u[i] = -alpha;
  for(int i=0;i<N;i++) l[i] = -alpha;
  
//   for(int i=1;i<N;i++){
//     d[i]=d[i]*d[i-1]-u[i-1]*l[i];
//     u[i]*=d[i-1];
//     f0[i]=f0[i]*d[i-1]-f0[i-1]*l[i];
//   }
//   f1[N-1]=f0[N-1]/d[N-1];
//   for (int i=N-2;i>=0;i--){
//     f1[i]=(f0[i]-u[i]*f1[i+1])/d[i];
//   }


  
  u[0] = u[0] / d[0]; 
  f0[0] = f0[0] / d[0];
  d[0] = 1.0;
  for(int i=1;i<N;i++){
    cmplx faktor = 1.0/(d[i] -  (-alpha)/d[i-1]*u[i-1]);
    d[i]  = 1;
       
    u[i] = u[i] * faktor;
    f0[i] = ( f0[i] - (-alpha)/d[i-1]*f0[i-1] )  * faktor;
  }
  
  
   f1[N-1]=f0[N-1];
    for (int i=N-2;i>=0;i--){
      f1[i]=f0[i]-u[i]*f1[i+1];
   }
  
  delete[] d;
  delete[] u;
  delete[] l;
}

//-----------------------------------
void stepN(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int N)
{
  for(int i = 0; i < N; i++){
    f1[i] = f1[i] * exp(cmplx(0., - norm(f1[i]) * dt));
  }
}

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double t)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx + t*5;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}

//-----------------------------------
void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 2*f/cosh(eta * x);
	}
}
