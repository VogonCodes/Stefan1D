#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <thread>
#include <functional>
#include <future>
#include <cstdio>

using namespace std;

using Matrix = vector<vector<double>>;
using Vector = vector<double>;

// Define constants
const double k = 2.0; // J/(m*s*K)
const double c = 4.0e3; // J*(kg*K)
const double rho = 1.0e3; // kg/m^3
const double L = 320.0e3; // J/kg

const double Ttop = 253.15; // K; = -20 C
const double Tbot = 273.15; // K; =  0C
const double Tmelt = 273.15; // K; =  0C

void printVector(Vector vector, int size){
	for(const auto& value : vector){
		cout << value << " ";
	}
	cout << endl;
}

void printMatrix(Matrix matrix, int rows, int cols){
	for(const auto& value : matrix){
		for (auto v: value){
			cout << v << " ";
		}
		cout << endl;
	}
	cout << endl;
}

// Useful stuff

void dumpTempVector(string filepath, Vector& vector, int size, double dx) {
	ofstream file(filepath);
	for (int i=0; i<size; i++) {
		file << i*dx << " " << vector[i] << endl;
	}
}

void forwardEuler2(Vector& y, Vector& y_old, size_t N, double a) {
	// y(x_i, t_{n+1}) = y(x_{i+1}, t_n - delta t) + a * [ y(x_{i+1}, t_{n-1}) - 2*y(x_i, t_{n-1}) + y(x_{i-1}, t_{n-1}) ]
	for (int i=1; i<N-1; i++) {
		y[i] = y_old[i] + a*(y_old[i+1] - 2*y_old[i] + y_old[i-1]);
	}
}

double findLambda(double dT, double latentHeat, double heatCap, double initGuess, double TOL) {
	double RHS = latentHeat * sqrt(M_PI) / (heatCap * dT);
	int i = 0;

	double lambda = initGuess, lambda_old = lambda + TOL;
	// Newton
	double f,fprime;
	double ee, erflm;

	f = lambda - RHS;
	//cout << i << " " << lambda << " " << lambda_old << " " << f << " " << fprime <<  endl;
	while (abs(f) > TOL) {
		ee = exp(-lambda*lambda);
		erflm = erf(lambda);

		f = ee/lambda/erflm - RHS;
		fprime = -ee/lambda/lambda/erflm - 2*ee/erflm - 2*ee*ee/sqrt(M_PI)/lambda/erflm/erflm;

		lambda_old = lambda;
		lambda -= f/fprime;

		//i++;
		//if (i<10) cout << i << " " << lambda << " " << lambda_old << " " << f << " " << fprime <<  endl;
		//if (i % 10000000 ==0) cout << i << " " << lambda << " " << lambda_old << " " << f << " " << fprime <<  endl;
	}

	return lambda;
}

void StefanProblem1D(Vector& numSol, Vector& analSol, Vector Tinit, size_t& numSolSize, size_t& analSolSize, double maxTime, double& dt, double& dx, double ymInit) {
	double time = dt;
	double ym = ymInit, ymOld = ymInit;
	double kappa = k/rho/c;
	double lambda = findLambda(20, L, c, 0.3, 1e-15);
	double maxdx = 0.005;
	double dxOld = dx;
	dx = 0.001;

	Vector oldSol(numSolSize);
	Vector tmp;
	int numSolSizeOld = numSolSize;
	oldSol = Tinit;
	numSol = Tinit;

	// Numerical solution + solidification interface depth
	while (time < maxTime) {
		ymOld = ym;
		oldSol = numSol;
		dxOld = dx;
		numSolSizeOld = numSolSize;
		// update solidification interface depth
		ym = ymOld + lambda * sqrt(kappa/time)*dt;


		// solve HEq for current time
		forwardEuler2(numSol, oldSol, numSolSize, kappa*dt/dx/dx);
		
		time += dt;
	}

	// Analytical solution
	double dT = Tmelt-Ttop;
	double erfLambda = erf(lambda);
	double sqr = sqrt(kappa*time);
	analSol.resize(numSolSize);
	analSolSize = numSolSize;
	for (int i=0; i<numSolSize; i++) {
		// for y<ym: theta = erf(eta)/erf(lambda) -> T = T0 + (Tm-T0)*erf(y/(2*sqrt(kappa*t)))/erf(lambda)
		if (ym > i*dx) analSol[i] = Ttop + dT*erf(i*dx/2/sqr)/erfLambda;
		else analSol[i] = Tmelt;
	}
}

// Useful stuff
int main() {
	size_t N = 3;
	size_t M = 0;
	double dx = 0.001;
	double dt = 0.1;
	double ymInit = N*dx;
	Vector solution(N,0.0), solAn;
	// Linear initial condition
	Vector Tinit(N, 273.15);
	for (int i=0; i<N; i++) Tinit[i] = Ttop + i*(Tmelt-Ttop)/M;

	StefanProblem1D(solution, solAn, Tinit, N, M, 5.0, dt, dx, ymInit);
	printVector(solAn,M);
	printVector(solution, N);
	dumpTempVector("initialCond.dat", Tinit, N, 0.001);
	dumpTempVector("numSol.dat", solution, N, 0.001);
	dumpTempVector("analSol.dat", solAn, N, 0.001);
}
