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

void StefanAnal1D(Vector& solution, size_t size, double dx, double time, double etam) {
	double deltaT = Tmelt - Ttop;
	double erfEtam = erf(etam);
	double eta;
	double kappa = k/rho/c;
	for (int i=0; i<size; i++) {
		// for y<ym: theta = erf(eta)/erf(lambda) -> T = T0 + (Tm-T0)*erf(y/(2*sqrt(kappa*t)))/erf(lambda)
		eta = i*dx / 2 / sqrt(kappa*time);
		solution[i] = Ttop + deltaT*erf(eta)/erf(etam);
	}
}

void StefanProblem1D(Vector& numSol, Vector& analSol, Vector Tinit, size_t& numSolSize, size_t& analSolSize, double maxTime, double& dt, double& dx, double ymInit, double lambdaAnal) {
	double time = dt;
	double ym = ymInit, ymOld = ymInit;
	double kappa = k/rho/c;
	double maxdx = 0.005;
	double dxOld = dx;
	double CFL = 0.5;
	//double lambda = ymInit / 2 / sqrt(kappa * time);
	double lambda = lambdaAnal;

	Vector oldSol(numSolSize);
	Vector tmp;
	int numSolSizeOld = numSolSize;
	oldSol = Tinit;
	numSol = Tinit;

	// Numerical solution + solidification interface depth
	int j = 0;
	while (time < maxTime) {
		// update old values
		ymOld = ym;
		oldSol = numSol;
		dxOld = dx;
		numSolSizeOld = numSolSize;

		// update solidification interface depth
		ym = ymOld + lambda * sqrt(kappa/time)*dt;
		//lambda = ym / 2 / sqrt(kappa * time);

		/* update grid :( */

		double TOL = 1e-14;
		if (ym/numSolSize > maxdx) numSolSize++;
		// update spatial step -- equidistant grid
		dx = ym/numSolSize;
		// resize the vector
		tmp.resize(numSolSize);
		tmp[numSolSize-1] = Tmelt;
		tmp[0] = Ttop;
		// linear interpolation to new grid
		for (int i=1; i<numSolSize-1; i++) {
			double xtilde = i*dx;
			if (xtilde/dxOld < TOL) tmp[i] = oldSol[i];
			else {
				int ix0 = static_cast<int>(xtilde/dxOld);
				if (ix0 >= (numSolSizeOld - 1)) tmp[i] = 273.15;
				else {
					int ix1 = ix0 + 1;
					if (ix1 >= (numSolSizeOld-1)) ix1 = numSolSizeOld - 1;
					// Ttilde = (T0*(x1 - xtilde) + T1*(xtilde - x0)) / (x1 - x0)
					tmp[i] = (oldSol[ix0]*(ix1*dxOld - xtilde) + oldSol[ix1]*(xtilde - ix0*dxOld)) / ((ix1-ix0)*dxOld);
				}
			}
		}
		oldSol.resize(numSolSize);
		numSol.resize(numSolSize);
		oldSol = tmp;


		// solve HEq for current time
		forwardEuler2(numSol, oldSol, numSolSize, kappa*dt/dx/dx);
		dumpTempVector("tmp/dataset_"+to_string(j)+".dat", numSol, numSolSize, dx);
		// CFL condition
		if (dt > dx*dx/2/kappa) dt = dx*dx/2/kappa;
		time += dt;
		j++;
	}

	// Analytical solution
	double dT = Tmelt-Ttop;
	analSol.resize(numSolSize);
	analSolSize = numSolSize;
	double eta, etam = ym / 2 / sqrt(kappa*time);
	//StefanAnal1D(analSol, analSolSize, dx, time, ym/2/kappa*time);
	for (int i=0; i<numSolSize; i++) {
		// for y<ym: theta = erf(eta)/erf(lambda) -> T = T0 + (Tm-T0)*erf(y/(2*sqrt(kappa*t)))/erf(lambda)
		eta = i*dx / 2 / sqrt(kappa*time);
		analSol[i] = Ttop + dT*erf(eta)/erf(etam);
	}
}

int main() {
	size_t N = 30;
	size_t M = 0;
	double dx = 0.0001;
	double ymInit = N*dx;
	double lambda = findLambda(Tmelt-Ttop, L, c, 0.5, 1e-15);
	double dt = ymInit*ymInit/4/lambda/lambda/(k/rho/c);;
	Vector solution(N,0.0), solAn;
	// Linear initial condition
	Vector Tinit(N, 273.15);
	//for (int i=0; i<N-1; i++) Tinit[i] = Ttop + i*(Tmelt-Ttop)/(N-1);
	//StefanAnal1D(Tinit, N, dxInit, dt, ymInit);
	// double lambda = findLambda(20.0, L, c, 0.3, 1e-15);
	for (int i=0;i<N-1; i++) Tinit[i] = Ttop + (Tmelt - Ttop)*erf(i*dx*sqrt(k/rho/c * dt))/erf(lambda);

	dumpTempVector("initialCond.dat", Tinit, N, dx);

	StefanProblem1D(solution, solAn, Tinit, N, M, 1000.0, dt, dx, ymInit, lambda);
	printVector(solAn,M);
	printVector(solution, N);
	dumpTempVector("numSol.dat", solution, N, dx);
	dumpTempVector("analSol.dat", solAn, N, dx);
}
