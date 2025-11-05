#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>

using namespace std;

using Matrix = vector<vector<double>>;
using Vector = vector<double>;

// Define constants
const double k = 2.0; // J/(m*s*K)
const double c = 4.0e3; // J*(kg*K)
const double rho = 1.0e3; // kg/m^3
const double L = 320.0e3; // J/kg
const double kappa = k/c/rho;

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

void dumpMatrix(string filepath, Matrix& matrix, size_t rows, size_t cols) {
	ofstream file(filepath);
	for (int i=0; i<rows; i++) {
		for (int j=0; j<cols; j++) {
			file << matrix[i][j] << " ";
		}
		file << endl;
	}
}

void dumpTempVector(string filepath, Vector& vector, int size, double dx) {
	ofstream file(filepath);
	for (int i=0; i<size; i++) {
		file << i*dx << " " << vector[i] << endl;
	}
}

void centralDifferences(Vector& y, Vector& y_old, size_t N, double a) {
	// y(x_i, t_{n+1}) = y(x_{i+1}, t_n - delta t) + a * [ y(x_{i+1}, t_{n-1}) - 2*y(x_i, t_{n-1}) + y(x_{i-1}, t_{n-1}) ]
	y[0] = y_old[0];
	y[N-1] = y_old[N-1];
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
	while (abs(f) > TOL) {
		ee = exp(-lambda*lambda);
		erflm = erf(lambda);

		f = ee/lambda/erflm - RHS;
		fprime = -ee/lambda/lambda/erflm - 2*ee/erflm - 2*ee*ee/sqrt(M_PI)/lambda/erflm/erflm;

		lambda_old = lambda;
		lambda -= f/fprime;
	}

	return lambda;
}

void StefanAnal1D(Vector& solution, size_t size, double dx, double ym, double time, double lambdaAnal) {
	double deltaT = Tmelt - Ttop;
	double erfEtam = erf(lambdaAnal);
	double eta;
	cout << deltaT << "," << erfEtam << "," << Ttop << "," << time << endl;
	for (int i=0; i<size; i++) {
		// for y<ym: theta = erf(eta)/erf(lambda) -> T = T0 + (Tm-T0)*erf(y/(2*sqrt(kappa*t)))/erf(lambda)
		if (i*dx > ym) solution[i] = Tmelt;
		eta = i*dx / 2 / sqrt(kappa*time);
		cout << eta << "," << erf(eta) << "," << erf(eta)/erfEtam << "," << deltaT*erf(eta)/erfEtam << endl;
		solution[i] = Ttop + deltaT*erf(eta)/erfEtam;
	}
}

void HeatEqAnal(Vector& solution, Vector Tinit, size_t size, double time, double dx) {
	// tepelne jadro: G(x,t) = 1/sqrt(4*pi*kappa*t) * exp(-x^2/(4*kappa*t))*Heaviside(t)
	// analyticke reseni: u(x,t) = convolution(G,u0)
	int padding = 2*size;
	size_t greenSize = 2*padding + 1;
	size_t fullSize = greenSize + size - 1;
	double x = 0.0;
	Vector Green(greenSize, 0.0);
	Vector fullSol(fullSize, 0.0);
	for (int i=0; i<greenSize; i++) {
		// cout << i << ": " << (i-static_cast<int>(size)) << " " << (i-size)*dx << " " << exp(-x*x/4/time) << " " << sqrt(4*M_PI*time) << endl;
		x = (i - padding) * dx;
		Green[i] = exp(-x*x/4/kappa/time) / sqrt(4*M_PI*kappa*time);
	}
	for (int i=0; i<size; i++) {
		for (int j=0; j<greenSize; j++) {
			fullSol[i+j] += Tinit[i]*Green[j]*dx;
		}
	}
	for (int i = 0; i< size; i++) solution[i] = fullSol[i+padding];
	dumpTempVector("green.dat", Green, greenSize, dx);
	dumpTempVector("fullSol.dat", fullSol, fullSize, dx);
}

void HeatEq1D(Vector& numSol, Vector Tinit, size_t numSolSize, double maxTime, double& time, double dx) {
	double dt = dx*dx/2/kappa;
	Vector oldSol(numSolSize);
	oldSol = Tinit;
	numSol = Tinit;

	while ((time+dt) >= maxTime) dt /= 10;

	time += dt;
	int j=0;
	while (time < maxTime) {
		// update old values
		oldSol = numSol;

		// solve HEq for current time
		//dumpTempVector("tmp/dataset_"+to_string(j)+".dat", numSol, numSolSize, dx);
		centralDifferences(numSol, oldSol, numSolSize, kappa*dt/dx/dx);

		dumpTempVector("tmp/dataset_"+to_string(j)+".dat", numSol, numSolSize, dx);

		// update time
		time += dt;
		j++;
	}
	cout << j << endl;
}

int main() {
	size_t N = 100;
	size_t M = 0;
	double dx = 0.01;

	double time = 0.0;

	Vector solution(N,0.0), solAn(M);
	// Initial temperature field
	Vector Tinit(N, Tmelt);
	Tinit[0] = Ttop;
	// /// Linear initial condition
	// for (int i=0; i<N-1; i++) Tinit[i] = Ttop + i*(Tmelt-Ttop)/(N-1);
	/// erf initial temp
	// double lambda = findLambda(Tmelt-Ttop, L, c, 0.5, 1e-15);
	// double time = N*dx*N*dx/4/lambda/lambda/(kappa);;
	// StefanAnal1D(Tinit, N, dx, N*dx, time, lambda);

	dumpTempVector("initialCond.dat", Tinit, N, dx);

	cout << time << endl;
	HeatEq1D(solution, Tinit, N, 5000.0, time, dx);
	//printVector(solAn,M);
	printVector(solution, N);
	dumpTempVector("numSol.dat", solution, N, dx);
	HeatEqAnal(solution, Tinit, N, 5000.0, dx);
	dumpTempVector("analSol.dat", solution, N, dx);
}
