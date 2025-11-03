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

bool updateGrid(Vector& solution, size_t& numSolSize, double ymCurrent, double& dx, double maxdx) {
	size_t numSolSizeOld = numSolSize;
	double dxOld = dx;
	double TOL = 1e-14;
	Vector tmp;

	if (ymCurrent/numSolSize > maxdx) numSolSize++;
	// update spatial step -- equidistant grid
	dx = ymCurrent/numSolSize;
	// resize the vector
	tmp.resize(numSolSize);
	tmp[numSolSize-1] = Tmelt;
	tmp[0] = Ttop;
	// linear interpolation to new grid
	for (int i=1; i<numSolSize-1; i++) {
		double xtilde = i*dx;
		if (xtilde/dxOld < TOL) tmp[i] = solution[i];
		else {
			int ix0 = static_cast<int>(xtilde/dxOld);
			if (ix0 >= (numSolSizeOld - 1)) tmp[i] = 273.15;
			else {
				int ix1 = ix0 + 1;
				if (ix1 >= (numSolSizeOld-1)) ix1 = numSolSizeOld - 1;
				// Ttilde = (T0*(x1 - xtilde) + T1*(xtilde - x0)) / (x1 - x0)
				tmp[i] = (solution[ix0]*(ix1*dxOld - xtilde) + solution[ix1]*(xtilde - ix0*dxOld)) / ((ix1-ix0)*dxOld);
			}
		}
	}
	solution.resize(numSolSize);
	for (int i=0; i<numSolSize; i++) solution[i] = tmp[i];

	if (numSolSizeOld != numSolSize) return true;
	return false;
}

void solidify(Matrix& ymNumeric, Matrix& ymAnal, size_t& ymSize, double lambdaAnal, double time, double dt) {
	double ym = ymNumeric[ymSize-1][1];
	double lambda = ym / 2 / sqrt(kappa * time);

	ymNumeric.push_back(Vector(2));
	ymAnal.push_back(Vector(2));
	ymNumeric[ymSize][0] = time;
	ymNumeric[ymSize][1] = ym + lambda * sqrt(kappa/time)*dt;
	ymAnal[ymSize][0] = time;
	ymAnal[ymSize][1] = ymAnal[ymSize-1][1]+lambdaAnal*sqrt(kappa/time)*dt;
	ymSize++;
}

void StefanProblem1D(Vector& numSol, Vector& analSol, Vector Tinit, size_t& numSolSize, size_t& analSolSize, Matrix& ymNumeric, Matrix& ymAnal, size_t& ymSize, double maxTime, double& time, double& dx) {
	double dt = dx*dx/2/kappa;
	double maxdx = 0.001; // TODO
	double TOL = 1e-14;
	double lambda = ymNumeric[ymSize-1][1] / 2 / sqrt(kappa * time);
	double lambdaAnal = findLambda(Tmelt-Ttop, L, c, 0.5, 1e-15);

	Vector oldSol(numSolSize);
	oldSol = Tinit;
	numSol = Tinit;

	// Numerical solution + solidification interface depth
	time += dt;
	int j = static_cast<int>(ymSize);
	while (time < maxTime) {
		// update old values
		oldSol = numSol;

		// update timestep
		// if (dt > dx*dx/2/kappa) dt = dx*dx/2/kappa;
		dt = dx*dx/2/kappa;

		// update solidification interface depth
		solidify(ymNumeric, ymAnal, ymSize, lambdaAnal, time, dt);

		/* update grid :( */
		if (updateGrid(oldSol,numSolSize, ymNumeric[ymSize-1][1], dx, maxdx)) numSol.resize(numSolSize);

		// solve HEq for current time
		//dumpTempVector("tmp/dataset_"+to_string(j)+".dat", numSol, numSolSize, dx);
		centralDifferences(numSol, oldSol, numSolSize, kappa*dt/dx/dx);

		// update time
		time += dt;
		j++;
	}

	// Analytical solution
	analSol.resize(numSolSize);
	analSolSize = numSolSize;

	cout << lambda << endl;
	cout << ymAnal[ymSize-1][1] << "," << dx << endl;
	StefanAnal1D(analSol, analSolSize, dx, ymAnal[ymSize-1][1], time, lambdaAnal);
}

int main() {
	size_t N = 3;
	size_t M = 0;
	double dx = 0.0001;

	double lambda = findLambda(Tmelt-Ttop, L, c, 0.5, 1e-15);
	double time = N*dx*N*dx/4/lambda/lambda/(kappa);;

	Vector solution(N,0.0), solAn(M);
	// Initial temperature field
	Vector Tinit(N, 273.15);
	// /// Linear initial condition
	// for (int i=0; i<N-1; i++) Tinit[i] = Ttop + i*(Tmelt-Ttop)/(N-1);
	/// erf initial temp
	StefanAnal1D(Tinit, N, dx, N*dx, time, lambda);

	size_t ymSize = 1;
	Matrix ym(ymSize,Vector(2, 0.0));
	Matrix ymAnal(ymSize,Vector(2, 0.0));
	ym[0][0] = time;
	ym[0][1] = N*dx;
	ymAnal[0][0] = time;
	//ymAnal[0][1] = lambda*sqrt(kappa/time)*time; // toto neni hezke
	ymAnal[0][1] = N*dx; // toto take neni hezke

	dumpTempVector("initialCond.dat", Tinit, N, dx);

	StefanProblem1D(solution, solAn, Tinit, N, M, ym, ymAnal, ymSize, 1000.0, time, dx);
	printVector(solAn,M);
	printVector(solution, N);
	dumpTempVector("numSol.dat", solution, N, dx);
	dumpTempVector("analSol.dat", solAn, N, dx);
	dumpMatrix("ym.dat", ym, ymSize, 2);
	dumpMatrix("ymAnal.dat", ymAnal, ymSize, 2);
}
