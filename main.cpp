#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <chrono>

using namespace std;

using Matrix = vector<vector<double>>;
using Vector = vector<double>;

// Define constants
const double k = 2.0; // J/(m*s*K)
const double c = 4.0e3; // J*(kg*K)
const double rho = 1.0e3; // kg/m^3
const double L = 320.0e3; // J/kg
const double kappa = k/c/rho;

const double CFL = 0.5;
const double maxdx = 0.001; // TODO

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
	for (size_t i=0; i<rows; i++) {
		for (int j=0; j<cols; j++) {
			file << matrix[i][j] << " ";
		}
		file << endl;
	}
}

void dumpTempVector(string filepath, Vector& vector, int size, double dx) {
	ofstream file(filepath);
	for (size_t i=0; i<size; i++) {
		file << i*dx << " " << vector[i] << endl;
	}
}

void centralDifferences(Vector& y, Vector& y_old, size_t N, double a) {
	// y(x_i, t_{n+1}) = y(x_{i+1}, t_n - delta t) + a * [ y(x_{i+1}, t_{n-1}) - 2*y(x_i, t_{n-1}) + y(x_{i-1}, t_{n-1}) ]
	y[0] = y_old[0];
	y[N-1] = y_old[N-1];
	for (size_t i=1; i<N-1; i++) {
		y[i] = y_old[i] + a*(y_old[i+1] - 2*y_old[i] + y_old[i-1]);
	}
}

Vector thomas(Vector& a, Vector& b, Vector& c, Vector& rhs, size_t N) {
	// Thomas algorithm for solving a system of linear equations with tridiagonal matrix
	// MODIFIES vectors b, c, rhs
	//
	// the problem:
	// ( b0 c0  0  0 0     )             ( d0 )
	// ( a1 b1 c1  0 0     )             ( d1 )
	// (  0 a2 b2 c2 0 ... ) * x = RHS = ( d2 )
	// (  .                )             ( .  )
	// (  .                )             ( .  )
	// (  .                )             ( .  )
	// 
	// a = (  0 a1 a2 ... aN )
	// b = ( b0 b1 b2 ... bN )
	// c = ( c0 c1 c2 ...  0 )
	
	// forward run: sort-of-Gauss
	// ( b0 c0  0  0 0     )    (  1 c'0   0  0   0     )
    // ( a1 b1 c1  0 0     )    (  0  1  c'1  0   0     )    a = (  0   0   0  ...  0 )
    // (  0 a2 b2 c2 0 ... ) -> (  0  0   1  c'2  0 ... ) -> b = (  1   1   1  ...  1 )
    // (  .                )    (  .          .         )    c = ( c'0 c'1 c'2 ...  0 )
    // (  .                )    (  .              .     )
    // (  .                )    (  .                 .  )

	c[0] = c[0]/b[0]; // c0' = c0/b0
	rhs[0] = rhs[0]/b[0]; // d0' = d0/b0
	for (size_t i=1; i<N; i++) {
		// ci' = ci/(bi - ai*c(i-1)')
		c[i] /= b[i]-a[i]*c[i-1];
		// di' = (di - ai*d(i-1)')/(bi - ai*c(i-1)')
		rhs[i] = (rhs[i] - a[i]*rhs[i-1]) / (b[i] - a[i]*c[i-1]);
	}

	// back substitution
	Vector x(N, 0.0);
	// !!! N refers to the size of rhs, which is Nx-2 -- it lacks boundary nodes !!!
	x[N-1] = rhs[N-1];
	for (int i=N-2; i>=0; --i) x[i] = rhs[i] - c[i]*x[i+1];
	return x;
}

void crankNicholson(Vector& y, Vector& yOld, size_t N, double mu) {
	// build rhs
	size_t M = N-2;
	Vector rhs(M, 0.0);
	for (size_t i=1; i<N-1; i++) {
		rhs[i-1] = mu*yOld[i-1] + (2-2*mu)*yOld[i] + mu*yOld[i+1];
	}
	rhs[0] += mu*Ttop;
	rhs[M-1] += mu*Tmelt;

	// build matrix
	Vector a(M,0.0), b(M, 0.0), c(M,0.0);
	for (size_t i=0; i<M; i++) {
		b[i] = 2+2*mu;
		if (i>0) a[i] = -mu;
		if (i<M-1) c[i] = -mu;
	}

	// solve system of linear equations
	Vector intSol = thomas(a, b, c, rhs, M);
	for (size_t i=0; i<M; i++) y[i+1] = intSol[i];
	y[0] = Ttop;
	y[N-1] = Tmelt;
}

double findLambda(double dT, double latentHeat, double heatCap, double initGuess, double TOL) {
	double RHS = latentHeat * sqrt(M_PI) / (heatCap * dT);
	int i = 0;

	double lambda = initGuess, lambda_old = lambda + TOL;
	// Newton
	double f,fprime;
	double ee, erflm;

	f=2*TOL;
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
	for (size_t i=0; i<size; i++) {
		// for y<ym: theta = erf(eta)/erf(lambda) -> T = T0 + (Tm-T0)*erf(y/(2*sqrt(kappa*t)))/erf(lambda)
		if (i*dx > ym) solution[i] = Tmelt;
		else {
			eta = i*dx / 2 / sqrt(kappa*time);
			solution[i] = Ttop + deltaT*erf(eta)/erfEtam;
		}
	}
}

bool updateGrid(Vector& solution, size_t& numSolSize, double ymCurrent, double& dx, double maxdx) {
	size_t numSolSizeOld = numSolSize;
	double dxOld = dx;
	double TOL = 1e-14;
	Vector tmp;

	if (ymCurrent/(numSolSize-1) > maxdx) numSolSize++;
	// update spatial step -- equidistant grid
	dx = ymCurrent/(numSolSize-1);
	// resize the temporary vector
	tmp.resize(numSolSize);
	tmp[numSolSize-1] = Tmelt;
	tmp[0] = Ttop;
	// linear interpolation to new grid
	for (size_t i=1; i<numSolSize-1; i++) {
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
	for (size_t i=0; i<numSolSize; i++) solution[i] = tmp[i];

	if (numSolSizeOld != numSolSize) return true;
	return false;
}

void solidify(Vector& numSol, size_t solSize, Matrix& ymNumeric, Matrix& ymAnal, size_t& ymSize, double lambdaAnal, double time, double dt, double dx) {
	ymNumeric.push_back(Vector(2));
	ymAnal.push_back(Vector(2));

	ymNumeric[ymSize][0] = time;
	ymNumeric[ymSize][1] = ymNumeric[ymSize-1][1] + k/rho/L * (numSol[solSize-1] - numSol[solSize-2])/dx * dt;

	ymAnal[ymSize][0] = time;
	ymAnal[ymSize][1] = ymAnal[ymSize-1][1] + lambdaAnal * sqrt(kappa/time) * dt;

	ymSize++;
}

void StefanProblem1D(Vector& numSol, Vector& analSol, Vector& Tinit, size_t& solSize, Matrix& ymNumeric, Matrix& ymAnal, size_t& ymSize, double maxTime, double& time, double& dx, unsigned int solver=0) {
	double dt = dx*dx/2/kappa;
	double TOL = 1e-14;
	double lambda = ymNumeric[ymSize-1][1] / 2 / sqrt(kappa * time);
	double lambdaAnal = findLambda(Tmelt-Ttop, L, c, 0.5, 1e-15);

	Vector oldSol(solSize);
	oldSol = Tinit;
	numSol = Tinit;

	// Numerical solution + solidification interface depth
	while (time < maxTime) {
		// update old values
		oldSol = numSol;

		// update time
		time += dt;
		if (time>maxTime) {
			dt = maxTime - time;
			time = maxTime;
		}

		// update timestep
		// if (dt > dx*dx/2/kappa) dt = dx*dx/2/kappa;
		dt = dx*dx/kappa*CFL/2;

		// update solidification interface depth
		solidify(numSol, solSize, ymNumeric, ymAnal, ymSize, lambdaAnal, time, dt, dx);

		/* update grid :( */
		if (updateGrid(oldSol,solSize, ymNumeric[ymSize-1][1], dx, maxdx)) numSol.resize(solSize); // analSol.resize(solSize);

		// solve HEq for current time
		//dumpTempVector("tmp/dataset_"+to_string(j)+".dat", numSol, solSize, dx);
		switch(solver) {
			case 1:
				crankNicholson(numSol, oldSol, solSize, kappa*dt/dx/dx);
				break;
			case 0:
			default:
				centralDifferences(numSol, oldSol, solSize, kappa*dt/dx/dx);
				break;
		}

	}

	// Analytical solution
	analSol.resize(solSize);

	StefanAnal1D(analSol, solSize, dx, ymAnal[ymSize-1][1], time, lambdaAnal);
}

int main() {
	auto start = chrono::high_resolution_clock::now();
	size_t N = 5;
	size_t M = 0;
	double dx = 0.001;
	double time;
	double lambda = findLambda(Tmelt-Ttop, L, c, 0.5, 1e-15);

	size_t ymSize = 1;
	Matrix ym(ymSize,Vector(2, 0.0));
	Matrix ymAnal(ymSize,Vector(2, 0.0));

	Vector Tinit(N, 273.15);
	Vector solution(N,0.0), solAn(M);

	/* Initial strength of solidified material */
	ym[0][1] = (N-1)*dx;
	ymAnal[0][1] = (N-1)*dx;

	time = (N-1)*dx*(N-1)*dx/4/lambda/lambda/kappa;
	ym[0][0] = time;
	ymAnal[0][0] = time;

	/* Initial temperature field */
	// /// Linear initial condition
	// for (size_t i=0; i<N-1; i++) Tinit[i] = Ttop + i*(Tmelt-Ttop)/(N-1);
	/// erf initial temp
	StefanAnal1D(Tinit, N, dx, ym[0][1], time, lambda);

	dumpTempVector("initialCond.dat", Tinit, N, dx);

	StefanProblem1D(solution, solAn, Tinit, N, ym, ymAnal, ymSize, 3*24*3600.0, time, dx, 1);
	dumpTempVector("numSol.dat", solution, N, dx);
	dumpTempVector("analSol.dat", solAn, N, dx);
	dumpMatrix("ym.dat", ym, ymSize, 2);
	dumpMatrix("ymAnal.dat", ymAnal, ymSize, 2);

	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "finished after " << duration.count()/1e3 << " s" << endl;
}
