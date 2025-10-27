#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <thread>
#include <functional>
#include <future>



using namespace std;

using Matrix = vector<vector<double>>;
using Vector = vector<double>;

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

void dumpMatrix(string filepath, Matrix& solution, int rows, int cols, double dt, double dx, bool bin=false){
	// output binary file (sort of*) readable by gnuplot
	// * skip first ?? bytes
	ofstream file;
	if (bin){
		double t, x, T;
		file.open(filepath, ios::binary);

		// no. of records
		file.write(reinterpret_cast<const char*>(&rows), sizeof(int));
		// no. of points per record
		file.write(reinterpret_cast<const char*>(&cols), sizeof(int));

		for (int j=0; j<rows; j++){
			//file.write(reinterpret_cast<const char*>(&cols), sizeof(static_cast<int32_t>(cols)));
			t = j*dt;
			for (int i=0; i<cols; i++) {
				x = i*dx;
				T = solution[j][i];
				file.write(reinterpret_cast<const char*>(&x), sizeof(double));
				file.write(reinterpret_cast<const char*>(&T), sizeof(double));
				file.write(reinterpret_cast<const char*>(&t), sizeof(double));
			}
		}
	}
	else{
		file.open(filepath);
		file << "# " << rows << " x " << cols << "; dt=" << dt << endl;
		for (int j=0; j<rows; j++) {
			for (int i=0; i<cols; i++) {
				file << i*dx << " " << solution[j][i] << " " << j*dt << endl;
			}
			file << endl << endl;
		}
	}
	file.close();
}

void dumpVector(string filepath, Vector& vector, int size, double dt, double dx) {
	ofstream file(filepath);
	for (int i=0; i<size; i++) {
		file << vector[i] << " " << 273.15 << " "  << i*dt << endl;
	}
}

void Heat1D(Matrix& solution, double dx, int noCells, double noSteps, double timestep, double diffCoef, double tempTop, double tempBottom, Vector tempInit) {
	// assumes solution and tempInit have the right dimensions (noCells x noSteps and noCells respectively)!
	double c = diffCoef * timestep / (dx * dx);
	solution[0][noCells-1] = tempBottom;
	solution[0][0] = tempTop;
	solution[0] = tempInit;

	// Central differences
	for (int j=1; j<noSteps; j++) {
		solution[j][0] = tempTop;
		solution[j][noCells-1] = tempBottom;
		for (int i=1; i<noCells-1; i++) {
			solution[j][i] = solution[j-1][i] + diffCoef*(solution[j-1][i+1] - 2*solution[j-1][i] + solution[j-1][i-1])*timestep/dx/dx;
		}
	}
}

double findLambda(double dT, double latentHeat, double heatCap, double initGuess) {
	const double TOL = 1.0e-15;
	double RHS = latentHeat * sqrt(M_PI) / (heatCap * dT);
	int i = 0;

	double lambda = initGuess, lambda_old = lambda + TOL;
	/*// use fixed-point iteration
	while (abs(lambda-lambda_old)>TOL) {
		lambda_old = lambda;
		lambda = exp(-lambda_old*lambda_old) / erf(lambda_old) / RHS;
		if (i % 10000000 ==0) cout << i << " " << lambda << " " << lambda_old << endl;
		i++;
	}*/

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

void StefanAnal(Matrix& solution, Vector& ym, int N, int T, double dx, double dt, double Ttop, double deltaT, double thermCond, double heatCap, double density, double latentHeat) {
	double lambda = findLambda(deltaT, latentHeat, heatCap, 0.3);
	double erflm = erf(lambda);
	double thermDiff = thermCond/density/heatCap;
	double eta;
	ym[0] = 0.0;

	for (int i=0; i<N; i++) solution[0][i] = Ttop + deltaT;
	for (int j=1; j<T; j++){
		ym[j] = ym[j-1] + lambda * sqrt(thermDiff/j*dt)*dt;
		for (int i=0; i<N; i++) {
			if (i*dx <= ym[j]){
				eta = i*dx/2/sqrt(thermDiff * j*dt);
				solution[j][i] = Ttop + deltaT*erf(eta)/erflm;
				//if (solution[j][i] > (Ttop + deltaT)) solution[j][i] = Ttop + deltaT;
			}
			else solution[j][i] = Ttop + deltaT;
		}
	}
}

void solidifBound(Vector& ym, int T, double dt, double deltaT, double thermCond, double heatCap, double density, double latentHeat){
	double lambda = findLambda(deltaT, latentHeat, heatCap, 0.3);
	double thermDiff = thermCond/density/heatCap;

	for (int j=1; j<T; j++) {
		ym[j] = ym[j-1] + lambda*sqrt(thermDiff/j*dt)*dt;
	}
}

int main(){
	/* Define stuff */
	const double k = 2.0; // J/(m*s*K)
	const double c = 4.0e3; // J*(kg*K)
	const double rho = 1.0e3; // kg/m^3
	const double L = 320.0e3; // J/kg

	const double Ttop = 253.15; // K; = -20 C
	const double Tbot = 273.15; // K; =  0C

	const double dt = 0.1; // s; temporal discretisation
	const double height = 0.05; // m
	const int T = 100;//8640000 s = 10 days; max time, non-dimensional
	const int N = 100; // spatial discretisation degree

	double dx = height / (N-1);

	Matrix solution(T,Vector(N,0.0));
	Matrix solutionAnal(T,Vector(N,0.0));
	Vector Tini(N,Tbot);
	Vector ym(T, 0.0);

	/* Create linear initial temperature profile */
	//double dT = 20.0/(N-1);
	//for (int i=0;i<N;i++) Tini[i] = 253.15 + i*dT;

	Tini[0] = Ttop;
	//printVector(Tini,N);

	// create threads to solve stuff
	auto thrHEq1D = async(launch::async, Heat1D, ref(solution),dx,N,T,dt,k/(rho*c), Ttop, Tbot, Tini);
	auto thrStefAn = async(launch::async, StefanAnal,ref(solutionAnal), ref(ym), N, T, dx, dt, Ttop, Tbot-Ttop, k, c, rho, L);
	// once computed, start saving to file
	auto thrHEq1DWrite = async(launch::async, [&]() {
			thrHEq1D.get();
			dumpMatrix("data.dat", ref(solution), T, N, dt, dx, true);
			});
	auto thrStefAnWrite = async(launch::async, [&]() {
			thrStefAn.get();
			dumpMatrix("anal.dat", ref(solutionAnal), T, N, dt, dx, true);
			dumpVector("ym.dat", ym, T, dt, dx);
			});
}
