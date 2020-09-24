#include <iostream>
#include <fstream>
#include <vector>
#include <utility>	
#include <stdarg.h>
#include <fstream>
#include <string>
#include <math.h>

#include "Eigen/Eigenvalues"
#include "Eigen/Dense"

using namespace std;

double Hnm(double x, double y);
double Snm(double x, double y);


double dx = 8 / 0.05292;				//krok siatki
double m = 0.067;						//masa
double w = 10/27211.6;					// omega
int n = 5;
double alpha = (m*w / 2.);
int i = 0;

vector<double> E;		

int main() {

	string filename = "energies_dx_min.txt";
	ofstream file(filename);

	Eigen::MatrixXd H(n, n),
		S(n, n),
		c,
		eigenValues(n, 1);
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver;


	for (dx = 0; dx <= 0.1 / 0.05292; dx += 0.001 / 0.05292) {		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				S(i, j) = Snm((i - n / 2)*dx, (j - n / 2)*dx); 
				H(i, j) = Hnm((i - n / 2)*dx, (j - n / 2)*dx);
			}
		}

		solver.compute(H, S);
		file << dx << "  ";
		file << solver.eigenvalues().transpose()*27211.6;
		file << endl;

	}
	file.close();

	cout << "end";
	int g;
	cin >> g;


}

//zwraca wartosc hamiltonianu
double Hnm(double xn, double xm) {
	return exp(-alpha*(xm*xm + xn*xn - 2.0*xm*xn) / 2) / m*sqrt(alpha)*sqrt(2.0)*sqrt(0.3141592653589793E1) / 4 
		 - exp(-alpha*(xm*xm + xn*xn - 2.0*xm*xn) / 2)*sqrt(alpha*alpha*alpha)
		 * pow(xm + xn, 2.0)*sqrt(2.0)*sqrt(0.3141592653589793E1) / m / 4 + exp(-alpha*(xm*xm + xn*xn - 2.0*xm*xn) / 2)
		 * sqrt(alpha*alpha*alpha)*(xm + xn)*xn*sqrt(2.0)*sqrt(0.3141592653589793E1) / m 
		 - exp(-alpha*(xm*xm + xn*xn - 2.0*xm*xn) / 2) / m*sqrt(alpha*alpha*alpha)*xn*xn*sqrt(2.0)*sqrt(0.3141592653589793E1)
		 + exp(-alpha*(xm*xm + xn*xn - 2.0*xm*xn) / 2)*m*w*w / sqrt(alpha)*pow(xm + xn, 2.0)*sqrt(2.0)
		 * sqrt(0.3141592653589793E1) / 16 + exp(-alpha*(xm*xm + xn*xn - 2.0*xm*xn) / 2)*m*w*w*sqrt(2.0) 
		 / sqrt(alpha*alpha*alpha)*sqrt(0.3141592653589793E1) / 16;
}

double Snm(double xn, double xm) {
	return  exp(-alpha*(xm*xm + xn*xn - 2 * xm*xn) / 2)*sqrt(2.E0) / sqrt(alpha)*sqrt(0.3141593E1) / 2.;

}


