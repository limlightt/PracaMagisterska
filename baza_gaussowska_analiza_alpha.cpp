#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <complex>
#include <cmath>
#include <stdarg.h>
#include <fstream>
#include <string>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Eigen/Eigenvalues"
#include "Eigen/Dense"

using namespace std;


complex<double> Sij(complex<double> ix, complex<double> iy, complex<double> jx, complex<double> jy);
complex<double> Piij(complex<double> ix, complex<double> iy, complex<double> jx, complex<double> jy);
complex<double> Piji(complex<double> ix, complex<double> iy, complex<double> jx, complex<double> jy);
complex<double> Va(complex<double> ix, complex<double> iy, complex<double> jx, complex<double> jy);
complex<double> Vb(complex<double> ix, complex<double> iy, complex<double> jx, complex<double> jy);


double vf = 0.3746;				//predkosc fermiego
double m = 0.067;				//masa
double w = 10 / 27211.6;		// omega
int N;
int i = 0;
double Vg = 200 / 27211.6;
double R = 40 / 0.05292;
double alpha = 0;

int eta = -1;

int main() {	
	double alpha_min;
	double alpha_max;
	double alpha_delta;
	int get_N;
	cin >> alpha_min;
	cin >> alpha_max;
	cin >> alpha_delta;
	cin >> get_N;
	N = get_N;
	double dx = R / (((N - 1) / 2.));		
	double alpha_base = 0.01 / dx;
	alpha_min *= alpha_base;
	alpha_max *= alpha_base;
	alpha_delta *= alpha_base;

	cout << "all read: from " << alpha_min << " to " << alpha_max << " where dx = " << dx << " and N = " << N << endl;

	complex<double>** Rij = new complex<double>*[N*N * 2];
	for (int i = 0; i < N*N * 2; ++i)
		Rij[i] = new complex<double>[3];

	int n = 0;
	for (int k = 0; k < 2; k++) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				Rij[n][0] = (i - (N - 1) / 2.) * dx;
				Rij[n][1] = (j - (N - 1) / 2.) * dx;
				Rij[n][2] = k;
				n++;
			}
		}
	}

	for (double a = alpha_min; a < alpha_max; a += alpha_delta){		
		cout << a << "\t" << to_string(a) << endl << endl;
		alpha = a;

		string filename = "d:/agh/magisterka/docelowy_modyfikowany/Projekt1/Projekt1/wyniki/alpha/alpha/eta" + to_string(eta) + "/N_" + to_string(N) + "/energie_alpha" + to_string(alpha) + ".txt";
		ofstream file(filename);

		Eigen::MatrixXcf H(2 * N * N, 2 * N * N),
							S(2 * N * N, 2 * N * N);
		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcf> solver;

		for (int n1 = 0; n1 < N*N * 2; n1++) {
			for (int n2 = 0; n2 < N*N * 2; n2++) {
				if (Rij[n1][2] == Rij[n2][2]) {
					S(n1, n2) = Sij(Rij[n1][0], Rij[n1][1], Rij[n2][0], Rij[n2][1]);
				}
				else
				{
					S(n1, n2) = 0;
				}
				//macierz H
				if (Rij[n1][2] == 0.0 && Rij[n2][2] == 0.0)//obie sa podsiec a
					H(n1, n2) = Va(Rij[n1][0], Rij[n1][1], Rij[n2][0], Rij[n2][1]);
				else if (Rij[n1][2] == 1.0 && Rij[n2][2] == 1.0)//obie sa podsiec b
					H(n1, n2) = Vb(Rij[n1][0], Rij[n1][1], Rij[n2][0], Rij[n2][1]);
				else if (Rij[n1][2] == 0.0 && Rij[n2][2] == 1.0)//podsiec a, podsiec b
					H(n1, n2) = Piij(Rij[n1][0], Rij[n1][1], Rij[n2][0], Rij[n2][1]);
				else //podsiec b, podsiec a
					H(n1, n2) = Piji(Rij[n1][0], Rij[n1][1], Rij[n2][0], Rij[n2][1]);
			}
		}

		cout << "matrixes filled" << endl;
		solver.compute(H, S);
		file << "\t" <<  solver.eigenvalues()*27211.6 << endl;
		file.close();

		for (int i = 0; i < 2 * N*N; i++) {

			string filename2 = "d:/agh/magisterka/docelowy_modyfikowany/Projekt1/Projekt1/wyniki/alpha/alpha/eta" + to_string(eta) + "/N_" + to_string(N) + "/" + to_string(alpha) + "_"+ to_string(i) + ".txt";
			ofstream file2(filename2);
			for (int j = 0; j < 2 * N*N; j++) {
				file2 << solver.eigenvectors().col(i)[j].real() << "\t" << solver.eigenvectors().col(i)[j].imag() << "\t";
			}
			file2 << endl;
		}
		
	}
	
	cout << "done" << endl;
	int g;
	cin >> g;
}

complex<double> Sij(complex<double> xi, complex<double> yi, complex<double> xj, complex<double> yj) {
	return exp(-alpha*(xi*xi + yi*yi + xj*xj + yj*yj - 2.0*xi*xj - 2.0*yi*yj) / 2.);
}

complex<double> Piij(complex<double> xi, complex<double> yi, complex<double> xj, complex<double> yj) {
	complex<double> i(0, 1);
	return (i * alpha * vf * (xi + xj) - 2.0 * i * alpha * vf * xj + vf * eta * alpha * (yi + yj) - 2.0 * alpha * vf * eta * yj)
		   *exp(-alpha * (xi*xi + yi*yi + xj*xj + yj*yj - 2.0 * xi*xj - 2.0 * yi*yj) / 2.0);
}

complex<double> Piji(complex<double> xi, complex<double> yi, complex<double> xj, complex<double> yj) {
	complex<double> i(0, 1);
	return (i * alpha * vf * (xi + xj) - 2.0 * i * alpha * vf * xj - vf * eta * alpha * (yi + yj) + 2.0 * alpha * vf * eta * yj)
		*exp(-alpha * (xi*xi + yi*yi + xj*xj + yj*yj - 2.0 * xi*xj - 2.0 * yi*yj) / 2.0);


}
complex<double> Va(complex<double> xi, complex<double> yi, complex<double> xj, complex<double> yj) {

	return (exp(-alpha*(xi*xi + yi*yi + xj*xj + yj*yj - 2.0*xi*xj - 2.0*yi*yj) / 2.0)*Vg
		- 2.0*exp(-alpha*(alpha*xi*xi*R*R + alpha*yi*yi*R*R + alpha*xj*xj*R*R + alpha*yj*yj*R*R - 2.0*R*R*alpha*xi*xj
			- 2.0*R*R*alpha*yi*yj + xi*xi + yi*yi + xj*xj + yj*yj) / (2.0*alpha*R*R + 1.0))*alpha*Vg*R*R / (2.0*alpha*R*R + 1.0));
}

complex<double> Vb(complex<double> xi, complex<double> yi, complex<double> xj, complex<double> yj) {
	return -(exp(-alpha*(xi*xi + yi*yi + xj*xj + yj*yj - 2.0*xi*xj - 2.0*yi*yj) / 2.0)*Vg
		- 2.0*exp(-alpha*(alpha*xi*xi*R*R + alpha*yi*yi*R*R + alpha*xj*xj*R*R + alpha*yj*yj*R*R - 2.0*R*R*alpha*xi*xj
			- 2.0*R*R*alpha*yi*yj + xi*xi + yi*yi + xj*xj + yj*yj) / (2.0*alpha*R*R + 1.0))*alpha*Vg*R*R / (2.0*alpha*R*R + 1.0));
}