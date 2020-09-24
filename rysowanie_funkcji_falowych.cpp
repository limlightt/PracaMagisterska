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

complex<double> gaussian(complex<double> x_n, complex<double>y_n, double x, double y);


double vf = 0.3746;						//predkosc fermiego
double m = 0.067;						//masa
double w = 10 / 27211.6;				// omega
int N = 21;
int i = 0;
double Vg = 200 / 27211.6;
double R = 40 / 0.05292;
double dx = 1.5*R / (((N - 1) / 2.));	//krok siatki
double alpha = 0.01 / dx;

double Nprim = 81;
double dx_prim = R / ((Nprim - 1) / 2), dy_prim = R / ((Nprim - 1) / 2);	

double eta = 1;

complex<double> j(0, 1);
complex<double> f_1A = 0, f_2A = 0, f_1B = 0, f_2B = 0;
double phi = 0.0;
double dphi = M_PI/50;

int main() {

	complex<double> PsiA = 0;
	complex<double> PsiB = 0;
	complex<double>* C = new complex<double>[N*N * 2];
	complex<double> mA, mB;

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

	alpha = 0.00009; 
	string filename = "d:/agh/magisterka/docelowy_rysowanie_funkcji/Projekt1/Projekt1/wyniki/N_" + to_string(N) + "/a_" + to_string(alpha) + "_m_eta1.txt";
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
	
	for (int C_param = (2*N*N/2.)-4; C_param < (2 * N*N/2)+4; C_param++) {
		cout << C_param << endl;
		for (int i = 0; i < 2 * N*N; i++) {
			C[i] = { solver.eigenvectors().col(C_param)[i].real(), solver.eigenvectors().col(C_param)[i].imag() };
		}
		
		string filenamePsi = "psiAB_N" + to_string(C_param) + ".txt";
		ofstream filePsi(filenamePsi);
		for (double x = -R; x <= R; x += dx_prim) {
			for (double y = -R; y <= R; y += dy_prim) {

				PsiA = 0.0;
				PsiB = 0.0;
				for (int n = 0; n < 2 * N*N; n++) {
					if (Rij[n][2] == 0.0) PsiA += C[n] * gaussian(Rij[n][0], Rij[n][1], x, y);
					else PsiB += C[n] * gaussian(Rij[n][0], Rij[n][1], x, y);
				}
				filePsi << x*0.05292 << "\t" << y*0.05292 << "\t" << PsiA.real() << "\t" << PsiA.imag() << "\t" << PsiB.real() << "\t" << PsiB.imag() << endl;

			}
			filePsi << endl;
		}
		
		filePsi.close();
	}
	
}

complex<double> gaussian(complex<double> x_n, complex<double>y_n, double x, double y) {
	return sqrt(alpha*2/M_PI)*exp(-alpha*(pow(x - x_n, 2) + pow(y - y_n, 2)));
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