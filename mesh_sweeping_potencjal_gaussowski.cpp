#pragma once

#include <iostream>
#include <vector>
#include <utility>		
#include <stdarg.h>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;


int eta;		
int m;																		

double v_f = 0.3746;								//predkosc fermiego
double v_g = 200 / 27211.6;							//napiecie na bramce

double h = 1;										//h kreslone
double e = 15 / 27211.6;							//ladunek elektronu
double R = 40 / 0.05292;							//R max								
int N = 400;										//ilosc pktow
double dr = 4 * R / N;

double E_min = -400 / 27211.6;
double E_max = 600 / 27211.6;
double iter = 5000;
double dE = abs(E_max - E_min) / iter;
double E_exact = 1.70548 * v_f / R;	
pair<double, double> x1;
pair<double, double> x2;


double a, b;

vector<pair<double, double>> fun1_vals;
vector<pair<double, double>> fun2_vals;


double Ua(double r);
double Ub(double r);
pair<double, double> fun1(double r, double dr, pair<double, double> old_fun2, pair<double, double> old_fun1, double E);
pair<double, double> fun2(double r, double dr, pair<double, double> old_fun1, pair<double, double> old_fun2, double E);
double countModule(pair<double, double> a);

void algorithm_for_wave_functions(double E);


void deleteFile(string filename);
void pushChartData(string filename, int n, ...);


vector<pair<double, double>> f1_course;
vector<pair<double, double>> f2_course;

vector<double> f1_solution;
vector<double> f2_solution;

vector<double> f2_temp;
vector<double> f1_temp;

vector<double> E_temp;


int main() {

	string filenameeta1 = "wyniki/minima_eta1.txt";
	ofstream fileout1(filenameeta1);

	string filenameetaminus1 = "wyniki/minima_eta-1.txt";
	ofstream fileout2(filenameetaminus1);


	for (double it_m = -3; it_m <=3; it_m++){
		for (double it_eta = -1; it_eta <= 1; it_eta += 2) {
			if (eta == 1) cout << "JEDEN" << endl;
			m = it_m;
			eta = it_eta;
	
			f1_temp.push_back(0);
			f1_temp.push_back(0);
			f1_temp.push_back(0);

			f2_temp.push_back(0);
			f2_temp.push_back(0);
			f2_temp.push_back(0);

			E_temp.push_back(0);
			E_temp.push_back(0);
			E_temp.push_back(0);
			
	
			for (int i = 0; i < iter; i++) {			//petla po energii
				double E = dE*i + E_min;

				// war. poczatkowe:
				// zigzag boundary: 1, 0
				// inf mass boundary: 0, i
				x1.first = 0;
				x1.second = 0;
				fun1_vals.push_back(x1);

				x2.first = 0;
				x2.second = 1;
				fun2_vals.push_back(x2);


				algorithm_for_wave_functions(E);


				f1_temp[0] = f1_temp[1];
				f1_temp[1] = fun1_vals[N - 1].first;

				f2_temp[0] = f2_temp[1];
				f2_temp[1] = fun2_vals[N - 1].second;

				E_temp[0] = E_temp[1];
				E_temp[1] = E;

				//rozpoznanie miejsca zerowego
	
				if (abs(m + eta) < abs(m) && ( (f1_temp[0] > 0 && f1_temp[1] < 0) || (f1_temp[1] > 0 && f1_temp[0] < 0) ) ) {
					a = (f1_temp[0] - f1_temp[1]) / (E_temp[0] - E_temp[1]);
					b = f1_temp[0] - (f1_temp[0] - f1_temp[1])*E_temp[0] / (E_temp[0] - E_temp[1]);
					if (eta == 1) { fileout1 << m << "\t" << eta << "\t" << (-b / a)*27211.6 << endl; }
					else { fileout2 << m << "\t" << eta << "\t" << (-b / a)*27211.6 << endl; }

				}

				if (abs(m + eta) > abs(m) && ( (f2_temp[0] > 0 && f2_temp[1] < 0) || (f2_temp[1] > 0 && f2_temp[0] < 0) ) ) {
					a = (f2_temp[0] - f2_temp[1]) / (E_temp[0] - E_temp[1]);
					b = f2_temp[0] - (f2_temp[0] - f2_temp[1])*E_temp[0] / (E_temp[0] - E_temp[1]);
					if (eta == 1) { fileout1 << m << "\t" << eta << "\t" << (-b / a)*27211.6 << endl; }
					else { fileout2 << m << "\t" << eta << "\t" << (-b / a)*27211.6 << endl;}
				}


				// kod do danych do wykresu odtwarzania wynikow
				////if (abs(m + eta) < abs(m) && ((f1_temp[0] > 0 && f1_temp[1] < 0) || (f1_temp[1] > 0 && f1_temp[0] < 0))) {
				//	//f1 - y
				//	//E - x
				//if ((f1_temp[0] > 0 && f1_temp[1] < 0) || (f1_temp[1] > 0 && f1_temp[0] < 0)){
				//	a = (f1_temp[0] - f1_temp[1]) / (E_temp[0] - E_temp[1]);
				//	b = f1_temp[0] - (f1_temp[0] - f1_temp[1])*E_temp[0] / (E_temp[0] - E_temp[1]);
				//	pushChartData(filename, 5, (-b / a)*27211.6, 0.001, 0.001, 0., (f2_temp[0] + f2_temp[1]) / 2.);
				//	cout << (f2_temp[0] + f2_temp[1]) / 2. << endl;
				//}

				////if (abs(m + eta) > abs(m) && ((f2_temp[0] > 0 && f2_temp[1] < 0) || (f2_temp[1] > 0 && f2_temp[0] < 0))) {
				//if ((f2_temp[0] > 0 && f2_temp[1] < 0) || (f2_temp[1] > 0 && f2_temp[0] < 0)){
				//	a = (f2_temp[0] - f2_temp[1]) / (E_temp[0] - E_temp[1]);
				//	b = f2_temp[0] - (f2_temp[0] - f2_temp[1])*E_temp[0] / (E_temp[0] - E_temp[1]);
				//	pushChartData(filename, 5, (-b / a)*27211.6, (f1_temp[0] + f1_temp[1]) / 2., 0., 0.001, 0.001);
				//}
				//pushChartData(filename, 5, E*27211.6, fun1_vals[N - 1].first, fun1_vals[N - 1].second, fun2_vals[N - 1].first, fun2_vals[N - 1].second);


				fun1_vals.clear();
				fun2_vals.clear();

			}
				
			f1_temp.clear();
			f2_temp.clear();
			E_temp.clear();

			cout << m << " " << eta << " done" << endl;
		}
	}
	cout << "end";
	int j;
	cin >> j;
}

void algorithm_for_wave_functions(double E) {
	for (int n = 1; n <= N; n++) {
		fun1_vals.push_back(fun1(4 * R - (n)*dr, dr, fun2_vals[n - 1], fun1_vals[n - 1], E));		
		fun2_vals.push_back(fun2(4 * R - (n)*dr, dr, fun1_vals[n - 1], fun2_vals[n - 1], E));   	
	}
}


pair<double, double> fun1(double r, double dr, pair<double, double> old_fun2, pair<double, double> old_fun1, double E) {

	double Re, Im;

	Re = -old_fun2.second * dr * Ub(r + dr) / (h * v_f) + old_fun2.second * dr * E / (h * v_f) + old_fun1.first - old_fun1.first * eta *m *dr / (r + dr);
	Im = old_fun2.first * dr * Ub(r + dr) / (h * v_f) - old_fun2.first * dr * E / (h * v_f) + old_fun1.second - old_fun1.second *eta * m * dr / (r + dr);
	return pair<double, double>(Re, Im);
}

pair<double, double> fun2(double r, double dr, pair<double, double> old_fun1, pair<double, double> old_fun2, double E) {

	double Re, Im;

	Re = -old_fun1.second * dr * Ua(r + dr) / (h * v_f) + old_fun1.second * dr * E / (h * v_f) + old_fun2.first + old_fun2.first * eta * m * dr / (r + dr) + old_fun2.first * eta * eta * dr / (r + dr);
	Im = old_fun1.first * dr * Ua(r + dr) / (h * v_f) - old_fun1.first * dr * E / (h * v_f) + old_fun2.second + old_fun2.second * eta * dr * m / (r + dr) + old_fun2.second * eta *eta *dr / (r + dr);

	return pair<double, double>(Re, Im);

}

double Ua(double r) {
	return v_g*(1.0 - exp(-r*r / (R*R)));
}

double Ub(double r) {
	return -v_g*(1.0 - exp(-r*r / (R*R)));
}


void deleteFile(string filename) {
	const char* c = filename.c_str();
	if (remove(c) != 0) cout << ("Error deleting file, could not have been created or smth \n");
	else
		cout << ("File " + filename + " successfully deleted \n");
}

void pushChartData(string filename, int n, ...)
{
	va_list values;
	va_start(values, n);

	fstream file(filename, ios::out | ios::app);
	double val;
	if (file.good())
	{
		for (int i = 0; i < n; i++)
		{
			val = va_arg(values, double);
			file << val << "\t";
		}
		file << "\n";
	}
	else cout << "problem z plikiem \n";
	va_end(values);
	file.close();

}

double countModule(pair<double, double> a) {
	return sqrt(pow(a.first, 2) + pow(a.second, 2));
}