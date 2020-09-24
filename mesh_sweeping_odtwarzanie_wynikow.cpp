#include "incl.h"

using namespace std;


double eta = 1;		 
double m = 0;

double v_f = 0.3746;								//predkosc fermiego
double h = 1;										//h kreslone
double e = 15 / 27211.6;							//ladunek elektronu
double R = 80 / 0.05292;							//R max						
int N = 1600;										//ilosc pktow
double dr = R / N;

double E_min = -400 / 27211.6;
double E_max = 600 / 27211.6;
double iter = 50000;
double dE = abs(E_max - E_min) / iter;
double E_exact = 7.044 * v_f / R;
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


	for (double it_m = -3; it_m <= 3; it_m++) {
		for (double it_eta = -1; it_eta <= 1; it_eta += 2) {
			m = it_m;
			eta = it_eta;

			string filename2 = "ostateczne/miejsca_zerowe_m_" + to_string(m) + "_eta_" + to_string(eta) + "_N_1600.txt";
			deleteFile(filename2);

			f1_temp.push_back(0);
			f1_temp.push_back(0);
			f1_temp.push_back(0);

			f2_temp.push_back(0);
			f2_temp.push_back(0);
			f2_temp.push_back(0);

			E_temp.push_back(0);
			E_temp.push_back(0);
			E_temp.push_back(0);



	for (int i = 0; i < iter; i++) {			
		double E = dE*i + E_min;

		// war. poczatkowe:
		// zigzag boundary: 1, 0
		// inf mass boundary: 0, i
		x1.first = 1;
		x1.second = 0;
		fun1_vals.push_back(x1);

		x2.first = 0;
		x2.second = 0;
		fun2_vals.push_back(x2);


		algorithm_for_wave_functions(E);		//E przy peti po energii, E_exact przy f(r)

												//przy petli po E: zapis do pliku gdy liczy dla r = dr -> ostatnia komorka:
												// energia, Re(f1), Im(f1), Re(f2), Im(f2)

		//petla zapisu do pliku gdy beda wykresy f(r), zapisz do pliku po kolei wartosci funkcji dla kazdego r
		//for (int i = 0; i < N-1; i++) {
		//promie�, Re(f1), Im(f1), Re(f2), Im(f2)	
		//	pushChartData(filename, 5, (R - i*dr)*0.05292, fun1_vals[i].first, fun1_vals[i].second, fun2_vals[i].first, fun2_vals[i].second);
		//}

		//zapisywanie f2 i E do tablicy
		f1_temp[0] = f1_temp[1];
		f1_temp[1] = fun1_vals[N - 1].first;

		f2_temp[0] = f2_temp[1];
		f2_temp[1] = fun2_vals[N - 1].second;

		E_temp[0] = E_temp[1];
		E_temp[1] = E*27211.6;

		//rozpoznanie miejsca zerowego
		
		if (abs(m + eta) < abs(m) && ((f1_temp[0] > 0 && f1_temp[1] < 0) || (f1_temp[1] > 0 && f1_temp[0] < 0))) {
			 a = (f1_temp[0] - f1_temp[1]) / (E_temp[0] - E_temp[1]);
			 b = f1_temp[0] - (f1_temp[0] - f1_temp[1])*E_temp[0] / (E_temp[0] - E_temp[1]);
			pushChartData(filename2, 5, (-b / a)*R / (v_f*27211.6), 0.001, 0.001, 0, (f2_temp[0] + f2_temp[1]) / 2.);
		}

		if (abs(m + eta) > abs(m) && ((f2_temp[0] > 0 && f2_temp[1] < 0) || (f2_temp[1] > 0 && f2_temp[0] < 0))) {
			 a = (f2_temp[0] - f2_temp[1]) / (E_temp[0] - E_temp[1]);
			 b = f2_temp[0] - (f2_temp[0] - f2_temp[1])*E_temp[0] / (E_temp[0] - E_temp[1]);
			pushChartData(filename2, 5, (-b / a)*R / (v_f*27211.6), (f1_temp[0] + f1_temp[1]) / 2., 0, 0.001, 0.001);
		}

		fun1_vals.clear();
		fun2_vals.clear();

	}
	f1_temp.clear();
	f2_temp.clear();
	E_temp.clear();

	cout << m << " " << eta << " done" << endl;
		}
	}


	int j;
	cin >> j;

}


void algorithm_for_wave_functions(double E) {
	for (int n = 1; n <= N; n++) {
		fun1_vals.push_back(fun1(R - (n)*dr, dr, fun2_vals[n - 1], fun1_vals[n - 1], E));		//r(liczone od srodka), fun2(r+dr)
		fun2_vals.push_back(fun2(R - (n)*dr, dr, fun1_vals[n - 1], fun2_vals[n - 1], E));   	
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
	return 0;
}

double Ub(double r) {
	return 0;
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

