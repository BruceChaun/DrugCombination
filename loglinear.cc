/*
 * atc only, 8 dim
418	22	154	155
precision;	0.95
recall:	0.729494
accuracy:	0.763685
F-measure:	0.825271

point 12 dim
400	21	155	173
precision;	0.950119
recall:	0.69808
accuracy:	0.740988
F-measure:	0.804829


19 dim
416	22	154	157
precision;	0.949772
recall:	0.726003
accuracy:	0.761015
F1-measure:	0.822948

22 dim
56	2	16	1
precision;	0.965517
recall:	0.982456
accuracy:	0.96
F1-measure:	0.973913


 */ 

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <sstream>

#define eps 0.01
#define DIM 22
#define C1 1e-3
#define C2 0.1
#define MAX_ITER 1e3
#define AMIN 0.05
#define AMAX 0.2
#define BOUND 5
#define SAMPLE 113	//training
#define TEST 75	//testing
#define LAMBDA 1

using namespace std;

class drug{
public:
	string dcc1;
	string atc1;
	string dcc2;
	string atc2;
	double enz;
	double tar;
	bool isCombine;
	bool predict;
};

drug drugs[SAMPLE];

double v_f(double v[], bool y, int i)
{
	double d = 0.0;
		if(drugs[i].atc1.substr(0,1) == drugs[i].atc2.substr(0,1) 
		&& drugs[i].atc1.substr(1,6) != drugs[i].atc2.substr(1,6) 
		&& y == true )
		{
			d += v[0];
		}
		if(drugs[i].atc1.substr(0,3) == drugs[i].atc2.substr(0,3) 
		&& drugs[i].atc1.substr(3,4) != drugs[i].atc2.substr(3,4) 
		&& y == true )
		{
			d += v[1];
		}
		if(drugs[i].atc1.substr(0,4) == drugs[i].atc2.substr(0,4) 
		&& drugs[i].atc1.substr(4,3) != drugs[i].atc2.substr(4,3) 
		&& y == true )
		{
			d += v[2];
		}
		if(drugs[i].atc1.substr(0,5) == drugs[i].atc2.substr(0,5) 
		&& drugs[i].atc1.substr(5,2) != drugs[i].atc2.substr(5,2) 
		&& y == true )
		{
			d += v[3];
		}
		if(drugs[i].atc1.substr(0,1) == drugs[i].atc2.substr(0,1) 
		&& drugs[i].atc1.substr(1,6) != drugs[i].atc2.substr(1,6) 
		&& y == false )
		{
			d += v[4];
		}
		if(drugs[i].atc1.substr(0,3) == drugs[i].atc2.substr(0,3) 
		&& drugs[i].atc1.substr(3,4) != drugs[i].atc2.substr(3,4) 
		&& y == false )
		{
			d += v[5];
		}
		if(drugs[i].atc1.substr(0,4) == drugs[i].atc2.substr(0,4) 
		&& drugs[i].atc1.substr(4,3) != drugs[i].atc2.substr(4,3) 
		&& y == false )
		{
			d += v[6];
		}
		if(drugs[i].atc1.substr(0,5) == drugs[i].atc2.substr(0,5) 
		&& drugs[i].atc1.substr(5,2) != drugs[i].atc2.substr(5,2) 
		&& y == false )
		{
			d += v[7];
		}
		/*if(y == true){
			d += drugs[i].enz * v[8];
			d += drugs[i].tar * v[9];
		}
		if(y == false){
			d += drugs[i].enz * v[10];
			d += drugs[i].tar * v[11];
		}*/
		double enz = drugs[i].enz;
		double tar = drugs[i].tar;
		if(enz <= 0.1 && enz > 0 && y == true){
			d += v[8];
		}
		if(enz <= 0.25 && enz > 0.1 && y == true){
			d += v[9];
		}
		if(enz > 0.25 && y == true){
			d += v[10];
		}
		if(enz == 0 && y == true){
			d += v[11];
		}
		if(tar <= 0.1 && tar > 0 && y == true){
			d += v[12];
		}
		if(tar <= 0.25 && tar > 0.1 && y == true){
			d += v[13];
		}
		if(tar > 0.25 && y == true){
			d += v[14];
		}
		if(tar == 0 && y == true){
			d += v[15];
		}
		if(enz <= 0.2 && enz > 0 && y == false){
			d += v[16];
		}
		if(enz <= 0.3 && enz > 0.2 && y == false){
			d += v[17];
		}
		if(enz > 0.3 && y == false){
			d += v[18];
		}
		if(enz == 0 && y == false){
			d += v[19];
		}
		if(tar > 0 && y == false){
			d += v[20];
		}
		if(tar == 0 && y == false){
			d += v[21];
		}
		
	return d;
}

void init(double x[])
{
	for(int i=0;i<DIM;i++)
		x[i] = 0.0;
}

double f(double x[])
{
	/*
	 * y = sum(v.f) - sum log sum(e^{v.f})
	 */
	//double y =-0.1*x[0]*x[0]-0.1*x[1]*x[1]+sin(x[0])+sin(x[1]);//0.5*x[0]+sin(x[0])-sin(2*x[0]);//x[0]* x[0] - 3*x[0];//0.5*x[0]+sin(x[0])-sin(2*x[0])+2*sin(3*x[0]);
	double L = 0.0;
	for(int i=0; i<SAMPLE; i++){
		L += v_f(x, drugs[i].isCombine, i);
		double sum = 0.0;
		sum = exp(v_f(x, true, i)) + exp(v_f(x, false, i));
		L -=  log(sum);
	}
	for(int i=0; i<DIM; i++){
		L -= x[i] * x[i] / 2  * LAMBDA;
	}
	return L;
}

double* gradient_f(double x[])
{
	static double grad[DIM];
	for(int i=0;i<DIM;i++){
		x[i] += eps;
		double diff1 = f(x);
		x[i] -= eps*2;
		double diff2 = f(x);
		x[i] += eps;
		double diff = diff1 - diff2;;
		grad[i] = diff / eps / 2;
	}
	return grad;
}

double h(double a, double x[])
{
	double grads[DIM];
	double *grad = gradient_f(x);
	for(int i=0;i<DIM;i++){
		grads[i] = *(grad+i) * a + x[i];
	}
	return f(grads);
}

double dh(double a, double x[])
{
	return ( h(a+eps, x)-h(a-eps, x) ) / eps / 2;
}

double multiVector(double *x1, double *x2)
{
	double s = 0.0;
	for(int i=0;i<DIM;i++)
	{
		s += *(x1 + i) * *(x2 + i);
	}
	return s;
}

double quadraticInterpolation(double a, double x[])
{
	return a * 2;
	double numerator = dh(0, x)*a*a;
	double denominator = h(0, x) + dh(0, x)*a - h(a, x);
	if(fabs(denominator) < eps*eps)
		return a;
	return numerator / denominator / 2;
}

double zoom(double a_low, double a_high, double x[])
{
	int k = 0;
	double h_low = h(a_low, x);
	double h_high = h(a_high, x);
	while( (k < MAX_ITER) && fabs(a_high - a_low) > eps){
		double a_mid = a_low * 1.0+ (a_high - a_low)/2.0;
		double h_mid = h(a_mid, x);
		/*double *d = gradient_f(x);
		double inner = multiVector(d, d);
		double f0 = inner * C1 * a_mid + f(x);*/
		double inner = dh(0, x);
		double f0 = inner * C1 * a_mid + f(x);
		if( (h_mid < f0) || (h_mid < h_low) ){
			a_high = a_mid;
			h_high = h_mid;
		}
		else {
			double g_mid = dh(a_mid, x);
			if( fabs(g_mid) <= fabs(C2 * inner))
				return a_mid;
			if(g_mid <= 0){
				a_high = a_mid;
				h_high = h_mid;
			}
			else {
				a_low = a_mid;
				h_low = h_mid;
			}
		}
		k ++;
	}
	return a_low;//ahigh;
}

double lineSearch(double a0, double x[])
{
	int k = 0;
	double a_pre = 0.01;
	double a_cur = a0;
	double h_pre = f(x);
	while ( (k < MAX_ITER) && fabs(a_cur - a_pre) > eps){
		double h_cur = h(a_cur, x);
		/*double *d = gradient_f(x);
		double inner = multiVector(d, d);
		double f0 = inner * C1 * a_cur + f(x);*/
		double inner = dh(0, x);
		double f0 = inner * C1 * a_cur + f(x);
		if( (h_cur < f0) || ((h_cur <= h(a_pre, x)) && k>0) )
			return zoom(a_pre, a_cur, x);
		double g_cur = dh(a_cur, x);
		if( fabs(g_cur) <= fabs(C2 * inner) )
			return a_cur;
		if(g_cur <= 0)
			return zoom(a_cur, a_pre, x);
		double a_new = quadraticInterpolation(a_cur, x);
		a_pre = a_cur;
		a_cur = a_new;
		h_pre = h_cur;
		k ++;
	}
	return AMIN;//a_cur;
}

void trainingFile()
{
	ifstream ifs("smallTraining.txt");
	string s;
	int i = 0;
	while(getline(ifs, s)){
		istringstream sin(s);
		string dcc1, dcc2, atc1, atc2;
		double enz, tar;
		bool result;
		sin >> dcc1 >> atc1 >> dcc2 >> atc2 >> enz >> tar >> result;
		drugs[i].dcc1 = dcc1;
		drugs[i].atc1 = atc1;
		drugs[i].dcc2 = dcc2;
		drugs[i].atc2 = atc2;
		drugs[i].enz = enz;
		drugs[i].tar = tar;
		drugs[i].isCombine = result;
		i++;
	}
}

void testFile(){
	ifstream ifs("smallTesting.txt");
	string s;
	int i = 0;
	while(getline(ifs, s)){
		istringstream sin(s);
		string dcc1, dcc2, atc1, atc2;
		double enz, tar;
		bool result;
		sin >> dcc1 >> atc1 >> dcc2 >> atc2 >> enz >> tar >> result;
		drugs[i].dcc1 = dcc1;
		drugs[i].atc1 = atc1;
		drugs[i].dcc2 = dcc2;
		drugs[i].atc2 = atc2;
		drugs[i].enz = enz;
		drugs[i].tar = tar;
		drugs[i].isCombine = result;
		i++;
	}
}

bool isConverge(double d[])
{
	for(int i=0; i<DIM; i++){
		if(fabs(*(d+i)) > eps)
			return false;
	}
	return true;
}

int main()
{
	//drug tDrugs[TEST];
	//testFile();
	trainingFile();
	/*for(int i=0; i<TEST; i++){
		cout << drugs[i].dcc1 << " " << drugs[i].atc1 << " " << drugs[i].dcc2 << " " << drugs[i].atc2 << " " << drugs[i].isCombine << endl;
	}
	return 0;*/
	double x[DIM];
	init(x);
	double y = f(x);
	while(true) {
		double *grad = gradient_f(x);
		if(isConverge(grad))
			break;
		double a = lineSearch(AMAX, x);
		for(int i=0;i<DIM;i++){
			x[i] = *(grad + i) * a + x[i];
			/*if(x[i] > BOUND)
				x[i] = BOUND;
			if(x[i] < -BOUND)
				x[i] = -BOUND;*/
		}
		double y_new = f(x);
		
		if(fabs(y - y_new) < eps){
			break;
		}
		if(y < y_new)
			y = y_new;
	}
	//cout << y << endl;
	
	//drug tDrugs[TEST];
	testFile();
	for(int i=0; i<TEST; i++){
		double d1 = v_f(x, true, i);
		double d2 = v_f(x, false, i);
		if(d1 > d2)
			drugs[i].predict = true;
		else
			drugs[i].predict = false;
	}
	int TP = 0;
	int FP = 0;
	int TN = 0;
	int FN = 0;
	for(int i=0; i<TEST; i++){
		if(drugs[i].isCombine == true && drugs[i].predict == true)
			TP ++;
		else if(drugs[i].isCombine == true && drugs[i].predict == false)
			FN ++;
		else if(drugs[i].isCombine == false && drugs[i].predict == true)
			FP ++;
		else 
			TN ++;
	}
	cout << TP << "\t" << FP << "\t" << TN << "\t" << FN << endl;
	cout << "precision;\t" << 1.0 * TP / (TP + FP) << endl;
	cout << "recall:\t" << 1.0 * TP / (TP + FN) << endl;
	cout << "accuracy:\t" << 1.0 * (TP + TN) / TEST << endl;
	cout << "F1-measure:\t" << 2.0 * TP / (2.0 * TP + FN + FP) << endl;
	return 0;
}
