/*
 * kernel	a1	a2	p	precision	recall	accuracy	F
 * sigmoid	1	0	/	0.652		0.682	0.625		0.667
 * poly		0.5	0	3	0.553		0.955	0.55
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std;

#define SIZE 40			//test size
#define SSIZE 60		//training size
#define w1 -1.5138
#define w2 1.31277
double b, a1, a2, p;

int main()
{
	cout << "b:\t";
	cin >> b;
	cout << "a1:\t";
	cin >> a1;
	cout << "a2\t",
	cin >> a2;
	cout << "power:\t";
	cin >> p;
	double x1[SSIZE];
	double x2[SSIZE];
	int y[SSIZE];
	double alpha[SSIZE];
	ifstream ifs("xaa.txt");
	string s;
	int m = 0;
	while (getline(ifs, s))
	{
		istringstream sin(s);
		double d1, d2;
		int d3;
		sin >> d1 >> d2 >> d3;
		x1[m] = d1;
		x2[m] = d2;
		y[m++] = d3;
		
	}
	
	ifstream ifs2("ay.txt");
	m = 0;
	while (getline(ifs2, s))
	{
		istringstream sin(s);
		double d;
		sin >> d;
		alpha[m++] = d;
		
	}
	
	double tx1[SIZE];
	double tx2[SIZE];
	int ty[SIZE];
	ifstream ifs3("xab.txt");
	m = 0;
	while(getline(ifs3, s))
	{
		istringstream sin(s);
		double d1, d2;
		int d3;
		sin >> d1 >> d2 >> d3;
		tx1[m] = d1;
		tx2[m] = d2;
		ty[m++] = d3;
	}
	
	int TP = 0;
	int FP = 0;
	int TN = 0;
	int FN = 0;
	for(int i=0;i<SIZE;i++)
	{
		double u = -b;
		//u += w1*x1[i] + w2*x2[i];
		for(int j=0;j<SSIZE;j++){
			double kernel = tx1[i]*x1[j]+tx2[i]*x2[j];
			//kernel = pow(a1*kernel + a2, p);
			kernel = tanh(a1*kernel+a2);
			u += kernel*y[j]*alpha[j];
		}
		//u = pow(0.5*u+2, 5.0) + b;
		if(u > 0 && ty[i] > 0)
			TP++;
		else if (u < 0 && ty[i] > 0)
			FN++;
		else if (u < 0 && ty[i] < 0)
			TN++;
		else 
			FP++;
	}
	cout << "precision;\t" << 1.0 * TP / (TP + FP) << endl;
	cout << "recall:\t" << 1.0 * TP / (TP + FN) << endl;
	cout << "accuracy:\t" << 1.0 * (TP + TN) / SIZE << endl;
	cout << "F-measure:\t" << 2.0 * TP / (2.0 * TP + FN + FP) << endl;
	return 0;
}
