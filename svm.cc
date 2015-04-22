#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <sstream>

#define SIZE 60
#define WSIZE 2
#define C 1
#define eps 0.001

using namespace std;

class point {
public:
	double x1;
	double x2;
	int y;				//real value
	double u;			//predict value
	double alpha;
};

point points[SIZE+1]; 
double w[WSIZE] = {0.0, 0.0};
double b;
double a1, a2, p;

void readfile()
{
	ifstream ifs("xaa.txt");
	string s;
	int m = 0;
	while (getline(ifs, s))
	{
		istringstream sin(s);
		double d1, d2;
		int d3;
		sin >> d1 >> d2 >> d3;
		points[m].x1 = d1;
		points[m].x2 = d2;
		points[m++].y = d3;
	}
	return;
}

bool is0C(double alpha)
{
	if(alpha == 0 || alpha == C)
		return true;
	return false;
}


bool eq(point p1, point p2)
{
	if(p1.x1 == p2.x1 && p1.x2 == p2.x2)
		return true;
	return false;
}

bool KKT(point p)
{
	double alpha = p.alpha;
	int y = p.y;
	double u = p.u;
	/*double E = u - y;
	double tol = 0.01;
	double r = y * E;
	if( (r < -tol && alpha < C) || (r > tol && alpha > 0) )
		return false;
	return true;*/
	
	if(alpha == 0 && y*u >= 1.0)
		return true;
	if(alpha == C && y*u <= 1.0)
		return true;
	if(!is0C(alpha) && fabs(y*u-1.0) < eps)
		return true;
	return false;
}

int find2nd(point p1)
{
	int n = -1;
	double E1 = p1.u - p1.y;
	double optE = -9999;
	for (int i=0;i<SIZE;i++){
		double E2 = points[i].u - points[i].y;
		double det = fabs(E2-E1);
		if(det > optE)
		{
			optE = det;
			n = i;
		}
	}
	/*if(n == -1){
		n = rand() % SIZE;
		while(eq(points[n], p1)){
			n = rand() % SIZE;
		}
	}*/
	return n;
}

double max(double d1, double d2)
{
	if(d1 > d2)
		return d1;
	return d2;
}

double min(double d1, double d2)
{
	if(d1 > d2)
		return d2;
	return d1;
}

double kernel(point p1, point p2)
{
	/*double d1 = p1.x1 - p2.x1;
	double d2 = p1.x2 - p2.x2;
	double d = (d1*d1 + d2*d2) / 0.02;
	return pow(2.71828183, -d);*/
	double inner = p1.x1*p2.x1 + p1.x2*p2.x2;
	return tanh(a1*inner+a2);
	return pow(a1*inner+a2, p);
}

double W(point p1, point p2, double alpha2)
{
	double alpha1 = p1.alpha;
	double k11 = kernel(p1, p1);
	double k22 = kernel(p2, p2);
	double k12 = kernel(p1, p2);
	int y1 = p1.y;
	int y2 = p2.y;
	int s = y1 * y2;
	
	double W = k11*alpha1*alpha1/2;
	W = k22*alpha2*alpha2/2 + W;
	W = s*k12*alpha1*alpha2 + W;
	double v1 = p1.u+b-y1*alpha1*k11-y2*alpha2*k12;
	double v2 = p2.u+b-y1*alpha1*k12-y2*alpha2*k22;
	W = y1*alpha1*v1+y2*alpha2*v2 + W;
	W = W - alpha1 - alpha2;
	return W;
}

void calc_u()
{
	for(int i=0;i<SIZE;i++){
		double d = 0.0;
		for(int j=0;j<SIZE;j++){
			point p = points[j];
			d += kernel(p, points[i])*p.alpha*p.y;
		}
		points[i].u = d - b;
	}
}

bool opt(int first, int second)
{
	if(eq(points[first], points[second]))
		return false;
	double alpha1 = points[first].alpha;
	double alpha2 = points[second].alpha;
	double old_alpha1 = alpha1;
	double old_alpha2 = alpha2;
	int y1 = points[first].y;
	int y2 = points[second].y;
	int s = y1 * y2;
	double E1 = points[first].u - points[first].y;
	double E2 = points[second].u - points[second].y;
	double L, H;
	if(y1 == y2){
		L = max(0, alpha2+alpha1-C);
		H = min(C, alpha2+alpha1);
	}
	else {
		L = max(0, alpha2-alpha1);
		H = min(C, C+alpha2-alpha1);
	}
	if(fabs(L-H) < eps)
		return false;
		
	double k11 = kernel(points[first], points[first]);
	double k22 = kernel(points[second], points[second]);
	double k12 = kernel(points[first], points[second]);
	double eta = k11 + k22 - 2*k12;
	if (eta > 0){
		alpha2 = (E1 - E2) * y2 / eta + alpha2;
		if (alpha2 > H)
			alpha2 = H;
		else if (alpha2 < L)
			alpha2 = L;
	}
	else {
		if(W(points[first], points[second], L) < W(points[first], points[second], H))
			alpha2 = L;
		else 
			alpha2 = H;
	}
	if(fabs(old_alpha2-alpha2) < eps)
		return false;
	alpha1 = (old_alpha2-alpha2)*s+alpha1;
	points[first].alpha = alpha1;
	points[second].alpha = alpha2;
	
	double b1 = b + E1;
	b1 += (alpha1-old_alpha1)*y1*k11;
	b1 += (alpha2-old_alpha2)*y2*k12;
	double b2 = b + E2;
	b2 += (alpha1-old_alpha1)*y1*k12;
	b2 += (alpha2-old_alpha2)*y2*k22;
	if(!is0C(alpha1))
		b = b1;
	else if(!is0C(alpha2))
		b = b2;
	else
		b = (b1+b2)/2;
		
	/*double w1 = (alpha1 - old_alpha1)*y1*points[first].x1;
	double w2 = (alpha2 - old_alpha2)*y2*points[second].x1;
	w[0] += w1 + w2;
	w1 = (alpha1 - old_alpha1)*y1*points[first].x2;
	w2 = (alpha2 - old_alpha2)*y2*points[second].x2;
	w[1] += w1 + w2;*/
	calc_u();
	return true;
}

int examine(int first)
{
	if(KKT(points[first]))
		return 0;
	int n = 0;
	for(int i=0;i<SIZE;i++){
		if(!is0C(points[i].alpha))
			n++;
	}
	if(n > 0){
		int second = find2nd(points[first]);
		if (opt(first, second))
			return 1;
	}
	else {
		n = rand() % SIZE;
		while(eq(points[n], points[first])){
			n = rand() % SIZE;
		}
		if (opt(first, n))
			return 1;
	}
	return 0;
}

void init()
{
	for(int i=0;i<SIZE;i++)
	{
		points[i].alpha = 0.0;
		//w[0] += points[i].alpha*points[i].y*points[i].x1;
		//w[1] += points[i].alpha*points[i].y*points[i].x2;
	}
	b = 0.0;
	calc_u();
	return;
}

void outputA()
{
	ofstream ofs("ay.txt");
	for (int i=0;i<SIZE;i++)
	{
		ofs << points[i].alpha << endl;
	}
}

int main()
{
	cout << "a1:\t";
	cin >> a1;
	cout << "a2\t",
	cin >> a2;
	cout << "power:\t";
	cin >> p;
	readfile();	
	init();
	int numChanged = 0;
	bool examineAll = true;
	while (numChanged > 1 || examineAll)
	{
		numChanged = 0;
		if(examineAll){
			for(int i=0;i<SIZE;i++)
			{
				numChanged += examine(i);
			}
		}
		else {
			for(int i=0;i<SIZE;i++)
			{
				if( !is0C(points[i].alpha) ){
					numChanged += examine(i);
				}
			}
		}
		if(examineAll)
			examineAll = false;
		else if (numChanged == 0)
			examineAll = true;
	}
	
	outputA();
	for(int i=0;i<SIZE;i++)
	{
		w[0] += points[i].alpha*points[i].y*points[i].x1;
		w[1] += points[i].alpha*points[i].y*points[i].x2;
	}
	cout << "--------------------------result--------------------------\n";
	cout << w[0] << "\t" << w[1] << "\t" << b << endl;
	return 0;
}
