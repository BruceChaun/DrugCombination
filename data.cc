#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

int main()
{
	/*double a[60];
	double x1[60];
	double x2[60];
	ifstream ifs("xaa.txt");
	string s;
	int m = 0;
	while (getline(ifs, s))
	{
		int n = s.length();
		int loc = s.find(" ");
		string s1 = s.substr(0, loc);
		string s2 = s.substr(loc+1, n);
		n  = n - loc - 1;
		loc = s2.find(" ");
		string s3 = s2.substr(0, loc);
		string s4 = s2.substr(loc+1, n);
		double d1 = atof(s1.c_str());
		double d2 = atof(s3.c_str());
		double d3 = atof(s4.c_str());
		x1[m] = d1;
		x2[m] = d2;
		a[m++] = d3;
	}
	for (int i=0;i<m;i++)
	{
		cout << x1[i] << "\t" << x2[i] << "\t" << a[i] << endl;
	}
	return 0;*/
	
	ofstream ofs("data.txt");
	for (int i=0;i<100;i++)
	{
		double x = (rand() % 101) / 100.0;
		double fluc = rand() % 101 / 100.0;
		double y = x*x;
		if(rand() %2 == 0)
			y += fluc;
		else
			y -= fluc;
		ofs << x << " " << y << " ";
		if(rand() %2 == 0)
			ofs << 1 << "\n";
		else
			ofs << -1 << "\n";
	}
	return 0;
}
