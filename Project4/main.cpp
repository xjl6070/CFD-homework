#include <iostream>
#include <cmath>
using namespace std;
double delta_t(double Q[][3], double dx, int N,double cfl);
void init(double Q[][3], double gamma, int N, int value);
void bound(double Q[][3], int N);
void TVD(double Q[][3], double F[][3], int N, double cfl, double dt, double dx);
int main()
{
	int i;
	int value = 1;  //选择case
	double cfl = 0.1; //CFL值
	double T_end = 3.5;   //时间长度
	double L = 20;    //空间长度
	double dx = 0.05;  //网格长度
	const int N = 400;   //网格数  N=L/dx
	double T = 0, dt, gamma = 1.4;
	double Q[N + 2][3];
	double F[N + 2][3];
	double count[N + 1];
	double k=1.0/3.0;
	double b=2.0;
	init(Q, gamma, N, value);
	while (T <= T_end)
	{
		dt = delta_t(Q, dx, N, cfl);
		T += dt;
		TVD(Q, F, N, cfl, dt, dx);
		bound(Q, N);
	}
	for (i = 0; i <= N; i++)
		count[i] = i * dx-L/2;
	for (i = 0; i <= N; i++)
	{
		cout << count[i] << "\t";
		cout << Q[i][0] << "\t";
		cout << Q[i][1] / Q[i][0] << "\t";
		cout << (gamma - 1) * (Q[i][2] - 0.5 * Q[i][1] * Q[i][1] / Q[i][0]) << endl;
	}
	return 0;
}







