#include <iostream>
#include <cmath>
using namespace std;
double delta_t(double Q[][82][4], double dx, int N,double cfl);
void init(double Q[][82][4], double gamma, int N, int value);
void bound(double Q[][82][4], int N ,double dx);
void TVD(double Q[][82][4], double F[][82][4], int N, double cfl, double dt, double dx);
int main()
{
	int i;
	int value = 1;  //ѡ��case
	double cfl = 0.1; //CFLֵ
	double T_end = 0.85;   //ʱ�䳤��
	double L_x = 2, L_y = 1;    //�ռ䳤��
	double dx = 0.025, dy =dx;  //���񳤶�
	const int N = 80;   //������  N=L/dx
	double T = 0, dt, gamma = 1.4;
	double Q[N + 2][N + 2][4];
	double F[N + 2][N + 2][4];
	double G[N + 2][N + 2][4];
	double count[N + 1];
	init(Q, gamma, N, value);
	while (T <= T_end)
	{
		dt = delta_t(Q, dx, N, cfl);
		T += dt;
		TVD(Q, F, N, cfl, dt, dx);
		bound(Q, N, dx);
	}
	//for (i = 0; i <= N; i++)
	//	count[i] = i * dx-L/2;
	/*for (i = 0; i <= N; i++)
	{
		cout << count[i] << "\t";
		cout << Q[i][0] << "\t";
		cout << Q[i][1] / Q[i][0] << "\t";
		cout << (gamma - 1) * (Q[i][2] - 0.5 * Q[i][1] * Q[i][1] / Q[i][0]) << endl;
	}*/
	return 0;
}







