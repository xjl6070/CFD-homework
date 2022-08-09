#include <iostream>
#include <cmath>
using namespace std;
double delta_t(double Q[][3], double dx, int N,double cfl);
void init(double Q[][3], double gamma, int N, int value);
void bound(double Q[][3], int N);
void MacCormack(double Q[][3], double F[][3], double Filter[][3], double Predictor[][3],
	 int N, double cfl, double beta, double dt, double dx, int Iswitch);
int main()
{
	int i;
	int value = 4;  //ѡ��case
	int Iswitch = 1; //ѡ��theta��ⷽʽ
	double cfl = 0.1; //CFLֵ
	double beta = 1.0;   //��ֵ �ı��˹�ճ��
	double T_end = 3.5;   //ʱ�䳤��
	double L = 40;    //�ռ䳤��
	double dx = 0.05;  //���񳤶�
	const int N = 800;   //������  N=L/dx
	double T = 0, dt, gamma = 1.4;
	double Q[N + 2][3];
	double F[N + 2][3];
	double Filter[N + 2][3];
	double Predictor[N + 2][3];
	double count[N + 1];
	init(Q, gamma, N, value);
	while (T <= T_end)
	{
		dt = delta_t(Q, dx, N, cfl);
		T += dt;
		MacCormack(Q, F, Filter, Predictor, N, cfl, beta, dt, dx, Iswitch);
		bound(Q, dx);
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







