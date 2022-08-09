#include <cmath>
#include<vector>
using namespace std;
void MacCormack(double Q[][3], double F[][3], double Filter[][3], double Predictor[][3], 
	int N, double cfl,double beta, double dt, double dx, int Iswitch)
{
	int i, j;
	vector<double> rou(N + 2);
	vector<double> u(N + 2);
	vector<double> p(N + 2);
	double r = dt / dx;
	double niu, theta_rou, theta_u, theta_p;
	const double gamma = 1.4;
	const double eps = 1e-6;
	for (i = 0; i <= N + 1; i++)
	{
		rou[i] = Q[i][0];
		u[i] = Q[i][1] / Q[i][0];
		p[i] = (gamma - 1) * (Q[i][2] - 0.5 * Q[i][1] * u[i]);
	}
	for (i = 1; i <= N; i++)  //选择开关函数
	{
		theta_rou = fabs(rou[i + 1] - 2 * rou[i] + rou[i - 1]) / (fabs(rou[i + 1] - rou[i]) + fabs(rou[i] - rou[i - 1]) + eps);
		theta_u = fabs(u[i + 1] - 2 * u[i] + u[i - 1]) / (fabs(u[i + 1] - u[i]) + fabs(u[i] - u[i - 1]) + eps);
		theta_p = fabs(p[i + 1] - 2 * p[i] + p[i - 1]) / (fabs(p[i + 1] - p[i]) + fabs(p[i] - p[i - 1]) + eps);
		switch (Iswitch)
		{
		case 1:
			niu = cfl * (1 - cfl) * beta * theta_rou;
			break;
		case 2:
			niu = cfl * (1 - cfl) * beta * theta_u;
			break;
		case 3:
			niu = cfl * (1 - cfl) * beta * theta_p;
			break;
		}
		for (j = 0; j <= 2; j++)  //滤波
			Filter[i][j] = Q[i][j] + 0.5 * niu * (Q[i + 1][j] - 2 * Q[i][j] + Q[i - 1][j]);
	}
	for (j = 0; j <= 2; j++)
	{
		for (i = 1; i <= N; i++)
			Q[i][j] = Filter[i][j];
	}
		for (i = 0; i <= N + 1; i++)//求滤波结果的 F 函数
		{
			u[i] = Q[i][1] / Q[i][0];
			p[i] = (gamma - 1) * (Q[i][2] - 0.5 * Q[i][1] * u[i]);
			F[i][0] = Q[i][1];
			F[i][1] = Q[i][1] * u[i] + p[i];
			F[i][2] = (Q[i][2] + p[i]) * u[i];
		}
		for (i = 1; i <= N+1; i++) //预测
		{
			for (j = 0; j <= 2; j++)
				Predictor[i][j] = Q[i][j] - r * (F[i][j] - F[i-1][j]);
		}
	for (i = 0; i <= N; i++)  //求预测结果的 F 函数
	{
		u[i] = Predictor[i][1] / Predictor[i][0];
		p[i] = (gamma - 1) * (Predictor[i][2] - 0.5 * Predictor[i][1] * u[i]);
		F[i][0] = Predictor[i][1];
		F[i][1] = Predictor[i][1] * u[i] + p[i];
		F[i][2] = (Predictor[i][2] + p[i]) * u[i];
	}
	for (i = 1; i <= N; i++)  //修正
	{
		for (j = 0; j <= 2; j++)
			Q[i][j] = 0.5 * (Q[i][j] + Predictor[i][j]) - 0.5 * r * (F[i+1][j] - F[i][j]);
	}
}