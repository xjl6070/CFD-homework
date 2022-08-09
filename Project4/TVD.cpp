#include <cmath>
#include<vector>
#include <iostream>
using namespace std;
double Qk(double x);
double minmod(double a, double b);
void TVD(double Q[][3], double F[][3], int N, double cfl, double dt, double dx)
{
	int i, j, k;
	vector<double> rou(N + 2);
	vector<double> u(N + 2);
	vector<double> p(N + 2);
	vector<double> H(N + 2);
	vector<double> c(N + 2);
	vector<double> rou_a(N + 1);
	vector<double> u_a(N + 1);
	vector<double> p_a(N + 1);
	vector<double> H_a(N + 1);
	vector<double> c_a(N + 1);
	vector<vector<double>> lamda(N + 2, vector<double>(3));
	vector<vector<vector<double> > > L(N + 2, vector<vector<double>>(3, vector<double>(3))); //vector定义三维矩阵
	vector<vector<vector<double> > > R(N + 2, vector<vector<double>>(3, vector<double>(3)));
	vector<vector<double>> IQ(N + 2, vector<double>(3));
	vector<vector<double>> g_a(N + 2, vector<double>(3));
	vector<vector<double>> g(N + 2, vector<double>(3));
	vector<vector<double>> GAMMA(N + 2, vector<double>(3));
	vector<vector<double>> psai(N + 2, vector<double>(3));
	double sum;
	vector<vector<double>> sigma(N + 2, vector<double>(3));
	vector<vector<double>> NF(N + 2, vector<double>(3));
	double r = dt / dx;
	const double gamma = 1.4;
	for (i = 0; i <= N + 1; i++)
	{
		rou[i] = Q[i][0];
		if (rou[i] == 0)
			rou[i] = 1.0e-8;
		u[i] = Q[i][1] / rou[i];
		p[i] = (gamma - 1) * (Q[i][2] - 0.5 * Q[i][1] * u[i]);
		c[i] = sqrt(gamma * p[i] / rou[i]);
		H[i] = c[i] * c[i] / (gamma - 1) + 0.5 * u[i] * u[i];
	}
	for (i = 0; i <= N; i++)
	{
		rou_a[i] = sqrt(rou[i] * rou[i + 1]);
		u_a[i] = (u[i] * sqrt(rou[i]) + u[i + 1] * sqrt(rou[i + 1])) / (sqrt(rou[i]) + sqrt(rou[i + 1]));
		H_a[i] = (H[i] * sqrt(rou[i]) + H[i + 1] * sqrt(rou[i + 1])) / (sqrt(rou[i]) + sqrt(rou[i + 1]));
		c_a[i] = sqrt((gamma - 1) * (H_a[i] - 0.5 * u_a[i] * u_a[i]));
		p_a[i] = rou_a[i] * c_a[i] * c_a[i] / gamma;
	}
	for (i = 0; i <= N; i++)
	{
		lamda[i][0] = u_a[i] - c_a[i];
		lamda[i][1] = u_a[i];
		lamda[i][2] = u_a[i] + c_a[i];
		L[i][0][0] = (0.5 * u_a[i] * u_a[i] + u_a[i] * c_a[i] / (gamma - 1)) * (gamma - 1) / (2 * c_a[i] * c_a[i]);
		L[i][0][1] = (-u_a[i] - c_a[i] / (gamma - 1)) * (gamma - 1) / (2 * c_a[i] * c_a[i]);
		L[i][0][2] = (gamma - 1) / (2 * c_a[i] * c_a[i]);
		L[i][1][0] = (-u_a[i] * u_a[i] + 2 * c_a[i] * c_a[i] / (gamma - 1)) * (gamma - 1) / (2 * c_a[i] * c_a[i]);
		L[i][1][1] = 2 * u_a[i] * (gamma - 1) / (2 * c_a[i] * c_a[i]);
		L[i][1][2] = -2 * (gamma - 1) / (2 * c_a[i] * c_a[i]);
		L[i][2][0] = (0.5 * u_a[i] * u_a[i] - u_a[i] * c_a[i] / (gamma - 1)) * (gamma - 1) / (2 * c_a[i] * c_a[i]);
		L[i][2][1] = (-u_a[i] + c_a[i] / (gamma - 1)) * (gamma - 1) / (2 * c_a[i] * c_a[i]);
		L[i][2][2] = (gamma - 1) / (2 * c_a[i] * c_a[i]);
		R[i][0][0] = 1;
		R[i][0][1] = 1;
		R[i][0][2] = 1;
		R[i][1][0] = u_a[i] - c_a[i];
		R[i][1][1] = u_a[i];
		R[i][1][2] = u_a[i] + c_a[i];
		R[i][2][0] = c_a[i] * c_a[i] / (gamma - 1) + 0.5 * u_a[i] * u_a[i] - u_a[i] * c_a[i];
		R[i][2][1] = 0.5 * u_a[i] * u_a[i];
		R[i][2][2] = c_a[i] * c_a[i] / (gamma - 1) + 0.5 * u_a[i] * u_a[i] + u_a[i] * c_a[i];
	}
	for (i = 0; i <= N; i++)//求I*deltaQ
	{
		for (j = 0; j <= 2; j++)
		{
			sum = 0;
			for (k = 0; k <= 2; k++)
				sum = sum + L[i][j][k] * (Q[i + 1][k] - Q[i][k]);
			IQ[i][j] = sum;
		}
	}
	for (i = 0; i <= N; i++)
	{
		for (j = 0; j <= 2; j++)
			g_a[i][j] = 0.5 * (Qk(r * lamda[i][j]) - r * lamda[i][j] * r * lamda[i][j]) * IQ[i][j];
	}
	for (i = 1; i <= N; i++)
	{
		for (j = 0; j <= 2; j++)
			g[i][j] = minmod(g_a[i][j], g_a[i - 1][j]);
	}
	for (i = 1; i <= N - 1; i++)
	{
		for (j = 0; j <= 2; j++)
		{
			if (IQ[i][j] == 0)
				GAMMA[i][j] = 0;
			else
				GAMMA[i][j] = (g[i + 1][j] - g[i][j]) / IQ[i][j];
		}
	}
	for (i = 1; i <= N - 1; i++)
	{
		for (j = 0; j <= 2; j++)
		{
			psai[i][j] = ((g[i][j] + g[i + 1][j] - Qk(r * lamda[i][j] + GAMMA[i][j]) * IQ[i][j])) / r;
		}
	}
	for (i = 0; i <= N + 1; i++)//求 F 函数
	{
		F[i][0] = rou[i] * u[i];
		F[i][1] = rou[i] * u[i] * u[i] + p[i];
		F[i][2] = (p[i] / (gamma - 1) + 0.5 * rou[i] * u[i] * u[i] + p[i]) * u[i];
	}
	for (i = 1; i <= N - 1; i++)
	{
		for (j = 0; j <= 2; j++)
		{
			sum = 0;
			for (k = 0; k <= 2; k++)
				sum = sum + psai[i][k] * R[i][j][k];   //向量相加，对应元素相加
			NF[i][j] = 0.5 * (F[i][j] + F[i + 1][j] + sum);
		}
	}
	for (i = 2; i <= N - 1; i++)
	{
		for (j = 0; j <= 2; j++)
		{
			Q[i][j] = Q[i][j] - r * (NF[i][j] - NF[i - 1][j]);
		}
	}
}