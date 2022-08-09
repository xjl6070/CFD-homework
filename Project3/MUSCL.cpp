#include <cmath>
#include<vector>
using namespace std;
void Roe(double Q[][3], double F[][3], int N, double cfl, double dt, double dx)
{
	int i, j;
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
	vector<vector<double>> niu(N + 2);
	for (i = 0; i < niu.size(); i++)
		niu[i].resize(3);
	vector<vector<double>> NF(N + 2);
	for (i = 0; i < NF.size(); i++)
		NF[i].resize(3);
	double r = dt / dx;
	double alpha1, alpha2, alpha3, alpha4, alpha5;
	const double gamma = 1.4;
	for (i = 0; i <= N + 1; i++)
	{
		rou[i] = Q[i][0];
		u[i] = Q[i][1] / Q[i][0];
		p[i] = (gamma - 1) * (Q[i][2] - 0.5 * Q[i][1] * u[i]);
		c[i] = sqrt(gamma * p[i] / rou[i]);
		H[i] = c[i] * c[i] / (gamma - 1) + 0.5 * u[i] * u[i];
	}
	for (i = 1; i <= N; i++)
	{
		rou_a[i] = sqrt(rou[i] * rou[i + 1]);
		u_a[i] = (u[i] * sqrt(rou[i]) + u[i + 1] * sqrt(rou[i + 1])) / (sqrt(rou[i]) + sqrt(rou[i + 1]));
		H_a[i] = (H[i] * sqrt(rou[i]) + H[i + 1] * sqrt(rou[i + 1])) / (sqrt(rou[i]) + sqrt(rou[i + 1]));
		c_a[i] = sqrt((gamma - 1) * (H_a[i] - 0.5 * u_a[i] * u_a[i]));
		p_a[i] = rou_a[i] * c_a[i] * c_a[i] / gamma;
		alpha1 = fabs(u_a[i]) * (rou[i + 1] - rou[i] - (p[i + 1] - p[i]) / (c_a[i] * c_a[i]));
		alpha2 = 0.5 / (c_a[i] * c_a[i]) * fabs(u_a[i] + c_a[i]) * (p[i + 1] - p[i] + rou_a[i] * c_a[i] * (u[i + 1] - u[i]));
		alpha3 = 0.5 / (c_a[i] * c_a[i]) * fabs(u_a[i] - c_a[i]) * (p[i + 1] - p[i] - rou_a[i] * c_a[i] * (u[i + 1] - u[i]));
		alpha4 = alpha1 + alpha2 + alpha3;
		alpha5 = c_a[i] * (alpha2 - alpha3);
		niu[i][0] = alpha4;
		niu[i][1] = u_a[i] * alpha4 + alpha5;
		niu[i][2] = H_a[i] * alpha4 + u_a[i] * alpha5 - c_a[i] * c_a[i] * alpha1 / (gamma - 1);
	}
	for (i = 0; i <= N + 1; i++)//Çó F º¯Êý
	{
		F[i][0] = rou[i] * u[i];;
		F[i][1] = rou[i] * u[i] * u[i] + p[i];
		F[i][2] = (p[i] / (gamma - 1) + 0.5 * rou[i] * u[i] * u[i] + p[i]) * u[i];
	}
	for (j = 0; j <= 2; j++)
	{
		for (i = 0; i <= N; i++)
			NF[i][j] = 0.5 * (F[i][j] + F[i + 1][j]) - 0.5 * niu[i][j];
	}
	for (j = 0; j <= 2; j++)
	{
		for (i = 1; i <= N; i++)
		{
			Q[i][j] = Q[i][j] - r * (NF[i][j] - NF[i - 1][j]);
		}
	}
}