#include <cmath>
double delta_t(double Q[][3], double dx, int N,double cfl)
{
	int i;
	double max_v, p, u, a, gamma = 1.4;
	max_v = 1e-100;
	for (i = 1; i <= N; i++)
	{
		u = Q[i][1] / Q[i][0];
		p = (gamma - 1) * (Q[i][2] - 0.5 * Q[i][1] * u);
		a = sqrt(gamma * p / Q[i][0]);
		if (a + fabs(u) > max_v)
			max_v = a + fabs(u);
	}
	return cfl * dx / max_v;
}