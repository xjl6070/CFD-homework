#include<cmath>
double Qk(double x)
{
	const double epsilon = 0.1;
	double Qk;
	if (fabs(x) >= epsilon)
		Qk = fabs(x);
	else
		Qk = (x * x + epsilon * epsilon) / (2 * epsilon);
	return Qk;
}