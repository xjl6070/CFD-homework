#include <iostream>
using namespace std;
#define pi 3.1415926
#include <math.h>
int main()
{
	int i, j;
	double vp = 400.0; // �ɻ��ٶ�
	double vm = 600.0;//�����ٶ�
	double D = 10.0;//ɱ�˰뾶
	double yp = 10000 * sin(pi / 3);//�ɻ�������
	double xm[2000];
	double ym[2000];
	double xp[2000];
	xm[0] = 0;//��ʼ����������
	ym[0] = 0;//��ʼ����������
	xp[0] = 10000 * cos(pi / 3);//��ʼ�ɻ�������
	double T = 0.02;
	double q = atan((yp - ym[0]) / (xp[0] - xm[0]));
	double vmy = vm * sin(q), vmx = vm * cos(q);

	for (j = 1; ;j++)
	{
		xp[j] = xp[j - 1] + vp * T;
		xm[j] = xm[j - 1] + vmx * T;
		ym[j] = ym[j - 1] + vmy * T;
		q = atan((yp - ym[j]) / (xp[j] - xm[j]));
		vmy = 600 * sin(q);
		vmx = sqrt(600.0 * 600.0 - vmy * vmy);
		if (sqrt((xm[j] - xp[j]) * (xm[j] - xp[j]) + (ym[j] - yp) * (ym[j] - yp)) <= D)
		{
			cout << "q=" << q << endl;
			cout << "hit targets!!!" << endl;
			cout << "ʱ��Ϊ" << (xp[j] - 10000 * cos(pi / 3)) / vp << "s" << endl;
			break;
		}
	}
	for (int i = 0; i <= j; i++)
	{
		cout << xp[i] << "\t";
		cout << xm[i] << "\t";
		cout << ym[i] << "\t";
		cout << yp << endl;
	}
	return 0;
}