#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
int main()
{
	int a = 1;
	int T;   //初始函数分布形式
	int N;   //选择格式
	double x = 8.0;  //空间长度
	double t = 1.0;  //时间长度
	double c;        //CFL
	double delta_x;
	cout << "初始函数分布形式：";
	cin >> T;
	cout << "格式：";
	cin >> N;
	cout << "delta_x:";
	cin >> delta_x;
	cout << "CFL=";
	cin >> c;
	int xn = x / delta_x;
	double delta_t = c * delta_x / a;
	double tn = t / delta_t;
	vector<double>u(xn + 1);
	vector<double>u_init(xn + 1);
	vector<double>uo(xn + 1);
	vector<double>u_end(xn + 1);
	vector<double>X(xn + 1);
	vector<vector<double>>data(xn+1, vector<double>(4));
	X[0] = - 1;
	double pi = 3.1415926;

	for (int i = 0; i <= xn - 1; i++)
	{
		if (i >= 1.0 / delta_x && i <= 2.0 / delta_x)
			switch (T)
			{
			case 1:
				u[i] = 1; //方波
				break;
			case 2:
				u[i] = exp(-16 * (i * delta_x - 1 - 0.5) * (i * delta_x - 1 - 0.5))
					* sin(40 * pi * (i * delta_x - 1));//振荡 
				break;
			case 3:
				u[i] = exp(-16 * (i * delta_x - 1 - 0.5) * (i * delta_x - 1 - 0.5))
					* (-64 * (i * delta_x - 1) * (i * delta_x - 1) * (i * delta_x - 1)
						* (i * delta_x - 2) * (i * delta_x - 2) * (i * delta_x - 2));//光滑
				break;
			default:cout << "Wrong.\n";
			}
		else
			u[i] = 0;
	}
	for (int i = 0; i <= xn - 1; i++)
		u_init[i] = u[i];

	//t=2精确波形
	for (int i = 0; i <= xn - 1; i++)
	{
		if (i >= (t + 1) / delta_x && i <= (t + 2) / delta_x)
		{
			switch (T)
			{
			case 1:
				u_end[i] = 1;
				break;
			case 2:
				u_end[i] = exp(-16 * (i * delta_x - t - 1 - 0.5) * (i * delta_x - t - 1 - 0.5)) 
					* sin(40 * pi * (i * delta_x - t - 1));
				break;
			case 3:
				u_end[i] = exp(-16 * (i * delta_x - (t + 1) - 0.5) * (i * delta_x - (t + 1) - 0.5))
					* (-64 * (i * delta_x - (t + 1)) * (i * delta_x - (t + 1)) * (i * delta_x - (t + 1))
						* (i * delta_x - (t + 1) - 1) * (i * delta_x - (t + 1) - 1) * (i * delta_x - (t + 1) - 1));
				break;
			default:cout << "Wrong.\n";
			}
		}
		else
			u_end[i] = 0;
	}
	for (int j = 1; j <= tn; j++)
	{
		for (int n = 0; n <= xn - 1; n++)
			uo[n] = u[n];
		switch (N)
		{
		case 1:
			for (int i = 1; i <= xn - 1; i++)
				u[i] = (1 - c) * uo[i] + c * uo[i - 1];  //FTBS格式
			break;
		case 2:
			for (int i = 1; i <= xn - 1; i++)
				u[i] = (0.5 - 0.5 * c) * uo[i + 1] + (0.5 + 0.5 * c) * uo[i - 1];  //Lax格式
			break;
		case 3:
			for (int i = 1; i <= xn - 1; i++)
				u[i] = (0.5 * c * c - 0.5 * c) * uo[i + 1] + (1 - c * c) * uo[i]
				+ (0.5 * c * c + 0.5 * c) * uo[i - 1];  //LaxWendroff格式
			break;
		case 4:
			for (int i = 2; i <= xn - 1; i++)
				u[i] = (1 - 1.5 * c + 0.5 * c * c) * uo[i] + (2 * c - c * c) * uo[i - 1]
				+ (0.5 * c * c - 0.5 * c) * uo[i - 2];  //Warming-beam格式
			break;
		default:cout << "Wrong.\n";
		}
	}
	for (int i = 0; i <= xn - 1; i++)
		X[i] = X[0] + i * delta_x;
	for (int j = 0; j <= 3; j++)
	{
		for (int i = 0; i <= xn - 1; i++)
		{
			if (j == 0)
				data[i][j] = X[i];
			if (j == 1)
				data[i][j] = u[i];
			if (j == 2)
				data[i][j] = u_init[i];
			if (j == 3)
				data[i][j] = u_end[i];
		}
	}
	for (int i = 0; i <= xn -1; i++)
	{
		for (int j = 0; j <= 3; j++)
			cout << data[i][j] << "\t";
		cout << "\n";
	}
	system("pause");
	return 0;
}

