void init(double Q[][3], double gamma, int N,int value)
{
	int i;
	double rou1, u1, p1, rou2, u2, p2;
	switch (value)
	{
	case 1://�������-�Ӵ����-���Ͳ�����֮Sod����
		rou1 = 1.0; u1 = 0.0; p1 = 1.0;
		rou2 = 0.125; u2 = 0.0; p2 = 0.1;
		break;
	case 2://�������-�Ӵ����-���Ͳ�����֮Lax����
		rou1 = 0.445; u1 = 0.698; p1 = 3.528;
		rou2 = 0.5; u2 = 0.0; p2 = 0.571;
		break;
	case 3://˫���Ͳ�����֮�����ٺ�������
		rou1 = 1.0; u1 = -2.0; p1 = 4.0;
		rou2 = 1.0; u2 = 2.0; p2 = 4.0;
		break;
	case 4://˫���Ͳ�����֮Sjogreen supersonic Expansion����   T_end=1.5
		rou1 = 1.0; u1 = -2.0; p1 = 0.4;
		rou2 = 1.0; u2 = 2.0; p2 = 0.4;
		break;
	case 5://�Ӵ����-˫���Ͳ�����
		rou1 = 1.0; u1 = -0.2; p1 = 0.5;
		rou2 = 0.5; u2 = 0.5; p2 = 0.5;
		break;
	case 6://�Ӵ����-˫�����������
		rou1 = 0.4; u1 = 0.5; p1 = 1.0;
		rou2 = 1.0; u2 = -0.5; p2 = 0.9;
		break;
	case 7://�Ӵ��������
		rou1 = 10; u1 = 1.0; p1 = 2.0;
		rou2 = 1.0; u2 = 1.0; p2 = 2.0;
		break;
	}	
	for (i = 0; i <= N / 2; i++)
	{
		Q[i][0] = rou1;
		Q[i][1] = rou1 * u1;
		Q[i][2] = p1 / (gamma - 1) + rou1 * u1 * u1 / 2;
	}
	for (i = N / 2 + 1; i <= N + 1; i++)
	{
		Q[i][0] = rou2;
		Q[i][1] = rou2 * u2;
		Q[i][2] = p2 / (gamma - 1) + rou2 * u2 * u2 / 2;
	}
}