void init(double Q[][4], double gamma, int N,int value)
{
	int i;
	double rou2, u2, v2, p2;
	switch (value)
	{
	case 1:
		rou2 = 8.0 / 3.0; u2 = 1.25; v2 = 0; p2 = 45.0 / 14.0;
		break;
	case 2:
		rou2 = 27.0 / 7.0; u2 = 20.0 / 9.0; v2 = 0; p2 = 155.0 / 21.0;
		break;
	}	
	for (i = 0; i <= N / 2; i++)
	{
		
	}
	for (i = N / 2 + 1; i <= N + 1; i++)
	{
		
	}
}