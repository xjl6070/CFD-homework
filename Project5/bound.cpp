void bound(double Q[][3], int N ,double dx)
{
	int i, j, k;
	for (j = 0; j <= 1.15/dx; j++)
	{
		for (k = 0; k <= 3; k++)
		{
			Q[j + 1][k] = Q[j][k];
		}
	}
	for (j = 0; j <= 2; j++)
	{
		Q[N - 1][j] = Q[N - 2][j];
		Q[N][j] = Q[N - 1][j];

	}
}