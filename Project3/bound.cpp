void bound(double Q[][3], int N)
{
	int j;
	for (j = 0; j <= 2; j++)
		Q[0][j] = Q[1][j];
	for (j = 0; j <= 2; j++)
		Q[N + 1][j] = Q[N][j];
}