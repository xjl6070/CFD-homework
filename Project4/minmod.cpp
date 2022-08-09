double minmod(double a, double b)
{
	if (a > 0 && b > 0)
	{
		if (a > b)
			return b;
		else
			return a;
	}
	else if (a < 0 && b < 0)
	{
		if (a > b)
			return a;
		else
			return b;
	}
	else
		return 0;
}