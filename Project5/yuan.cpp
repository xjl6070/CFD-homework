#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define GAMA 1.4//���峣��
#define PI   3.141592654
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

#define Lx 4.0//��������
#define Ly 1.0

#define TT 3.0//��ʱ��
#define Sf 0.8//ʱ�䲽������

#define Jx 400//������
#define Jy 100

//ȫ�ֱ���
double U[Jx + 2][Jy + 2][4], Uf[Jx + 2][Jy + 2][4], EFf[Jx + 2][Jy + 2][4];

/*-------------------------------------------------------
����ʱ�䲽��
���: U�� ��ǰ��������dx��dy�� �����ȣ�
����: ʱ�䲽����
---------------------------------------------------------*/
double CFL(double U[Jx + 2][Jy + 2][4], double dx, double dy)
{
    int i, j;
    double maxvel, p, u, v, vel;
    maxvel = 1e-100;
    for (i = 1; i <= Jx; i++)
        for (j = 1; j <= Jy; j++)
        {
            u = U[i][j][1] / U[i][j][0];
            v = U[i][j][2] / U[i][j][0];
            p = (GAMA - 1) * (U[i][j][3] - 0.5 * U[i][j][0] * (u * u + v * v));
            vel = sqrt(GAMA * p / U[i][j][0]) + sqrt(u * u + v * v);
            if (vel > maxvel)maxvel = vel;
        }
    return Sf * MIN(dx, dy) / maxvel;
}

/*-----------------------------------------------------------------
��ʼ��
���: �ޣ�
����: U�� �Ѿ������ĳ�ʼֵ��dx��dy�� �����ȡ�
--------------------------------------------------------------------*/
void Init(double U[Jx + 2][Jy + 2][4], double& dx, double& dy)//��ʼ��
{
    int i, j;
    double rou1 = 1.0, u1 = 2.9, v1 = 0, p1 = 0.71429; //��ʼ����
    double rou2 = 1.69997, u2 = 2.61934, v2 = -0.50632, p2 = 1.52819;
    double theta, lx;

    dx = Lx / Jx;
    dy = Ly / Jy;

    theta = 29 * PI / 180;//���伤���Ƕ� 

    for (j = 0; j <= Jy + 1; j++)
        for (i = 0; i <= Jx + 1; i++)
        {
            lx = (1 - j * dy) / tan(theta);
            if (i * dx <= lx)
            {
                U[i][j][0] = rou1;
                U[i][j][1] = rou1 * u1;
                U[i][j][2] = rou1 * v1;
                U[i][j][3] = p1 / (GAMA - 1) + rou1 * (u1 * u1 + v1 * v1) / 2;
            }
            else
            {
                U[i][j][0] = rou2;
                U[i][j][1] = rou2 * u2;
                U[i][j][2] = rou2 * v2;
                U[i][j][3] = p2 / (GAMA - 1) + rou2 * (u2 * u2 + v2 * v2) / 2;
            }
        }
}



/*-------------------------------------------------------
�߽�����
���: dx��dy�� �����ȣ�
����: U�� �Ѿ������߽硣
---------------------------------------------------------*/
void bound(double U[Jx + 2][Jy + 2][4], double dx, double dy)
{
    int i, j, k;
    double rou1 = 1.0, u1 = 2.9, v1 = 0, p1 = 0.71429; //��ʼ����
    double rou2 = 1.69997, u2 = 2.61934, v2 = -0.50632, p2 = 1.52819;

    //��߽�
    for (j = 0; j <= Jy + 1; j++)
    {
        U[0][j][0] = rou1;
        U[0][j][1] = rou1 * u1;
        U[0][j][2] = rou1 * v1;
        U[0][j][3] = p1 / (GAMA - 1) + rou1 * (u1 * u1 + v1 * v1) / 2;
    }

    //�ұ߽� 
    for (j = 0; j <= Jy + 1; j++)
        for (k = 0; k < 4; k++)
        {
            U[Jx + 1][j][k] = U[Jx][j][k];
        }

    //�ϱ߽�
    for (i = 0; i <= Jx + 1; i++)
    {
        U[i][Jy + 1][0] = rou2;
        U[i][Jy + 1][1] = rou2 * u2;
        U[i][Jy + 1][2] = rou2 * v2;
        U[i][Jy + 1][3] = p2 / (GAMA - 1) + rou2 * (u2 * u2 + v2 * v2) / 2;
    }

    //�±߽�
    for (i = 0; i <= Jx + 1; i++)
    {
        U[i][0][0] = U[i][1][0];
        U[i][0][1] = U[i][1][1];
        U[i][0][2] = 0;
        U[i][0][3] = U[i][1][3];
    }
}

/*-------------------------------------------------------
����U����E
���: U�� ��ǰUʸ����
����: E�� ����õ���Eʸ����
U��E�������Euler�����顣
---------------------------------------------------------*/
void U2E(double U[4], double E[4])
{
    double u, v, p;
    u = U[1] / U[0];
    v = U[2] / U[0];
    p = (GAMA - 1) * (U[3] - 0.5 * (U[1] * U[1] + U[2] * U[2]) / U[0]);
    E[0] = U[1];
    E[1] = U[0] * u * u + p;
    E[2] = U[0] * u * v;
    E[3] = (U[3] + p) * u;
}

/*-------------------------------------------------------
����U����F
���: U�� ��ǰUʸ����
����: E�� ����õ���Fʸ����
U��F�����Euler�����顣
---------------------------------------------------------*/
void U2F(double U[4], double F[4])//����U����F
{
    double u, v, p;
    u = U[1] / U[0];
    v = U[2] / U[0];
    p = (GAMA - 1) * (U[3] - 0.5 * (U[1] * U[1] + U[2] * U[2]) / U[0]);
    F[0] = U[2];
    F[1] = U[0] * u * v;
    F[2] = U[0] * v * v + p;
    F[3] = (U[3] + p) * v;
}

/*-------------------------------------------------------
��ά ��ָ�ʽ�� ����Ĳ������
���: U����һʱ��Uʸ����Uf��Ef ��ʱ����
     dx��x����������ȣ�dt�� ʱ�䲽����
����: U������õ��ĵ�ǰʱ��Uʸ����
---------------------------------------------------------*/
void LLx(double U[Jx + 2][Jy + 2][4], double Uf[Jx + 2][Jy + 2][4], double Ef[Jx + 2][Jy + 2][4], double dx, double dt)
{
    int i, j, k;
    double r, nu, q;

    r = dt / dx;
    nu = 3.0 * r * (1 - 3.0 * r);

    for (i = 1; i <= Jx; i++)
        for (j = 0; j <= Jy + 1; j++)
        {
            q = fabs(fabs(U[i + 1][j][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i - 1][j][0]))
                / (fabs(U[i + 1][j][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i - 1][j][0]) + 1e-100); //���غ���
            for (k = 0; k < 4; k++)
                Ef[i][j][k] = U[i][j][k] + 0.5 * nu * q * (U[i + 1][j][k] - 2 * U[i][j][k] + U[i - 1][j][k]);//�˹������ 
        }
    for (j = 0; j <= Jy + 1; j++)
        for (k = 0; k < 4; k++)
            for (i = 1; i <= Jx; i++)U[i][j][k] = Ef[i][j][k];

    for (i = 0; i <= Jx + 1; i++)
        for (j = 0; j <= Jy + 1; j++)
            U2E(U[i][j], Ef[i][j]);

    for (i = 0; i <= Jx; i++)
        for (j = 0; j <= Jy + 1; j++)
            for (k = 0; k < 4; k++)
                Uf[i][j][k] = 0.5 * (U[i + 1][j][k] + U[i][j][k]) - 0.5 * r * (Ef[i + 1][j][k] - Ef[i][j][k]); //U(n+1/2)(i+1/2)(j)

    for (i = 0; i <= Jx; i++)
        for (j = 0; j <= Jy + 1; j++)
            U2E(Uf[i][j], Ef[i][j]); //E(n+1/2)(i+1/2)(j)

    for (i = 1; i <= Jx; i++)
        for (j = 0; j <= Jy + 1; j++)
            for (k = 0; k < 4; k++)
                U[i][j][k] = U[i][j][k] - r * (Ef[i][j][k] - Ef[i - 1][j][k]); //U(n+1)(i)(j)  
}

/*-------------------------------------------------------
��ά ��ָ�ʽ��y�����ϵĲ������
���: U����һʱ��Uʸ����Uf��Ff����ʱ������
     Dy��y�����������ȣ�dt��ʱ�䲽����
����: U������õ��ĵ�ǰʱ��Uʸ����
---------------------------------------------------------*/
void LLy(double U[Jx + 2][Jy + 2][4], double Uf[Jx + 2][Jy + 2][4], double Ff[Jx + 2][Jy + 2][4], double dy, double dt)
{
    int i, j, k;
    double r, nu, q;

    r = dt / dy;

    nu = 3.0 * r * (1 - 3.0 * r);
    //nu=2.0;

    for (i = 0; i <= Jx + 1; i++)
        for (j = 1; j <= Jy; j++)
        {
            q = fabs(fabs(U[i][j + 1][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i][j - 1][0]))
                / (fabs(U[i][j + 1][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i][j - 1][0]) + 1e-100);
            for (k = 0; k < 4; k++)
                Ff[i][j][k] = U[i][j][k] + 0.5 * nu * q * (U[i][j + 1][k] - 2 * U[i][j][k] + U[i][j - 1][k]);//�˹������ 
        }
    for (i = 0; i <= Jx + 1; i++)
        for (k = 0; k < 4; k++)
            for (j = 1; j <= Jy; j++)U[i][j][k] = Ff[i][j][k];

    for (i = 0; i <= Jx + 1; i++)
        for (j = 0; j <= Jy + 1; j++)
            U2F(U[i][j], Ff[i][j]);

    for (i = 0; i <= Jx + 1; i++)
        for (j = 0; j <= Jy; j++)
            for (k = 0; k < 4; k++)
                Uf[i][j][k] = 0.5 * (U[i][j + 1][k] + U[i][j][k]) - 0.5 * r * (Ff[i][j + 1][k] - Ff[i][j][k]); //U(n+1/2)(i)(j+1/2)

    for (i = 0; i <= Jx + 1; i++)
        for (j = 0; j <= Jy; j++)
            U2F(Uf[i][j], Ff[i][j]); //F(n+1/2)(i)(j+1/2)

    for (i = 0; i <= Jx + 1; i++)
        for (j = 1; j <= Jy; j++)
            for (k = 0; k < 4; k++)
                U[i][j][k] = U[i][j][k] - r * (Ff[i][j][k] - Ff[i][j - 1][k]); //U(n+1)(i)(j)
}

/*-------------------------------------------------------
��ά ��ָ�ʽ�����
���: U����һʱ��Uʸ����Uf��EFf����ʱ������
     Dx��x����������ȣ�dy��y�����ϵ������ȣ�
     Dt�� ʱ�䲽����
����: U�� ����õ��ĵ�ǰʱ��Uʸ����
---------------------------------------------------------*/
void Lax_Wendroff_2DSolver(double U[Jx + 2][Jy + 2][4], double Uf[Jx + 2][Jy + 2][4], double EFf[Jx + 2][Jy + 2][4], double dx, double dy, double dt)//Lax_Wendroff_2D�����
{
    LLx(U, Uf, EFf, dx, dt / 2.0);
    bound(U, dx, dy);
    LLy(U, Uf, EFf, dy, dt / 2.0);
    bound(U, dx, dy);
    LLy(U, Uf, EFf, dy, dt / 2.0);
    bound(U, dx, dy);
    LLx(U, Uf, EFf, dx, dt / 2.0);
    bound(U, dx, dy);
}

/*-------------------------------------------------------
��� �ļ���ʽ���ݻ�ͼ
���: U����ǰʱ��Uʸ����
     dx�� x����������ȣ�dy�� y����������ȣ�
����: �ޡ�
---------------------------------------------------------*/
void Output(double U[Jx + 2][Jy + 2][4], double dx, double dy)
{
    int i, j;
    FILE* fp;
    double rou, u, v, p;

    fp = fopen("result.plt", "w");
    fprintf(fp, "TITLE     = \"Dataset\"\nVARIABLES = \"x\" \"y\" \"rou\" \"u\" \"v\" \"p\" \"E\"");
    fprintf(fp, "ZONE T=\"Zone 1\"\nI=%d J=%d K=%d ZONETYPE=Ordered\n", Jy + 2, Jx + 2, 1);
    fprintf(fp, "DATAPACKING=POINT\n");

    for (i = 0; i <= Jx + 1; i++)
        for (j = 0; j <= Jy + 1; j++)
        {
            rou = U[i][j][0];
            u = U[i][j][1] / rou;
            v = U[i][j][2] / rou;
            p = (GAMA - 1) * (U[i][j][3] - 0.5 * U[i][j][0] * (u * u + v * v));
            fprintf(fp, "%20f%20f%20.10e%20.10e%20.10e%20.10e%20.10e\n", i * dx, j * dy, rou, u, v, p, U[i][j][3]);
        }
    fclose(fp);

    fp = fopen("result.txt", "w");
    j = Jy / 2;
    i = Jx / 2;
    for (i = 0; i <= Jx + 1; i++)
        //for(j=0;j<=Jy+1;j++)
    {
        rou = U[i][j][0];
        u = U[i][j][1] / rou;
        v = U[i][j][2] / rou;
        p = (GAMA - 1) * (U[i][j][3] - 0.5 * U[i][j][0] * (u * u + v * v));
        fprintf(fp, "%20f%20.10e%20.10e%20.10e%20.10e%20.10e\n", i * dx, rou, u, v, p, U[i][j][3]);
    }
    fclose(fp);
}

/*-------------------------------------------------------
������
���: �ޣ�
����: �ޡ�
---------------------------------------------------------*/
void main()
{
    double T, dx, dy, dt;

    Init(U, dx, dy);
    T = 0;
    while (T < TT)
    {
        dt = CFL(U, dx, dy);
        T += dt;
        printf("T=%10g    dt=%10g\n", T, dt);
        Lax_Wendroff_2DSolver(U, Uf, EFf, dx, dy, dt);
    }
    Output(U, dx, dy);
}
