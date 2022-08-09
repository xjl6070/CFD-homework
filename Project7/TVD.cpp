/*----------------------------------------------------------------------------------------------
*Harten二阶TVD格式求解二维平面激波在刚性斜面反射问题*
-------------------------------------------------------------------------------------------------*/
//#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
#define Ms 2.0 //马赫数
#define gamma 1.4//气体常数
#define PI   3.141592654
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

#define Lx 2.0//计算区域
#define Ly 1.0

#define TT 0.85//总时间
#define cfl 0.1 //CFL数

#define Nx 400//网格数
#define Ny 200

//全局变量
double Q[Nx + 2][Ny + 2][4], F[Nx + 2][Ny + 2][4], G[Nx + 2][Ny + 2][4];

/*-------------------------------------------------------
计算时间步长
入口: Q， 当前物理量，dx、dy， 网格宽度；
返回: 时间步长。
---------------------------------------------------------*/
double CFL(double Q[Nx + 2][Ny + 2][4], double dx, double dy)
{
    int i, j;
    double maxvel, p, u, v, vel;
    maxvel = 1e-100;
    for (i = 1; i <= Nx + 1; i++)
        for (j = 1; j <= Ny + 1; j++)
        {
            u = Q[i][j][1] / Q[i][j][0];
            v = Q[i][j][2] / Q[i][j][0];
            p = (gamma - 1) * (Q[i][j][3] - 0.5 * Q[i][j][0] * (u * u + v * v));
            //vel = sqrt(gamma * p / Q[i][j][0]) + sqrt(u * u + v * v);
            vel = MAX(sqrt(gamma * p / Q[i][j][0]) + u, sqrt(gamma * p / Q[i][j][0]) + v);
            if (vel > maxvel)maxvel = vel;
        }
    return cfl * MIN(dx, dy) / maxvel;
}

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
/*-----------------------------------------------------------------
初始化
--------------------------------------------------------------------*/
void Init(double Q[Nx + 2][Ny + 2][4], double& dx, double& dy)
{
    int i, j;
    double rou1 = 1.0, u1 = 0, v1 = 0, p1 = 5.0 / 7.0; //初始条件
    double rou2 = rou1 * Ms * Ms / (1 + (gamma - 1) / (gamma + 1) * (Ms * Ms - 1));
    double u2 = u1 + 2 / (gamma + 1) * (Ms - 1.0 / Ms);
    double v2 = 0;
    double p2 = p1 * (1 + 2 * gamma / (gamma + 1) * (Ms * Ms - 1));
    dx = Lx / Nx;
    dy = Ly / Ny;
    for (j = 0; j <= Ny + 1; j++)
    {
        Q[0][j][0] = rou2;
        Q[0][j][1] = rou2 * u2;
        Q[0][j][2] = rou2 * v2;
        Q[0][j][3] = p2 / (gamma - 1) + rou2 * (u2 * u2 + v2 * v2) / 2;
        Q[1][j][0] = Q[0][j][0];
        Q[1][j][1] = Q[0][j][1];
        Q[1][j][2] = Q[0][j][2];
        Q[1][j][3] = Q[0][j][3];
    }
    for (j = 0; j <= Ny + 1; j++)
    {
        for (i = 2; i <= Nx + 1; i++)
        {
            Q[i][j][0] = rou1;
            Q[i][j][1] = rou1 * u1;
            Q[i][j][2] = rou1 * v1;
            Q[i][j][3] = p1 / (gamma - 1) + rou1 * (u1 * u1 + v1 * v1) / 2;
        }
    }
}
/*-------------------------------------------------------
边界条件
---------------------------------------------------------*/
void bound(double Q[Nx + 2][Ny + 2][4], double dx, double dy)
{
    int i, j, k;
    double rou1 = 1.0, u1 = 0, v1 = 0, p1 = 5.0 / 7.0; //初始条件
    double rou2 = rou1 * Ms * Ms / (1 + (gamma - 1) / (gamma + 1) * (Ms * Ms - 1));
    double u2 = u1 + 2 / (gamma + 1) * (Ms - 1.0 / Ms);
    double v2 = 0;
    double p2 = p1 * (1 + 2 * gamma / (gamma + 1) * (Ms * Ms - 1));
    //下边界
    for (i = 0; i <= (float)(1.15 / dx); i++)
    {
        Q[i][1][0] = Q[i][2][0];
        Q[i][1][1] = Q[i][2][1];
        Q[i][1][2] = 0;
        Q[i][1][3] = Q[i][2][3];
    }
    for (i = 0; i <= (float)(1.15 / dx); i++)
    {
        for (k = 0; k <= 3; k++)
            Q[i][0][k] = Q[i][1][k];
    }
    //楔面 
    for (j = 1; j <= (float)(0.85 / dy); j++)
    {
        i = j + (float)(1.15 / dx);
        k = (float)(1.15 / dx) + 1;
        //第一层
        Q[i][j - 1][0] = Q[i - 1][j][0];
        Q[i][j - 1][1] = Q[i - 1][j][2];
        Q[i][j - 1][2] = Q[i - 1][j][1];
        Q[i][j - 1][3] = Q[i - 1][j][3];
        Q[Nx + 1][k - 1][0] = Q[Nx][k][0];
        Q[Nx + 1][k - 1][1] = Q[Nx][k][2];
        Q[Nx + 1][k - 1][2] = Q[Nx][k][1];
        Q[Nx + 1][k - 1][3] = Q[Nx][k][3];
        //第二层
        Q[i + 1][j - 1][0] = Q[i - 1][j + 1][0];
        Q[i + 1][j - 1][1] = Q[i - 1][j + 1][2];
        Q[i + 1][j - 1][2] = Q[i - 1][j + 1][1];
        Q[i + 1][j - 1][3] = Q[i - 1][j + 1][3];
    }
    for (j = 0; j <= (float)(0.85 / dx) + 1; j++)
    {
        i = j + (float)(1.15 / dx);
        //壁面
        Q[i][j][0] = Q[i - 1][j + 1][0];
        Q[i][j][1] = (Q[i - 1][j + 1][1] + Q[i - 1][j + 1][2]) / 2;
        Q[i][j][2] = (Q[i - 1][j + 1][2] + Q[i - 1][j + 1][1]) / 2;
        Q[i][j][3] = Q[i - 1][j + 1][3];
    }
    //右上边界
    for (j = (float)(0.85 / dx) + 2; j <= Ny + 1; j++)
    {
        for (k = 0; k <= 3; k++)
            Q[Nx + 1][j][k] = Q[Nx][j][k];
    }
    //上边界
    for (i = 0; i <= Nx + 1; i++)
    {
        for (k = 0; k <= 3; k++)
            Q[i][Ny + 1][k] = Q[i][Ny][k];
    }
}
/*-------------------------------------------------------
二维差分格式
---------------------------------------------------------*/
void TVD(double Q[Nx + 2][Ny + 2][4], double dx, double dy, double dt)
{
    int i, j, k, m;
    double sum;
    double cx = dt / dx, cy = dt / dy;
    double rou[Nx + 2][Ny + 2];
    double u[Nx + 2][Ny + 2];
    double v[Nx + 2][Ny + 2];
    double p[Nx + 2][Ny + 2];
    double c[Nx + 2][Ny + 2];
    double H[Nx + 2][Ny + 2];
    double rou_ax[Nx + 2][Ny + 2], rou_ay[Nx + 2][Ny + 2];
    double u_ax[Nx + 2][Ny + 2], u_ay[Nx + 2][Ny + 2];
    double v_ax[Nx + 2][Ny + 2], v_ay[Nx + 2][Ny + 2];
    double p_ax[Nx + 2][Ny + 2], p_ay[Nx + 2][Ny + 2];
    double c_ax[Nx + 2][Ny + 2], c_ay[Nx + 2][Ny + 2];
    double H_ax[Nx + 2][Ny + 2], H_ay[Nx + 2][Ny + 2];
    double lamda_A[Nx + 2][Ny + 2][4];
    double lamda_B[Nx + 2][Ny + 2][4];
    double LA[Nx + 2][Ny + 2][4][4];
    double RA[Nx + 2][Ny + 2][4][4];
    double LB[Nx + 2][Ny + 2][4][4];
    double RB[Nx + 2][Ny + 2][4][4];
    double IAQ[Nx + 2][Ny + 2][4];   //I*Q
    double IBQ[Nx + 2][Ny + 2][4];   //I*Q
    double g_aB[Nx + 2][Ny + 2][4];
    double g_B[Nx + 2][Ny + 2][4];
    double g_aA[Nx + 2][Ny + 2][4];
    double g_A[Nx + 2][Ny + 2][4];
    double GAMMA_B[Nx + 2][Ny + 2][4];
    double GAMMA_A[Nx + 2][Ny + 2][4];
    double psai_A[Nx + 2][Ny + 2][4];
    double psai_B[Nx + 2][Ny + 2][4];
    double NF[Nx + 2][Ny + 2][4];
    double NG[Nx + 2][Ny + 2][4];
    for (i = 0; i <= Nx + 1; i++)
    {
        for (j = 0; j <= Ny + 1; j++)
        {
            rou[i][j] = Q[i][j][0];
            u[i][j] = Q[i][j][1] / Q[i][j][0];
            v[i][j] = Q[i][j][2] / Q[i][j][0];
            p[i][j] = (gamma - 1) * (Q[i][j][3] - (Q[i][j][1] * Q[i][j][1] + Q[i][j][2] * Q[i][j][2]) / (2 * Q[i][j][0]));
            c[i][j] = sqrt(gamma * p[i][j] / rou[i][j]);
            H[i][j] = c[i][j] * c[i][j] / (gamma - 1) + 0.5 * (u[i][j] * u[i][j] + v[i][j] * v[i][j]);
        }
    }
    for (i = 0; i <= Nx; i++)
    {
        for (j = 0; j <= Ny + 1; j++)//x方向roe平均
        {
            rou_ax[i][j] = sqrt(rou[i][j] * rou[i + 1][j]);
            u_ax[i][j] = (u[i][j] * sqrt(rou[i][j]) + u[i + 1][j] * sqrt(rou[i + 1][j])) / (sqrt(rou[i][j]) + sqrt(rou[i + 1][j]));
            v_ax[i][j] = (v[i][j] * sqrt(rou[i][j]) + v[i + 1][j] * sqrt(rou[i + 1][j])) / (sqrt(rou[i][j]) + sqrt(rou[i + 1][j]));
            H_ax[i][j] = (H[i][j] * sqrt(rou[i][j]) + H[i + 1][j] * sqrt(rou[i + 1][j])) / (sqrt(rou[i][j]) + sqrt(rou[i + 1][j]));
            c_ax[i][j] = sqrt((gamma - 1) * (H_ax[i][j] - 0.5 * (u_ax[i][j] * u_ax[i][j] + v_ax[i][j] * v_ax[i][j])));
            p_ax[i][j] = rou_ax[i][j] * c_ax[i][j] * c_ax[i][j] / gamma;
        }
    }
    for (i = 0; i <= Nx + 1; i++)//y方向roe平均
    {
        for (j = 0; j <= Ny; j++)
        {
            rou_ay[i][j] = sqrt(rou[i][j] * rou[i][j + 1]);
            u_ay[i][j] = (u[i][j] * sqrt(rou[i][j]) + u[i][j + 1] * sqrt(rou[i][j + 1])) / (sqrt(rou[i][j]) + sqrt(rou[i][j + 1]));
            v_ay[i][j] = (v[i][j] * sqrt(rou[i][j]) + v[i][j + 1] * sqrt(rou[i][j + 1])) / (sqrt(rou[i][j]) + sqrt(rou[i][j + 1]));
            H_ay[i][j] = (H[i][j] * sqrt(rou[i][j]) + H[i][j + 1] * sqrt(rou[i][j + 1])) / (sqrt(rou[i][j]) + sqrt(rou[i][j + 1]));
            c_ay[i][j] = sqrt((gamma - 1) * (H_ay[i][j] - 0.5 * (u_ay[i][j] * u_ay[i][j] + v_ay[i][j] * v_ay[i][j])));
            p_ay[i][j] = rou_ay[i][j] * c_ay[i][j] * c_ay[i][j] / gamma;
        }
    }
    for (i = 0; i <= Nx; i++)
    {
        for (j = 0; j <= Ny + 1; j++)
        {
            lamda_A[i][j][0] = u_ax[i][j];
            lamda_A[i][j][1] = u_ax[i][j];
            lamda_A[i][j][2] = u_ax[i][j] - c_ax[i][j];
            lamda_A[i][j][3] = u_ax[i][j] + c_ax[i][j];
            LA[i][j][0][0] = (-(u_ax[i][j] * u_ax[i][j] + v_ax[i][j] * v_ax[i][j]) + 2 * c_ax[i][j] * c_ax[i][j] / (gamma - 1)) * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][0][1] = 2 * u_ax[i][j] * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][0][2] = 2 * v_ax[i][j] * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][0][3] = -2 * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][1][0] = -v_ax[i][j] * (u_ax[i][j] * u_ax[i][j] + v_ax[i][j] * v_ax[i][j]) * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][1][1] = 2 * u_ax[i][j] * v_ax[i][j] * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][1][2] = (2 * c_ax[i][j] * c_ax[i][j] / (gamma - 1) + 2 * v_ax[i][j] * v_ax[i][j]) * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][1][3] = -2 * v_ax[i][j] * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][2][0] = (0.5 * (u_ax[i][j] * u_ax[i][j] + v_ax[i][j] * v_ax[i][j]) + u_ax[i][j] * c_ax[i][j] / (gamma - 1)) * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][2][1] = (-u_ax[i][j] - c_ax[i][j] / (gamma - 1)) * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][2][2] = -v_ax[i][j] * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][2][3] = (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][3][0] = (0.5 * (u_ax[i][j] * u_ax[i][j] + v_ax[i][j] * v_ax[i][j]) - u_ax[i][j] * c_ax[i][j] / (gamma - 1)) * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][3][1] = (-u_ax[i][j] + c_ax[i][j] / (gamma - 1)) * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][3][2] = -v_ax[i][j] * (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            LA[i][j][3][3] = (gamma - 1) / (2 * c_ax[i][j] * c_ax[i][j]);
            RA[i][j][0][0] = 1;
            RA[i][j][0][1] = 0;
            RA[i][j][0][2] = 1;
            RA[i][j][0][3] = 1;
            RA[i][j][1][0] = u_ax[i][j];
            RA[i][j][1][1] = 0;
            RA[i][j][1][2] = u_ax[i][j] - c_ax[i][j];
            RA[i][j][1][3] = u_ax[i][j] + c_ax[i][j];
            RA[i][j][2][0] = 0;
            RA[i][j][2][1] = 1;
            RA[i][j][2][2] = v_ax[i][j];
            RA[i][j][2][3] = v_ax[i][j];
            RA[i][j][3][0] = 0.5 * (u_ax[i][j] * u_ax[i][j] - v_ax[i][j] * v_ax[i][j]);
            RA[i][j][3][1] = v_ax[i][j];
            RA[i][j][3][2] = H_ax[i][j] - c_ax[i][j] * u_ax[i][j];
            RA[i][j][3][3] = H_ax[i][j] + c_ax[i][j] * u_ax[i][j];
        }
    }
    for (i = 0; i <= Nx + 1; i++)
    {
        for (j = 0; j <= Ny; j++)
        {
            lamda_B[i][j][0] = v_ay[i][j];
            lamda_B[i][j][1] = v_ay[i][j];
            lamda_B[i][j][2] = v_ay[i][j] - c_ay[i][j];
            lamda_B[i][j][3] = v_ay[i][j] + c_ay[i][j];
            LB[i][j][0][0] = -u_ay[i][j] * (u_ay[i][j] * u_ay[i][j] + v_ay[i][j] * v_ay[i][j]) * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][0][1] = (2 * c_ay[i][j] * c_ay[i][j] / (gamma - 1) + 2 * u_ay[i][j] * u_ay[i][j]) * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][0][2] = 2 * u_ay[i][j] * v_ay[i][j] * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][0][3] = -2 * u_ay[i][j] * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][1][0] = (-(u_ay[i][j] * u_ay[i][j] + v_ay[i][j] * v_ay[i][j]) + 2 * c_ay[i][j] * c_ay[i][j] / (gamma - 1)) * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][1][1] = 2 * u_ay[i][j] * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][1][2] = 2 * v_ay[i][j] * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][1][3] = -2 * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][2][0] = (0.5 * (u_ay[i][j] * u_ay[i][j] + v_ay[i][j] * v_ay[i][j]) + v_ay[i][j] * c_ay[i][j] / (gamma - 1)) * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][2][1] = -u_ay[i][j] * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][2][2] = (-v_ay[i][j] - c_ay[i][j] / (gamma - 1)) * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][2][3] = (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][3][0] = (0.5 * (u_ay[i][j] * u_ay[i][j] + v_ay[i][j] * v_ay[i][j]) - v_ay[i][j] * c_ay[i][j] / (gamma - 1)) * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][3][1] = -u_ay[i][j] * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][3][2] = (-v_ay[i][j] + c_ay[i][j] / (gamma - 1)) * (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            LB[i][j][3][3] = (gamma - 1) / (2 * c_ay[i][j] * c_ay[i][j]);
            RB[i][j][0][0] = 0;
            RB[i][j][0][1] = 1;
            RB[i][j][0][2] = 1;
            RB[i][j][0][3] = 1;
            RB[i][j][1][0] = 1;
            RB[i][j][1][1] = 0;
            RB[i][j][1][2] = u_ay[i][j];
            RB[i][j][1][3] = u_ay[i][j];
            RB[i][j][2][0] = 0;
            RB[i][j][2][1] = v_ay[i][j];
            RB[i][j][2][2] = v_ay[i][j] - c_ay[i][j];
            RB[i][j][2][3] = v_ay[i][j] + c_ay[i][j];
            RB[i][j][3][0] = u_ay[i][j];
            RB[i][j][3][1] = -0.5 * (u_ay[i][j] * u_ay[i][j] - v_ay[i][j] * v_ay[i][j]);
            RB[i][j][3][2] = H_ay[i][j] - c_ay[i][j] * v_ay[i][j];
            RB[i][j][3][3] = H_ay[i][j] + c_ay[i][j] * v_ay[i][j];
        }
    }
    for (i = 0; i <= Nx + 1; i++)
    {
        for (j = 0; j <= Ny; j++)
        {
            for (k = 0; k <= 3; k++)
            {
                sum = 0;
                for (m = 0; m <= 3; m++)
                    sum += LB[i][j][k][m] * (Q[i][j + 1][m] - Q[i][j][m]);
                IBQ[i][j][k] = sum;
            }
        }
    }
    for (i = 0; i <= Nx + 1; i++)
    {
        for (j = 0; j <= Ny; j++)
        {
            for (k = 0; k <= 3; k++)
                g_aB[i][j][k] = 0.5 * (Qk(cy * lamda_B[i][j][k]) -
                    cy * lamda_B[i][j][k] * cy * lamda_B[i][j][k]) * IBQ[i][j][k];
        }
    }
    for (i = 0; i <= Nx + 1; i++)
    {
        for (j = 1; j <= Ny + 1; j++)
        {
            for (k = 0; k <= 3; k++)
                g_B[i][j][k] = minmod(g_aB[i][j][k], g_aB[i][j - 1][k]);
        }
    }
    for (i = 0; i <= Nx + 1; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            for (k = 0; k <= 3; k++)
            {
                if (IBQ[i][j][k] == 0)
                    GAMMA_B[i][j][k] = 0;
                else
                    GAMMA_B[i][j][k] = (g_B[i][j + 1][k] - g_B[i][j][k]) / IBQ[i][j][k];
            }
        }
    }
    for (i = 0; i <= Nx; i++)
    {
        for (j = 0; j <= Ny + 1; j++)
        {
            for (k = 0; k <= 3; k++)
            {
                sum = 0;
                for (m = 0; m <= 3; m++)
                    sum += LA[i][j][k][m] * (Q[i + 1][j][m] - Q[i][j][m]);
                IAQ[i][j][k] = sum;
            }
        }
    }
    for (i = 0; i <= Nx; i++)
    {
        for (j = 0; j <= Ny + 1; j++)
        {
            for (k = 0; k <= 3; k++)
                g_aA[i][j][k] = 0.5 * (Qk(cy * lamda_A[i][j][k]) - cy * lamda_A[i][j][k] * cy * lamda_A[i][j][k]) * IAQ[i][j][k];
        }
    }
    for (i = 1; i <= Nx + 1; i++)
    {
        for (j = 0; j <= Ny + 1; j++)
        {
            for (k = 0; k <= 3; k++)
                g_A[i][j][k] = minmod(g_aA[i][j][k], g_aA[i - 1][j][k]);
        }
    }
    for (i = 1; i <= Nx; i++)
    {
        for (j = 0; j <= Ny + 1; j++)
        {
            for (k = 0; k <= 3; k++)
            {
                if (IAQ[i][j][k] == 0)
                    GAMMA_A[i][j][k] = 0;
                else
                    GAMMA_A[i][j][k] = (g_A[i + 1][j][k] - g_A[i][j][k]) / IAQ[i][j][k];
            }
        }
    }
    for (i = 1; i <= Nx; i++)
    {
        for (j = 0; j <= Ny + 1; j++)
        {
            for (k = 0; k <= 3; k++)
                psai_A[i][j][k] = ((g_A[i][j][k] + g_A[i + 1][j][k] -
                    Qk(cx * lamda_A[i][j][k] + GAMMA_A[i][j][k]) * IAQ[i][j][k])) / cx;
        }
    }
    for (i = 0; i <= Nx + 1; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            for (k = 0; k <= 3; k++)
                psai_B[i][j][k] = ((g_B[i][j][k] + g_B[i][j + 1][k] -
                    Qk(cy * lamda_B[i][j][k] + GAMMA_B[i][j][k]) * IBQ[i][j][k])) / cy;
        }
    }
    for (i = 0; i <= Nx + 1; i++)
        for (j = 0; j <= Ny + 1; j++)
        {
            F[i][j][0] = rou[i][j] * u[i][j];
            F[i][j][1] = rou[i][j] * u[i][j] * u[i][j] + p[i][j];
            F[i][j][2] = rou[i][j] * u[i][j] * v[i][j];
            F[i][j][3] = (p[i][j] / (gamma - 1) + 0.5 * rou[i][j] * (u[i][j] * u[i][j] + v[i][j] * v[i][j]) + p[i][j]) * u[i][j];
        }
    for (i = 1; i <= Nx; i++)
    {
        for (j = 0; j <= Ny; j++)
        {
            for (k = 0; k <= 3; k++)
            {
                sum = 0;
                for (m = 0; m <= 3; m++)
                    sum = sum + psai_A[i][j][m] * RA[i][j][k][m];
                NF[i][j][k] = 0.5 * (F[i][j][k] + F[i + 1][j][k] + sum); //向量相加，对应元素相加
            }
        }
    }
    for (i = 0; i <= Nx + 1; i++)
        for (j = 0; j <= Ny + 1; j++)
        {
            G[i][j][0] = rou[i][j] * v[i][j];
            G[i][j][1] = rou[i][j] * u[i][j] * v[i][j];
            G[i][j][2] = rou[i][j] * v[i][j] * v[i][j] + p[i][j];
            G[i][j][3] = (p[i][j] / (gamma - 1) + 0.5 * rou[i][j] * (u[i][j] * u[i][j] + v[i][j] * v[i][j]) + p[i][j]) * v[i][j];
        }
    for (i = 0; i <= Nx; i++)
    {
        for (j = 1; j <= Ny; j++)
        {
            for (k = 0; k <= 3; k++)
            {
                sum = 0;
                for (m = 0; m <= 3; m++)
                    sum = sum + psai_B[i][j][m] * RB[i][j][k][m];
                NG[i][j][k] = 0.5 * (G[i][j][k] + G[i][j + 1][k] + sum); //向量相加，对应元素相加
            }
        }
    }
    for (i = 2; i <= Nx; i++)
    {
        for (j = 2; j <= Ny; j++)
        {
            for (k = 0; k <= 3; k++)
                Q[i][j][k] = Q[i][j][k] - cx * (NF[i][j][k] - NF[i - 1][j][k]) - cy * (NG[i][j][k] - NG[i][j - 1][k]);
        }
    }
}
/*-------------------------------------------------------
输出数据,写入文件
---------------------------------------------------------*/
void Output(double Q[Nx + 2][Ny + 2][4], double dx, double dy)
{
    int i, j;
    FILE* fp;
    double rou, u, v, p;
    fp = fopen("result.plt", "w");
    fprintf(fp, "TITLE     = \"Dataset\"\nVARIABLES = \"x\" \"y\" \"rou\" \"u\" \"v\" \"p\" \"E\"\n");
    fprintf(fp, "ZONE I=%d,J=%d,K=%d,F=POINT\n", Ny + 1, Nx + 1, 1);
    for (j = 0; j <= 0.85 / dy + 1; j++)
        for (i = j + (float)(1.15 / dx + 1); i <= Nx + 1; i++)
        {
            Q[i][j][0] = 1e-10;
            Q[i][j][1] = 1e-100;
            Q[i][j][2] = 1e-100;
            Q[i][j][3] = 1e-10;
        }
    for (i = 1; i <= Nx + 1; i++)
        for (j = 1; j <= Ny + 1; j++)
        {
            rou = Q[i][j][0];
            u = Q[i][j][1] / rou;
            v = Q[i][j][2] / rou;
            p = (gamma - 1) * (Q[i][j][3] - 0.5 * Q[i][j][0] * (u * u + v * v));
            fprintf(fp, "%20f%20f%20.10e%20.10e%20.10e%20.10e%20.10e\n", (i - 1) * dx, (j - 1) * dy, rou, u, v, p, Q[i][j][3]);
        }
    fclose(fp);

    fp = fopen("result.txt", "w");
    j = Ny / 2;
    for (i = 0; i <= Nx + 1; i++)
    {
        rou = Q[i][j][0];
        u = Q[i][j][1] / rou;
        v = Q[i][j][2] / rou;
        p = (gamma - 1) * (Q[i][j][3] - 0.5 * Q[i][j][0] * (u * u + v * v));
        fprintf(fp, "%20f%20.10e%20.10e%20.10e%20.10e%20.10e\n", i * dx, rou, u, v, p, Q[i][j][3]);
    }
    fclose(fp);
}

/*-----------------------------------------
主函数
------------------------------------------*/
int main()
{
    double T, dx, dy, dt;
    int i, j, k;
    Init(Q, dx, dy);
    T = 0;
    while (T < TT)
    {
        dt = CFL(Q, dx, dy);
        T += dt;
        bound(Q, dx, dy);
        TVD(Q, dx, dy, dt);
        printf("T=%10g    dt=%10g\n", T, dt);
    }
    Output(Q, dx, dy);
}
