#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
#define Ms 2.0 //�����
#define gamma 1.4//���峣��
#define PI   3.141592654
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

#define Lx 2.0//��������
#define Ly 1.0

#define TT 0.1//��ʱ��
#define cfl 0.1 //CFL��

#define Nx 400//������
#define Ny 200
//ȫ�ֱ���
extern double Q[Nx + 2][Ny + 2][4], F[Nx + 2][Ny + 2][4], G[Nx + 2][Ny + 2][4];

double CFL(double dx, double dy);

double Qk(double x);

double minmod(double a, double b);

void Init(double& dx, double& dy);

void bound(double dx, double dy);

void TVD(double dx, double dy, double dt);

void Output(double dx, double dy);