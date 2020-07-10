#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>

double calculateD(double m_w, double m_b, double E_n, double V_0);

double* calculateHoleEnergy(double d_e, double V_0e, double E_nh, double m_wh, double m_bh, double stepSizeE_nh, double toleranz);

double* calculateHoleEnergy(double d_e, double V_0e, double E_nh, double m_wh, double m_bh, double toleranz);

void calcAndSave(FILE* file, double startE_ne, double endE_ne, double startE_nh, int steps, double stepSize, double toleranz);