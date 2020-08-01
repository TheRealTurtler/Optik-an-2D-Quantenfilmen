#define _CRT_SECURE_NO_WARNINGS

// ========================== Includes ==============================

#include <stdio.h>
#include <math.h>

// ======================== Deklarationen ===========================

double calculateD(double m_w, double m_b, double E_n, double V_0);
double* calculateHoleEnergy(double d_e, double V_0e, double E_nh, double m_wh, double m_bh, double stepSizeE_nh, double toleranz);
double* calculateHoleEnergy(double d_e, double V_0e, double E_nh, double m_wh, double m_bh, double toleranz);
void calcAndSave(FILE* file, double startE_ne, double endE_ne, double startE_nh, int steps, double stepSize, double toleranz);

// =========================== Macros ===============================

// Energie von eV zu J umrechnen
#define _EV(x) ((x) * eVtoJ)

// ====================== Phys. Konstanten ==========================

// reduziertes Planck'sches Wirkungsquantum
const double planckRed = 1.054571817e-34;

// Ruhemasse eines Elektrons
const double m_0 = 9.1093837015e-31;

// Umrechnungskonstante für eV zu J
const double eVtoJ = 1.602176e-19;

// ============================ Main ================================

int main()
{
	// Rechentoleranz
	const double toleranz = 1e-3;

	// Ausgabedatei öffnen
	FILE* file = fopen("Optik_an_2D_Quantenfilmen.dat", "w");

	// Spalten beschriften
	fprintf(file, "%s", "E_ne [eV]\t");
	fprintf(file, "%s", "d_e [nm]\t");
	fprintf(file, "%s", "E_nhh [eV]\t");
	fprintf(file, "%s", "d_hh [nm]\t");
	fprintf(file, "%s", "E_nlh [eV]\t");
	fprintf(file, "%s", "d_lh [nm]\n");

	// Werte berechnen und in Datei speichern
	calcAndSave(file, _EV(1e-3), _EV(1e-2), _EV(1e-5), 9, _EV(1e-7), toleranz);
	calcAndSave(file, _EV(1e-2), _EV(0.110), _EV(1e-5), 100, _EV(1e-6), toleranz);
	calcAndSave(file, _EV(0.110), _EV(0.160), _EV(1e-5), 25, _EV(1e-5), toleranz);

	// Datei schließen
	fclose(file);
}

// ========================= Funktionen =============================

// Dicke des Potenzialtops berechnen
double calculateD(double m_w, double m_b, double E_n, double V_0)
{
	return ((2.0 * planckRed) / (sqrt(2.0 * m_w * E_n)) * atan(sqrt((m_w * (V_0 - E_n)) / (m_b * E_n))));
}

// Energie der Löcher berechnen
double* calculateHoleEnergy(double d_e, double V_0e, double E_nh, double m_wh, double m_bh, double stepSizeE_nh, double toleranz)
{
	const double V_0h = 0.140 * eVtoJ;

	static double returnArray[2];

	double d_h = calculateD(m_wh, m_bh, E_nh, V_0h);

	double d_h_prev;

	while (!(fabs(d_e - d_h) < toleranz * d_e))
	{
		// Schrittweite
		E_nh += stepSizeE_nh;

		// Abbruchbedingung
		if (E_nh >= V_0e)
		{
			break;
		}

		// Vorherigen Wert für Breite speichern
		d_h_prev = d_h;

		// d für Löcher erneut berechnen
		d_h = calculateD(m_wh, m_bh, E_nh, V_0h);

		// Abbruchbedingung
		if ((d_h < d_e) && (d_h < d_h_prev))
		{
			printf("Toleranz nicht erreicht, Werte zu klein!\n");
			break;
		}
		else if ((d_h > d_e) && (d_h > d_h_prev))
		{
			printf("Toleranz nicht erreicht, Werte zu groß!\n");
			break;
		}
	}

	returnArray[0] = E_nh;
	returnArray[1] = d_h;

	return returnArray;
}

// Light und Heavy Hole Energie, sowie Dicke des Potenzialtopfs berechnen und in Datei speichern
void calcAndSave(FILE* file, double startE_ne, double endE_ne, double startE_nh, int stepsE_ne, double stepSizeE_nh, double toleranz)
{
	// Konstanten
	const double m_we = 0.067 * m_0;
	const double m_be = 0.09388 * m_0;
	const double m_wlh = 0.087 * m_0;
	const double m_blh = 0.10716 * m_0;
	const double m_whh = 0.34 * m_0;
	const double m_bhh = 0.396 * m_0;
	const double V_0e = 0.259 * eVtoJ;

	double E_ne = startE_ne;

	double* array;

	// Berechnung und Ausgabe in Datei
	for (E_ne; E_ne <= endE_ne; E_ne += (endE_ne - startE_ne) / stepsE_ne)
	{
		// Breite für Elektronen berechnen
		double d_e = calculateD(m_we, m_be, E_ne, V_0e);

		fprintf(file, "%lf\t", E_ne * 1 / eVtoJ);
		fprintf(file, "%lf\t", d_e * 1 / 1e-9);

		// Breite für Heavy Holes berechnen
		array = calculateHoleEnergy(d_e, V_0e, startE_nh, m_whh, m_bhh, stepSizeE_nh, toleranz);

		fprintf(file, "%lf\t", array[0] * 1 / eVtoJ);
		fprintf(file, "%lf\t", array[1] * 1 / 1e-9);

		// Breite für Light Holes berechnen
		array = calculateHoleEnergy(d_e, V_0e, startE_nh, m_wlh, m_blh, stepSizeE_nh, toleranz);

		fprintf(file, "%lf\t", array[0] * 1 / eVtoJ);
		fprintf(file, "%lf\n", array[1] * 1 / 1e-9);
	}
}