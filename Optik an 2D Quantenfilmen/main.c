#include "main.h"

// ================ Phys. Konstanten ======================
const double planckRed = 1.054571817e-34;
const double m_0 = 9.1093837015e-31;
const double eVtoJ = 1.602176e-19;

// ================ Main ==================================

int main()
{
	// Konstanten
	const double toleranz = 1e-3;
		
	const double m_we = 0.067 * m_0;
	const double m_be = 0.09388 * m_0;
	const double m_wlh = 0.087 * m_0;
	const double m_blh = 0.10716 * m_0;
	const double m_whh = 0.34 * m_0;
	const double m_bhh = 0.396 * m_0;

	const double V_0e = 0.259 * eVtoJ;

	// Startwert für Elektronenenergie
	double E_ne = 1e-2 * eVtoJ;

	// Startwert für Lochenergie
	double E_nh = 1e-5 * eVtoJ;

	double* array;


	FILE* file = fopen("Optik_an_2D_Quantenfilmen.dat", "w");

	fprintf(file, "%s", "E_ne [eV]\t");
	fprintf(file, "%s", "d_e [nm]\t");
	fprintf(file, "%s", "E_nhh [eV]\t");
	fprintf(file, "%s", "d_hh [nm]\t");
	fprintf(file, "%s", "E_nlh [eV]\t");
	fprintf(file, "%s", "d_lh [nm]\n");

	for (int i = 0; i < 100; i++)
	{
		// Breite für Elektronen berechnen
		double d_e = calculateD(m_we, m_be, E_ne, V_0e);

		fprintf(file, "%lf\t", E_ne * 1 / eVtoJ);
		fprintf(file, "%lf\t", d_e * 1 / 1e-9);

		// Breite für Heavy Holes berechnen
		array = calculateHoleEnergy(d_e, V_0e, E_nh, m_whh, m_bhh, toleranz);

		fprintf(file, "%lf\t", array[0]* 1 / eVtoJ);
		fprintf(file, "%lf\t", array[1] * 1 / 1e-9);

		// Breite für Light Holes berechnen
		array = calculateHoleEnergy(d_e, V_0e, E_nh, m_wlh, m_blh, toleranz);

		fprintf(file, "%lf\t", array[0] * 1 / eVtoJ);
		fprintf(file, "%lf\n", array[1] * 1 / 1e-9);

		E_ne += 1e-3 * eVtoJ;
	}
	//printf("\nErgebnis:\n");
	//printf("E_ne= %lf nm\n", E_ne * 1 / eVtoJ);
	//printf("d_e= %lf nm\n", d_e * 1 / 1e-9);

	

	//printf("E_nhh= %lf eV\n", array[0] * 1 / eVtoJ);
	//printf("d_hh= %lf nm\n", array[1] * 1 / 1e-9);

	
		
	//printf("E_nlh= %lf eV\n", array[0] * 1 / eVtoJ);
	//printf("d_lh= %lf nm\n", array[1] * 1 / 1e-9);

	//printf("Toleranz: %.3lf %%\n", toleranz * 100);

	// Tabelle
	// E_ne | d_e | E_nhh | d_hh | E_lh | d_lh

	fclose(file);
}

// ================ Funktionen ============================

double calculateD(double m_w, double m_b, double E_n, double V_0)
{
	return ((2.0 * planckRed) / (sqrt(2.0 * m_w * E_n)) * atan(sqrt((m_w * (V_0 - E_n)) / (m_b * E_n))));
}

double* calculateHoleEnergy(double d_e, double V_0e, double E_nh, double m_wh, double m_bh, double toleranz)
{
	const double V_0h = 0.140 * eVtoJ;

	static double returnArray[2];
		
	double d_h = calculateD(m_wh, m_bh, E_nh, V_0h);

	double d_h_prev;

	while (!(fabs(d_e - d_h) < toleranz * d_e))
	{
		//printf("%lf\t", E_nh * 1 / eVtoJ);
		//printf("%lf\t", d_e * 1 / 1e-9);
		//printf("%lf\n", d_h * 1 / 1e-9);

		// Schrittweite
		E_nh += 1e-6 * eVtoJ;

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

	//printf("\nErgebnis:\n");
	//printf("E_nh= %lf eV\n", E_nh * 1 / eVtoJ);
	//printf("d_e= %lf nm\n", d_e * 1 / 1e-9);
	//printf("d_h= %lf nm\n", d_h * 1 / 1e-9);
	//printf("Toleranz: %.3lf %%\n", toleranz * 100);

	return returnArray;
}