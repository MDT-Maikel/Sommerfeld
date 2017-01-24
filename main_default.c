/*--
	micrOMEGAs for Colored Dark Sectors

	The code add Sommerfeld corrections for annihilation of colored particles
	as presented in arXiv:1612.02825. These colored particles are identified
	by their PDG code, see the code for	more details. This code works with a
	FeynRules model file shipped with the paper arXiv:1701.abcde.

	Also included are bound state formation effects as described in the paper
	arxiv:1611.08133 (extended to color sextets in the paper arXiv:1701.abcde).

	This code needs a modified version of micrOMEGAs v4.3.2 and a diff can be
	found in this file "micromegas_4.3.2_bound_states.patch".

	Run the code as
		./main <file with parameters> <sommerfeld> <bound state formation>
	for example
		./main data.par off on
	to run without Sommerfeld corrections but with bound state corrections.
--*/


#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"
#include "lib/pmodel.h"
#include "stdbool.h"


// Variable which sets sommerfeld corrections on or off.
static bool sommerfeld_on;

// Variable which sets bound state formation on or off.
static bool bsf_on;

// Helper functions.
long color(long pdg);
long spin(long pdg);
double casimir2(int color);
double invexp(double x);
double alpha_strong(double q);

// Cross section functions.
double xx_to_qq(double alpha_s, double alpha_sommerfeld, int rep, int spin, double m, double v, bool sommerfeld);
double xx_to_gg(double alpha_s, double alpha_sommerfeld, int rep, int spin, double m, double v, bool sommerfeld);

double ss_to_qq(double alpha_s, double alpha_sommerfeld, int rep, double m, double v);
double ff_to_qq(double alpha_s, double alpha_sommerfeld, int rep, double m, double v);
double vv_to_qq(double alpha_s, double alpha_sommerfeld, int rep, double m, double v);

double ss_to_gg(double alpha_s, double alpha_sommerfeld, int rep, double m, double v);
double ff_to_gg(double alpha_s, double alpha_sommerfeld, int rep, double m, double v);
double vv_to_gg(double alpha_s, double alpha_sommerfeld, int rep, double m, double v);

double ss_to_qq_sommerfeld(double alpha_s, double alpha_sommerfeld, int rep, double m, double v);
double ff_to_qq_sommerfeld(double alpha_s, double alpha_sommerfeld, int rep, double m, double v);
double vv_to_qq_sommerfeld(double alpha_s, double alpha_sommerfeld, int rep, double m, double v);

double ss_to_gg_sommerfeld(double alpha_s, double alpha_sommerfeld, int rep, double m, double v);
double ff_to_gg_sommerfeld(double alpha_s, double alpha_sommerfeld, int rep, double m, double v);
double vv_to_gg_sommerfeld(double alpha_s, double alpha_sommerfeld, int rep, double m, double v);

// Bound state formation functions.
double *alpha_table;
double exp_cut(double x);
void read_table_alpha(void);
double alphaS_bs(int color, double m, double *alpha_table);
static void r_simpson(double(*func)(double, Parameters), Parameters pars, double * f, double a, double b, double eps, double *aEps, double *ans, double *aAns, int *deepness);
double simpsonArg( double (*func)(double, Parameters), Parameters pars, double a,double b, double  eps);
double zeta(int color, double m);
double zetap(int color, double m);
double kappa(int color, double m);
double BE(int color, double m);
double nu(int color, double m, double z, double u);
double bohr_radius(int color, double m);
int g_freedom(int spin, int color);
double fMB(double m, double T, double vrel);
double sigmaDiss(int spin,int color, double m, double T, double u);
double GammaBS(int spin, int color, double m, int spin_eta);
double GammaDissIntegrand(double u, Parameters pars);
double GammaDiss(int spin, int color, double m, double T);
double sigmaBSFaveraged(int spin, int color, double m, double T);
double GammaBSaveraged(int spin, int color, double m, int spin_eta, double T);
double bound_state_rate(int spin, int color, double m, double T);


/*-- Main Program --*/

// Main part of the program.
int main(int argc, char** argv)
{
	int err;
	char cdmName[10];
	int spin2, charge3, cdim;
	ForceUG = 0;  /* to Force Unitary Gauge assign 1 */

	if (argc == 1)
	{
		printf("Correct usage: ./main <file with parameters> <sommerfeld> <bound state formation>\n");
		printf("Example: ./main data.par\n");
		exit(1);
	}

	// Determine if sommerfeld corrections and bound state formation are enabled.
	sommerfeld_on = argc >= 3 && strcmp(argv[2], "off") != 0;
	bsf_on = argc >= 4 && strcmp(argv[3], "off") != 0;
	printf("Sommerfeld corrections enabled: %s\n", sommerfeld_on ? "true" : "false");
	printf("Bound state formation enabled: %s\n", bsf_on ? "true" : "false");

	// Read in parameter file.
	err = readVar(argv[1]);
	if (err == -1)
	{
		printf("Can not open the file\n");
		exit(1);
	}
	else if (err > 0)
	{
		printf("Wrong file contents at line %d\n", err);
		exit(1);
	}

	err = sortOddParticles(cdmName);
	if (err)
	{
		printf("Can't calculate %s\n", cdmName);
		return 1;
	}

	if (CDM1) 
	{ 
		qNumbers(CDM1, &spin2, &charge3, &cdim);
		printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n", CDM1, spin2, Mcdm1); 
		if (charge3)
			printf("Dark Matter has electric charge %d/3\n", charge3);
		if (cdim != 1)
			printf("Dark Matter is a color particle\n");
	}
	if (CDM2) 
	{ 
		qNumbers(CDM2, &spin2, &charge3, &cdim);
		printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n", CDM2, spin2, Mcdm2); 
		if (charge3)
			printf("Dark Matter has electric charge %d/3\n", charge3);
		if (cdim != 1)
			printf("Dark Matter is a color particle\n");
	}
	printMasses(stdout, 1);

	// Read the table with alpha strong for bound states if needed.
	if (bsf_on)
		read_table_alpha();

	// Calculate the relic density.
	int fast = 0;
	double Beps = 1.E-7;
	double cut = 0.0001;
	double Omega, OmegaFO;
	printf("\n==== Calculation of relic density =====\n");   

	double Xf, XfFO;
	Omega = darkOmega(&Xf, fast, Beps);
	OmegaFO = darkOmegaFO(&XfFO, fast, Beps);

	printf("Xf=%.4e Omega=%.4e\n", Xf, Omega);
	printf("Xf(FO)=%.4e Omega(FO)=%.4e\n", XfFO, OmegaFO);
	printChannels(XfFO, cut, Beps, 1, stdout);
	printf("omega_h^2 = %.4E\n", Omega);
	printf("omega_h^2(FO) = %.4E\n", OmegaFO);

	killPlots();
	return 0;
}


/*-- Cross Section Improvement --*/

void improveCrossSection(long n1, long n2, long n3, long n4, double pin, double *res)
{
	// Return zero for all process which do not have two colored X's.
	if (abs(n1) < 9000000 || abs(n2) < 9000000 || color(n1) < 3 || color(n2) < 3)
	{
		printf("WARNING: process %d %d -> %d %d is being ignored\n", (int)n1, (int)n2, (int)n3, (int)n4);
		*res = 0;
		return;	
	}

	// Get incoming particle masses.
	double m1 = pMass(pdg2name(n1));
	double m2 = pMass(pdg2name(n2));
	double m = (m1 + m2) / 2.0;
	if (m1 != m2)
		printf("WARNING: masses of incoming particles are are not equal: %f and %f", m1, m2);

	// Calculate beta and s.
	double v = pin / sqrt(pow(pin, 2.0) + pow(m1, 2.0));
	double s = pow(sqrt(pow(pin, 2.0) + pow(m1, 2.0)) + sqrt(pow(pin, 2.0) + pow(m2, 2.0)), 2.0);

	// Calculate alpha_sommerfeld at the scale of the soft gluons.
	double alpha_sommerfeld = alpha_strong(pin);
	// micrOMEGAs uses its own running for the hard process.
	double alpha_mo = parton_alpha(GGscale);

	// Determine color and spin of X.
	long color_x = color(n1);
	long spin_x = spin(n1);
	if (color_x != 3 && color_x != 6 && color_x != 8)
		printf("color of x is invalid: %d", (int)color_x);
	if (spin_x < 1 || spin_x > 6)
		printf("spin of x is invalid: %d", (int)spin_x);

	// Add sommerfeld factor for XX -> qq.
	if (abs(n1) == abs(n2) && abs(n3) >= 1 && abs(n3) <= 6 && abs(n4) >= 1 && abs(n4) <= 6)
	{
		double xsec_mo = *res;
		double xsec = xx_to_qq(alpha_mo, alpha_sommerfeld, color_x, spin_x, m, v, sommerfeld_on);
		// Safety check: xsec is not a number.
		if (!isfinite(xsec) || isnan(xsec))
		{
			printf("WARNING: xsec not a number (%e) for %d %d -> %d %d\n", xsec, (int)n1, (int)n2, (int)n3, (int)n4);
			printf("\tmass: %f, s: %f, v: %f, p: %f, alpha_s: %f, alpha_sommerfeld: %f\n", m, sqrt(s), v, pin, alpha_mo, alpha_sommerfeld);
		}
		// Safety check: micromegas vs. analytic cross section (0.1% agreement needed). 
		if (!sommerfeld_on && fabs(1000 * (xsec - xsec_mo) / xsec_mo) > 1)
		{
			printf("WARNING: xsec mismatch to analytic for %d %d -> %d %d\n", (int)n1, (int)n2, (int)n3, (int)n4);
			printf("\tmass: %f, s: %f, v: %f, p: %f, alpha_s: %f, alpha_sommerfeld: %f\n", m, sqrt(s), v, pin, alpha_mo, alpha_sommerfeld);
			printf("\txsec(mo): %.8e, xsec(analytic): %.8e, ratio(mo/analytic): %.6f\n", xsec_mo, xsec, xsec_mo / xsec);
		}
		*res = xsec;
		return;
	}

	// Add sommerfeld factor for XX -> gg.
	if (abs(n1) == abs(n2) && n3 == 21 && n4 == 21)
	{
		double xsec_mo = *res;
		double xsec = xx_to_gg(alpha_mo, alpha_sommerfeld, color_x, spin_x, m, v, sommerfeld_on);
		// Safety check: xsec is not a number.
		if (!isfinite(xsec) || isnan(xsec))
		{
			printf("WARNING: xsec not a number (%e) for %d %d -> %d %d\n", xsec, (int)n1, (int)n2, (int)n3, (int)n4);
			printf("\tmass: %f, s: %f, v: %f, p: %f, alpha_s: %f, alpha_sommerfeld: %f\n", m, sqrt(s), v, pin, alpha_mo, alpha_sommerfeld);
		}
		// Safety check: micromegas vs. analytic cross section (0.1% agreement needed). 
		if (!sommerfeld_on && fabs(1000 * (xsec - xsec_mo) / xsec_mo) > 1)
		{
			printf("WARNING: xsec mismatch to analytic for %d %d -> %d %d\n", (int)n1, (int)n2, (int)n3, (int)n4);
			printf("\tmass: %f, s: %f, v: %f, p: %f, alpha_s: %f, alpha_sommerfeld: %f\n", m, sqrt(s), v, pin, alpha_mo, alpha_sommerfeld);
			printf("\txsec(mo): %.8e, xsec(analytic): %.8e, ratio(mo/analytic): %.6f\n", xsec_mo, xsec, xsec_mo / xsec);
		}
		*res = xsec;
		return;
	}

	// Add sommerfeld factor for X1 X2 -> (q/g) (q/g).
	if ((n3 == 21 || (abs(n3) >= 1 && abs(n3) <= 6)) && (n4 == 21 || (abs(n4) >= 1 && abs(n4) <= 6)))
	{
		if (sommerfeld_on)
			printf("WARNING: no Sommerfeld corrections for %d %d -> %d %d\n", (int)n1, (int)n2, (int)n3, (int)n4);
		return;
	}

	// We only allow the channels XX -> gg and XX -> qq, in other cases the cross section is zero.
	printf("in: %d %d, out: %d %d set to zero\n", (int)n1, (int)n2, (int)n3, (int)n4);
	*res = 0;
	return; 
}


/*-- Helpers --*/

long color(long pdg)
{
	// The color is coded in the last digit of the PDG number.
	return abs(pdg) % 100;
}

long spin(long pdg)
{
	// The spin is coded in the 3rd to last digit of the PDG number.
	return abs(pdg / 100) % 100;
}

double casimir2(int color)
{
    switch (color)
    {
		case 3: return 4./3;
		case 6: return 10./3;
		case 8: return 3;
		default:
			printf("WARNING: casimir2 called with color: %d.\n", color);
			exit(42);
    }
}

double invexp(double x)
{
	// This function is defined to prevent numerical issues for low x.	
	return exp(x) / (exp(x) - 1);
}

// See the Mathematica notebook for the details of this definition of alpha_strong.
double alpha_strong(double q)
{
	// Implement a cut off for the momentum q of 1 GeV to not enter the non-perturbative regime.
	q = fmax(q, 1.0);
	// MSbar masses for the quarks.
	double mtop = 160;
	double mbottom = 4.18;
	double mcharm = 1.28;
   	// Determine number of active flavors.
	double nf = 6.0;	
	if (q < mcharm) nf = 3.0;
	else if (q < mbottom) nf = 4.0;
	else if (q < mtop) nf = 5.0;
	// Determine the threshold lambdas.
	double lambda = 0.08896768177299201;
	if (q < mcharm) lambda = 0.33348050663724466;
	else if (q < mbottom) lambda = 0.2913885366061117;
	else if (q < mtop) lambda = 0.20953346238097081;
	// Determine alpha_strong and return it.
	double t = log(pow(q / lambda, 2.0));
	double z3 = 1.202056903159594;
	double b0 = (33.0 - 2.0 * nf) / (12.0 * M_PI);
	double b1 = (153.0 - 19.0 * nf) / (24.0 * pow(M_PI, 2.0));
	double b2 = (2857.0 - 5033.0 / 9.0 * nf + 325.0 / 27.0 * pow(nf, 2.0)) / (128.0 * pow(M_PI, 3.0));
	double b3 = ((149753.0 / 6.0 + 3564.0 * z3) - (1078361.0 / 162.0 + 6508.0 / 27.0 * z3) * nf + (50065.0 / 162.0 + 6472.0 / 81.0 * z3) * pow(nf, 2.0) + 1093.0 / 729.0 * pow(nf, 3.0)) / (256.0 * pow(M_PI, 4.0));
	return 1.0 / (b0 * t) * (1.0 - b1 / pow(b0, 2.0) * log(t) / t + (pow(b1, 2.0) * (pow(log(t), 2.0) - log(t) - 1.0) + b0 * b2) / (pow(b0, 4.0) * pow(t, 2.0)) - 1.0 / (pow(b0, 6.0) * pow(t, 3.0)) * (pow(b1, 3.0) * (pow(log(t), 3.0) - 5.0 / 2.0 * pow(log(t), 2.0) - 2.0 * log(t) + 1.0 / 2.0) + 3.0 * b0 * b1 * b2 * log(t) - 0.5 * pow(b0, 2.0) * b3));
}


/*-- Cross Sections --*/

double xx_to_qq(double alpha_s, double alpha_sommerfeld, int rep, int spin, double m, double v, bool sommerfeld)
{
	if (!sommerfeld)
	{
		switch (spin)
		{
			case 1: case 2: return ss_to_qq(alpha_s, alpha_sommerfeld, rep, m, v);
			case 3: case 4: return ff_to_qq(alpha_s, alpha_sommerfeld, rep, m, v);
			case 5: case 6: return vv_to_qq(alpha_s, alpha_sommerfeld, rep, m, v);
			default: printf("WARNING: xx_to_qq called for invalid spin %d.\n", spin); return 0.0;
		}
	}
	else
	{
		switch (spin)
		{
			case 1: case 2: return ss_to_qq_sommerfeld(alpha_s, alpha_sommerfeld, rep, m, v);
			case 3: case 4: return ff_to_qq_sommerfeld(alpha_s, alpha_sommerfeld, rep, m, v);
			case 5: case 6: return vv_to_qq_sommerfeld(alpha_s, alpha_sommerfeld, rep, m, v);
			default: printf("WARNING: xx_to_qq called for invalid spin %d.\n", spin); return 0.0;
		}
	}
	return 0.0;
}

double xx_to_gg(double alpha_s, double alpha_sommerfeld, int rep, int spin, double m, double v, bool sommerfeld)
{
	if (!sommerfeld)
	{
		switch (spin)
		{
			case 1: case 2: return ss_to_gg(alpha_s, alpha_sommerfeld, rep, m, v);
			case 3: case 4: return ff_to_gg(alpha_s, alpha_sommerfeld, rep, m, v);
			case 5: case 6: return vv_to_gg(alpha_s, alpha_sommerfeld, rep, m, v);
			default: printf("WARNING: gg_to_qq called for invalid spin %d.\n", spin); return 0.0;
		}
	}
	else
	{
		switch (spin)
		{
			case 1: case 2: return ss_to_gg_sommerfeld(alpha_s, alpha_sommerfeld, rep, m, v);
			case 3: case 4: return ff_to_gg_sommerfeld(alpha_s, alpha_sommerfeld, rep, m, v);
			case 5: case 6: return vv_to_gg_sommerfeld(alpha_s, alpha_sommerfeld, rep, m, v);
			default: printf("WARNING: gg_to_qq called for invalid spin %d.\n", spin); return 0.0;
		}
	}
	return 0.0;
}

double ss_to_qq(double alpha_s, double alpha_sommerfeld, int rep, double m, double v)
{
	switch (rep)
	{
		case 3: return ID_S3S3_QQ_NS;
		case 6: return ID_S6S6_QQ_NS;
		case 8: return ID_S8S8_QQ_NS;
		default: printf("WARNING: ss_to_qq called for invalid representation %d.\n", rep); return 0.0;
	}
	return 0.0;
}

double ff_to_qq(double alpha_s, double alpha_sommerfeld, int rep, double m, double v)
{
	switch (rep)
	{
		case 3: return ID_F3F3_QQ_NS;
		case 6: return ID_F6F6_QQ_NS;
		case 8: return ID_F8F8_QQ_NS;
		default: printf("WARNING: ff_to_qq called for invalid representation %d.\n", rep); return 0.0;
	}
	return 0.0;
}

double vv_to_qq(double alpha_s, double alpha_sommerfeld, int rep, double m, double v)
{
	switch (rep)
	{
		case 3: return ID_V3V3_QQ_NS;
		case 6: return ID_V6V6_QQ_NS;
		case 8: return ID_V8V8_QQ_NS;
		default: printf("WARNING: vv_to_qq called for invalid representation %d.\n", rep); return 0.0;
	}
	return 0.0;
}

double ss_to_gg(double alpha_s, double alpha_sommerfeld, int rep, double m, double v)
{
	switch (rep)
	{
		case 3: return ID_S3S3_GG_NS;
		case 6: return ID_S6S6_GG_NS;
		case 8: return ID_S8S8_GG_NS;
		default: printf("WARNING: vv_to_gg called for invalid representation %d.\n", rep); return 0.0;
	}
	return 0.0;
}

double ff_to_gg(double alpha_s, double alpha_sommerfeld, int rep, double m, double v)
{
	switch (rep)
	{
		case 3: return ID_F3F3_GG_NS;
		case 6: return ID_F6F6_GG_NS;
		case 8: return ID_F8F8_GG_NS;
		default: printf("WARNING: vv_to_gg called for invalid representation %d.\n", rep); return 0.0;
	}
	return 0.0;
}

double vv_to_gg(double alpha_s, double alpha_sommerfeld, int rep, double m, double v)
{
	switch (rep)
	{
		case 3: return ID_V3V3_GG_NS;
		case 6: return ID_V6V6_GG_NS;
		case 8: return ID_V8V8_GG_NS;
		default: printf("WARNING: vv_to_gg called for invalid representation %d.\n", rep); return 0.0;
	}
	return 0.0;
}

double ss_to_qq_sommerfeld(double alpha_s, double alpha_sommerfeld, int rep, double m, double v)
{
	switch (rep)
	{
		case 3: return ID_S3S3_QQ_SO;
		case 6: return ID_S6S6_QQ_SO;
		case 8: return ID_S8S8_QQ_SO;
		default: printf("WARNING: ss_to_qq_sommerfeld called for invalid representation %d.\n", rep); return 0.0;
	}
	return 0.0;
}

double ff_to_qq_sommerfeld(double alpha_s, double alpha_sommerfeld, int rep, double m, double v)
{
	switch (rep)
	{
		case 3: return ID_F3F3_QQ_SO;
		case 6: return ID_F6F6_QQ_SO;
		case 8: return ID_F8F8_QQ_SO;
		default: printf("WARNING: ff_to_qq_sommerfeld called for invalid representation %d.\n", rep); return 0.0;
	}
	return 0.0;
}

double vv_to_qq_sommerfeld(double alpha_s, double alpha_sommerfeld, int rep, double m, double v)
{
	switch (rep)
	{
		case 3: return ID_V3V3_QQ_SO;
		case 6: return ID_V6V6_QQ_SO;
		case 8: return ID_V8V8_QQ_SO;
		default: printf("WARNING: vv_to_qq_sommerfeld called for invalid representation %d.\n", rep); return 0.0;
	}
	return 0.0;
}

double ss_to_gg_sommerfeld(double alpha_s, double alpha_sommerfeld, int rep, double m, double v)
{
	switch (rep)
	{
		case 3: return ID_S3S3_GG_SO;
		case 6: return ID_S6S6_GG_SO;
		case 8: return ID_S8S8_GG_SO;
		default: printf("WARNING: ss_to_gg_sommerfeld called for invalid representation %d.\n", rep); return 0.0;
	}
	return 0.0;
}

double ff_to_gg_sommerfeld(double alpha_s, double alpha_sommerfeld, int rep, double m, double v)
{
	switch (rep)
	{
		case 3: return ID_F3F3_GG_SO;
		case 6: return ID_F6F6_GG_SO;
		case 8: return ID_F8F8_GG_SO;
		default: printf("WARNING: ff_to_gg_sommerfeld called for invalid representation %d.\n", rep); return 0.0;
	}
	return 0.0;
}

double vv_to_gg_sommerfeld(double alpha_s, double alpha_sommerfeld, int rep, double m, double v)
{
	switch (rep)
	{
		case 3: return ID_V3V3_GG_SO;
		case 6: return ID_V6V6_GG_SO;
		case 8: return ID_V8V8_GG_SO;
		default: printf("WARNING: vv_to_gg_sommerfeld called for invalid representation %d.\n", rep); return 0.0;
	}
	return 0.0;
}


/*-- Improve Averaged Cross Section --*/

double improveAveragedCrossSection(long n1, long n2, double mdm, double T)
{
	if (!bsf_on)
		return 0.0;

	// Return zero for all process which do not have two equal colored X's.
	if (abs(n1) < 9000000 || abs(n2) < 9000000 || abs(n1) != abs(n2))
	{
		printf("WARNING: bound state corrections for process %d %d -> ??? are being ignored\n", (int)n1, (int)n2);
		return 0.0;	
	}

	// Get incoming particle masses.
	double m1 = pMass(pdg2name(n1));
	double m2 = pMass(pdg2name(n2));
	double m = (m1 + m2) / 2.0;
	if (m1 != m2)
		printf("WARNING: masses of incoming particles are are not equal: %f and %f", m1, m2);

	// Determine color and spin of X.
	long color_x = color(n1);
	long spin_x = spin(n1);

	// Calculate the bound state formation rate.
	double bsf_rate = bound_state_rate(spin_x, color_x, m, T);
	// Safety check: bound state formation rate is not a number.
	if (!isfinite(bsf_rate) || isnan(bsf_rate))
	{
		printf("WARNING: averaged BSF rate is not a number (%e) for %d %d -> ??? \n", bsf_rate, (int)n1, (int)n2);
		printf("\tmass: %f, T: %f\n", m, T);
	}
	return bsf_rate;
}


/*-- Bound State Formation --*/

// Helper function to compute an exponential without underflowing.
double exp_cut(double x)
{
	if (x < -50)
		return 0;
	return exp(x);
}

void read_table_alpha(void)
{
	// The table that is being read has 2001 times 3 entries.
	int nr_entries = 2001;
	alpha_table = (double*) calloc(3 * nr_entries, sizeof(double));
	FILE * falpha = fopen("alpha_strong_bsf.txt", "r");
	int m;
	int counter = 0;
	while (!feof(falpha))
	{
		if (counter/3 > nr_entries)
		{
			printf("Size of alpha table greater than N=%d\n", nr_entries);
			exit(10);
		}
		fscanf(falpha, "%d %lf %lf %lf\n", &m, &(alpha_table[counter]), &(alpha_table[counter + 1]), &(alpha_table[counter + 2]));
		counter += 3;
	}
	if (falpha == NULL)
		perror("The folloriwng error occured");
	//fclose(falpha);
	return;
} 

// Alpha strong for the bound state formation is determined from a recursive formula in Mathermatica and is read in from a table.
double alphaS_bs(int color, double m, double *alpha_table)
{
	// Min, max and step size of the table.
	int mmin = 0, mmax = 20000, mstep = 10;
	int color_index;
	if (color == 3) color_index = 0;
	if (color == 6) color_index = 1;
	if (color == 8) color_index = 2;
	int mbin = (int) floor((m - mmin) / mstep + 0.001);
	return alpha_table[mbin * 3 + color_index];
}

// Integration method from micrOMEGAs.
static void r_simpson(double(*func)(double, Parameters), Parameters pars, double *f, double a, double b, double eps, double *aEps, double *ans, double *aAns, int *deepness)
{
	double f1[5];
	int i;
	int d1 = *deepness + 1, d2 = *deepness + 1;
	double s1, s2, s3, e_err;

	s1 = (f[0] + 4 * f[4] + f[8])/6;
	s2 = (f[0] + 4 * f[2] + 2 * f[4] + 4 * f[6] + f[8]) / 12;
	s3 = (f[0] + 4 * f[1] + 2 * f[2] + 4 * f[3] + 2 * f[4] + 4 * f[5] + 2 * f[6] + 4 * f[7] + f[8]) / 24;

	if (!isfinite(s3))
	{
		*ans = s3;
		*aAns = s3;
		return;
	}

	e_err = eps * fabs(s3);
	i = 0;
	if ((fabs(s3 - s2) < e_err && fabs(s3 - s1) < 16 * e_err))
		i = 1;
	else if (fabs(s3 - s2) * (b - a) < 0.1 * (*aEps) && fabs(s3 - s1) * (b - a) < 1.6 * (*aEps)) 
	{
		i = 1;
		*aEps -= fabs((s3 - s2) * (b - a));
	}

	if (i || *deepness > 20)
	{
		*ans += s3 * (b - a);
		*aAns += (fabs(f[0]) + 4 * fabs(f[2]) + 2 * fabs(f[4]) + 4 * fabs(f[6]) + fabs(f[8])) * fabs(b - a) / 12;
		return;
	}

	for (i = 0; i < 5; i++)
		f1[i] = f[4 + i];
	for (i = 8; i > 0; i -= 2)
		f[i] = f[i / 2];
	for (i = 1; i < 8; i += 2)
		f[i] = (*func)(a + i * (b - a) / 16, pars);

	r_simpson(func, pars, f, a, (a + b) / 2, eps, aEps, ans, aAns, &d1);
	for (i = 0; i < 5; i++)
		f[2 * i] = f1[i];
	for (i = 1; i < 8; i+= 2)
		f[i] = (*func)((a + b) / 2 + i * (b - a) / 16, pars);
	r_simpson(func, pars, f, (a + b) / 2, b, eps, aEps, ans, aAns, &d2);
	if (d1 > d2)
		*deepness = d1;
	else
		*deepness = d2;
	return;   
}

// Integration method from micrOMEGAs.
double simpsonArg(double (*func)(double, Parameters), Parameters pars, double a, double b, double eps)
{
	double f[9];
	double aEps;
	int i, j;	

	aEps=0;
	if (a==b)
		return 0;
	for (i = 0; i < 9; i++)
	{
		f[i] = (*func) (a + i * (b - a) / 8, pars);
		aEps += fabs(f[i]);
	}
	if (aEps==0.)
		return 0;
	eps = eps / 2;
	aEps = eps * aEps * fabs(b - a) / 9;

	for (j = 0; ; j++)
	{
		double ans=0.0, aAns=0.0; 
		int deepness = 1;
		r_simpson(func, pars, f, a, b, eps, &aEps, &ans, &aAns, &deepness);
		if (5 * aAns * eps > aEps || j >= 2 )
			return ans;
		if (!isfinite(aAns))
			return aAns;
		for (i = 0; i < 9; i++)
			f[i]=(*func)(a + i * (b - a) / 8, pars);
		aEps = aAns * eps;
	}
}

double zeta(int color, double m)
{
	double alphaS_boundstate = alphaS_bs(color, m, alpha_table);
	return casimir2(color) * alphaS_boundstate;
}

double zetap(int color, double m)
{
	double alphaS_boundstate = alphaS_bs(color, m, alpha_table);
	return (casimir2(color) - 1.5) * alphaS_boundstate;
}

double kappa(int color, double m)
{
	return zeta(color, m) / fabs(zetap(color, m));
}

double BE(int color, double m)
{
	double z = zeta(color, m);
	return pow(z, 2.0) * m / 4.0;
}

double nu(int color, double m, double z, double u)
{
	return 1.0 / kappa(color, m) * sqrt(z / u);
}

double bohr_radius(int color, double m)
{
	return 2.0 / (zeta(color, m) * m);
}

int g_freedom(int spin, int color)
{
	int spinfact = (spin - 1) / 2 + 1;
	return spinfact * color;
}

double fMB(double m, double T, double vrel)
{
	return pow(m / (4 * M_PI * T), 1.5) * 4 * M_PI * pow(vrel, 2.0) * exp_cut(-m * pow(vrel, 2.0) / (4.0 * T));
}

double sigmaDiss(int spin, int color, double m, double T, double u)
{
	double alphaS = parton_alpha(GGscale);
	double zetp = zetap(color, m);
	double k = kappa(color, m);
	double E = BE(color, m);
	double a = bohr_radius(color, m); 
	double z = E / T;
	double omega = E * (1 + u / z);
	// Compute prefactor with gluon averaging and color.
	int symmetry_factor = (spin - 1) % 2 == 0 ? 2: 1;
	double prefact = 1.0 / 8 * casimir2(color) * symmetry_factor;
	// Compute common factor (for attractive and repulsive potentials).
	double coeff = pow(2, 9) * pow(M_PI, 2.0) / 3.0 * alphaS * pow(a, 2.0) * pow(E / omega, 4.0);
	double coeff2 = 1.0;
	if (u > 1e-6)
	{
		double v = nu(color, m, z, u);
		coeff *= (1 + pow(v, 2))/(1 + pow(k * v, 2));
		coeff /= k * (1 - exp_cut(-2 * M_PI * v));
		// Compute factor that depends on whether the potential is attractive or repulsive.
		if (zetp > 0)
			coeff2 = exp_cut(-4 * v * atan(1.0 / (k * v)));	
		else
			coeff2 = exp_cut(4 * v * atan(1.0 / (k * v)) - 2.0 * M_PI * v);
	}
	else
	{
		// Indeterminate form in 0, needs to be treated separately.
		if (zetp > 0)
			coeff2 = exp_cut(-4.0 / k) / pow(k, 3.0);
		else
			coeff2 = 0;
	}
	return prefact * coeff * coeff2;
}

double sigmaStimulatedBSF(double u, Parameters pars)
{
	double E = BE(pars.color, pars.m);
	double z = E / pars.T;
	double vrel = zeta(pars.color, pars.m) * sqrt(u / z);
	double omega = E + 0.25 * pars.m * pow(vrel, 2.0);
	int gg = 16;
	int gX = g_freedom(pars.spin, pars.color);
	int spinfact = (pars.spin - 1) / 2 + 1;
	int symmetry_factor = (pars.spin - 1) % 2 == 0 ? 2 : 1;
	double stimulated_emission = 1 + 1/(exp(omega/pars.T) - 1);
	double prefact = gg * pow(omega, 2.0) / (pow(gX, 2.0) * pow(spinfact, 2.0) * pow(0.5 * pars.m * vrel, 2.0));
	return prefact * sigmaDiss(pars.spin, pars.color, pars.m, pars.T, u) * symmetry_factor * stimulated_emission;
}

double GammaBS(int spin, int color, double m, int spin_eta)
{
	int spin_fact = (spin - 1) / 2 + 1;
	double alphaS = parton_alpha(GGscale);
	double zet = zeta(color, m);
	double col_factor = 0.0;
	if (color == 3) col_factor = 1.0 / 6;
	if (color == 6) col_factor = 25.0 / 12;
	if (color == 8) col_factor = 9.0 / 4;
	double symmetry_factor = (spin - 1) % 2 == 0 ? 0.5 : 1;
	double spin_eta_factor = 1.0;
	// Extra factor for a spin-2 bound state decaying to vectors.
	if (spin_eta > 0)
	{
		spin_eta_factor = 0.0;
		if (spin == 5 || spin == 6)
			spin_eta_factor = 16.0 / 3.0;
	}
	return col_factor * spin_fact * symmetry_factor * spin_eta_factor * m * pow(alphaS, 2.0) * pow(zet, 3.0);
}

double GammaDissIntegrand(double u, Parameters pars)
{
	double E = BE(pars.color, pars.m);
	int gg = 16;
	double prefact = gg * 4 * M_PI / pow(2 * M_PI, 3);
	double z = E / pars.T;
	double distr = pow(E, 3) * pow(1 + u / z, 2) / (z * (exp(z + u) - 1));
	double integrand = prefact * distr * sigmaDiss(pars.spin, pars.color, pars.m, pars.T, u);
	return integrand;
}

double GammaDiss(int spin, int color, double m, double T)
{
	double E = BE(color, m);
	double z = E / T;
	double zet = zeta(color, m);
	double upper_u = z / 4.0 / pow(zet, 2.0);
	Parameters pars = {spin, color, m, T};
	double gamma = simpsonArg(GammaDissIntegrand, pars, 0, upper_u, 1e-6);
	return gamma;
}

double sigmaBSFaveraged(int spin, int color, double m, double T)
{
	double E = BE(color, m);
	double z = E / T;
	double zet = zeta(color, m);
	Parameters pars = {spin, color, m, T};
	double sigma = simpsonArg(s_integrand_BSF, pars, 0, 1, 1e-6);
	return sigma;
}

double GammaBSaveraged(int spin, int color, double m, int spin_eta, double T)
{
	double E = BE(color, m);
	double mBS = 2 * m - E;
	double ratio = bessK1(mBS / T) / bessK2(mBS / T);
	if (isnan(ratio))
		ratio = 1;
	return GammaBS(spin, color, m, spin_eta) * ratio;
}

double bound_state_rate(int spin, int color, double m, double T)
{
	// This code includes the rate for the formation of spin-2 bound states for vectors.
	double gamma_bs_spin0 = GammaBSaveraged(spin, color, m, 0, T);
	double bound_state_rate_spin0 = sigmaBSFaveraged(spin, color, m, T) * gamma_bs_spin0 / (gamma_bs_spin0 + GammaDiss(spin, color, m, T));
	double bound_state_rate_spin2 = 0.0;
	if (spin == 5 || spin == 6)
	{
		double gamma_bs_spin2 = GammaBSaveraged(spin, color, m, 2, T);
		bound_state_rate_spin2 = 25 * sigmaBSFaveraged(spin, color, m, T) * gamma_bs_spin2 / (gamma_bs_spin2 + GammaDiss(spin, color, m, T));
	}
	return bound_state_rate_spin0 + bound_state_rate_spin2;
}

