#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Definitions of useful mathematical constants
 * M_E        -	e
 * M_LOG2E    -	log2(e)
 * M_LOG10E   -	log10(e)
 * M_LN2      -	ln(2)
 * M_LN10     -	ln(10)
 * M_PI       -	pi
 * M_PI_2     -	pi/2
 * M_PI_4     -	pi/4
 * M_1_PI     -	1/pi
 * M_2_PI     -	2/pi
 * M_2_SQRTPI -	2/sqrt(pi)
 * M_SQRT2    -	sqrt(2)
 * M_SQRT1_2  -	1/sqrt(2)
 */

#define M_E        2.71828182845904523536
#define M_LOG2E    1.44269504088896340736
#define M_LOG10E   0.434294481903251827651
#define M_LN2      0.693147180559945309417
#define M_LN10     2.30258509299404568402
#define M_PI       3.14159265358979323846
#define M_PI_2     1.57079632679489661923
#define M_PI_4     0.785398163397448309616
#define M_1_PI     0.318309886183790671538
#define M_2_PI     0.636619772367581343076
#define M_2_SQRTPI 1.12837916709551257390
#define M_SQRT2    1.41421356237309504880
#define M_SQRT1_2  0.707106781186547524401

/* ���������� ���������� � ����� ��������� */
double fStirling(long n)
/* ������������ ���������� ���������� �� ������� ��������� n! ~= sqrt(2.0*pi*n)*pow(n/e,n) */
{
	double d;
	d = sqrt(2.0*M_PI*n) * pow(n/M_E, n);
	if (d == 0.0) 
		d = 1.0;
	return d;
}
double fact(int n)
/* ��������� ��������� n */
{
	const long maxn = 13;
	if (n < 0)
		return 1.0;
	else if (n > maxn)
		return fStirling(n);
	double factn [maxn+1] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 
		3628800, 39916800, 479001600, 6227020800 };
	return factn[n];
}
double combination(int n, int m)
/* ��������� ����� ��������� n ��������� �� m */
{
	if (n < 0 || m < 0 || n < m)
		return 0.0;
	double fn, fm, fmn;
	fn = fact(n);
	fm = fact(m);
	fmn = fact(n-m);
	fn = fn/(fm*fmn);
	return fn;
}

/* ������� �������� � �������� */
double fBernulli(int n, int k, double p)
/* ������� ��������. ��������� ����������� ����������� ���������� ������� k ��� � 
����� �� n ���������, ����������� ����������� ������� ������� ��������� � ����� p */
{
	/* �������� ���������� */
	if (n <= 0 || k < 0 || k > n || p > 1 || p < 0)
		return -1.0;
	/* ���������� ���������� */
	double q, c, r, rp, rq;
	/* ���������� ����������� � ������� */
	q = 1.0 - p;
	rp = pow(p, k);			// rp = p^k
	rq = pow(q, n-k);		// rq = q^(n-k)
	c = combination(n, k);	// c = C(n,k)
	r = c*rp*rq;
	return r;
}
double fPuasson(long n, long k, double p, double lambda = 0)
/* ������� ��������. ��������� ����������� ����������� ������� k ��� � ����� �� n ���������.
����������� p ������ ���� ����, � ����� ��������� n ������ */
{
	if (n < 0 || p < 0.0 || p > 1.0)
		return 0.0;
	double r;
	if (lambda == 0)
		lambda = n*p;
	r = pow(lambda, k)/fact(k)*pow(M_E, -lambda);
	return r;
}

/* ���������� ������������ ����������� ���������� ��������� ������� � ����� �������� */
double calc_x(long n, long k, const double p)
/*	��������� �������� x ��� ��������� ������� ������-������� � ������������ ������� �������. 
	x = (k - n*p) / sqrt(n*p*q).
	������� ���������� �� �����.
*/
{
	double np, npq, sqrtnpq, knp, x, q;
	q = 1.0 - p;
	np = n*p;
	npq = np*q;
	sqrtnpq = sqrt(npq);
	knp = k - n*p;
	x = knp / sqrtnpq;
	printf("x = (k-n*p)/sqrt(n*p*q) = ");
	printf("(%d-%d*%.2f)/sqrt(%d*%.2f*%.2f) = \n = ", k, n, p, n, p, q);
	printf("(%d-%.2f)/sqrt(%.2f) = ", k, np, npq);
	printf("%.2f/%.2f = ", knp, sqrtnpq);
	printf("%.2f\n", x);	
	putchar('\n');
	return x;
}	
double phi_x(double x)
/* ���������� ������� �� �� x */
{
	x = exp(- (x * x) / 2) * 1 / sqrt(2 * M_PI);
	return x;
}
double LocalLaplas(long n, long k, const double p)
/* ��������� ��������������� ����������� �������� ��������� ������� ������-������� */
{
	if (n <= 0 || p <= 0 || p >= 1)
		return -1;
	double q;
	q = 1.0 - p;
	double npq;
	double x, phi, y;
	npq = sqrt(n * p * q);
	printf("sqrt(n*p*q) = sqrt(%d*%.2f*%.2f) = %.2f\n", n, p, q, npq);
	x = (k - n * p) / npq;
	printf("x = (k - n*p)/sqrt(n*p*q) = (%d - %d*%.2f)/%.2f = %.2f\n", k, n, p, npq, x);
	phi = phi_x(x);
	printf("phi(x) = phi(%.2f) = %.4f\n", x, phi);
	y = phi / npq;
	printf("P(n,k) = phi(x)/sqrt(n*p*q) = %.4f/%.2f = %.5f\n", phi, npq, y);	
	return y;
}
double IntegralLaplas(long n, long k1, long k2, const double p)
/* ��������� ����������� ����������� ������� �� k1 �� k2 ��� � ����� �� n ���������.
��������! �� ��������� �������� ������� �������. ������� ������� ������� */
{
	if (n <= 0 || k2 < k1 || k1 > n || k2 > n || p < 0 || p > 1)
		return -1;
	double q, r;
	double F2, F1, x1, x2;
	q = 1.0 - p;
	x1 = calc_x(n, k1, p);
	x2 = calc_x(n, k2, p); 
	F1 = -0.5;
	F2 = 0.5;
	printf("x1 = %.2f\t F1 = %.4f\n", x1, F1);
	printf("x2 = %.2f\t F2 = %.4f\n", x2, F2);
	putchar('\n');
	r = F2 - F1;
	printf("P%d(%d,%d) = %.3f - %.3f = %.4f\n", n, k1, k2, F2, F1, r);
	putchar('\n');
	return r;
}
double calc_epsx(long n, double p, double eps)
/* ��������� �������� ��������� ��� ������� ������� ��� ���������� ����������� ���������� 
������������� ������� ������� m/n �� ����������� ����������� ������� p �� ������� �� ����� eps */
{
	if (n <= 0 || p < 0 || p > 1 || eps < 0 || eps > 1)
		return -1;	
	double q, x;
	double npq, sqrtnpq, epsn;
	q = 1.0 - p;
	if (n == 0 || p == 0 || q == 0)
		return -1;
	epsn = eps*n;
	npq = (n*p*q);
	sqrtnpq = sqrt(npq);
	x = (eps*n)/sqrt(n*p*q);
	printf("x = (eps*n/sqrt(n*p*q) = (%.3f*%d)/sqrt(%d*%.3f*%.3f) = \n = ", eps, n, n, p, q);
	printf("%.3f/sqrt(%.3f) = ", epsn, npq);	
	printf("%.3f/%.3f = ", epsn, sqrtnpq);		
	printf("%.3f\n", x);	
	return x;
}

/* ���������� ������� ��� ���������� ����������� � ����� �������� */
double sBernulli(int n, int k, double p)
{
	double r;
	r = n > 20 ? LocalLaplas(n, k, p) : fBernulli(n, k, p);
	return r;
}

/* ������� ������ ������ */
int print_comb_table(int n)
/* �������� ������� ����� ��������� n ��������� �� m ������� � 0 �� n */  
{
	int m;
	long Cnm;
	for (m = 0; m <= n; m ++)
	{
		Cnm = (long) combination(n, m);
		printf("C(n,m) = C(%d,%d) = %d\n", n, m, Cnm);
	}	
	return m;
}
int print_bin_distr(int n, double p)
/* ���������� ������������� ����� ������������� ���������� ��������� �������� */
{
	if (n <= 0 || p < 0 || p > 1)
		return -1;
	int k;
	double r;
	printf("X\t");
	for (k = 0; k <= n; k++)
		printf("%d\t", k);
	printf("\np\t");
	for (k = 0; k <= n; k++)
	{
		r = sBernulli(n, k, p);
		printf("%.4f\t", r);
	}
	printf("\n");
	return k;
}
int print_Puasson_table(long n, long maxn, double p, double lambda = 0)
/* �������� ������� �������� ����������� ���������� ��������� ������� �� ������� ��������
�� ���������� ������� �� ����� maxn */
{
	int k;
	double fp;
	double sp = 0; // ������� ����������� �����������
	printf("k\tp\n");
	for (k = 0; k <= maxn; k++)
	{
		fp = fPuasson(n, k, p, lambda);
		printf("%d\t%.4f\n", k, fp);
		sp += fp;
	}
	/* ����� ����������� ����������� */
	printf("0..%d\t%.4f\n", k-1, sp); 
	/* ����� ������� ����������� */
	sp = 1.0 - sp;
	printf(">%d\t%.4f\n", k-1, sp);
	return k;
}

/* ��������� �������� � �� ��������� �������������� */
/* ����� ��� ������ � ���������� ��������� ��������� */
class discrete_distr
{
	int n;			// ���������� ��������� �������� ���������� ��������� ��������
	char c;			// ������ ��������� ��������
	double *v;		// ������ �������� ���������� ��������
	double *v2;		// ������ ��������� �������� ���������� ��������
	double *p;		// ������ ������������ ��������� �������� ���������� ��������
	double sp;		// ��������� �������� ���� ������������
	double M;		// �������������� �������� ��.��������
	double M2;		// �������������� �������� �������� ��.��������
	double D;		// ��������� ���������� ��������� ��������
	double D_2;		// ��������� ���������� ��������� ��������, ����������� �� �����������
	double sigma;	// ������� �������������� ����������
	double *dev;	// ������ ���������� �������� ���������� �������� �� �����������
	double *dev2;	// ������ ��������� ���������� �������� ���������� �������� �� �����������
public:
	discrete_distr(int np = 0, char nc = 'X');
	~discrete_distr();
	inline int getn() { return this->n; } // ������� ���������� ��������� ��������		
	double getv(int np); // ������� ������������ �������� 
	double setv(int nv, double newv); // ���������� ������������ ��������
	double getp(int np); // ������� �������� ����������� ��� ������������� ��������
	double setp(int np, double newp); // ���������� �������� ����������� ��� ������������� ��������
	void set_bin_distr(double prob); // ���������� ������������� �������������
	double calc_sump();	// ������� � ������� ����� ������������ ���� ��������� �������� 
	double calc_M();	// ������� � ������� ��������������� �������� � � ��������
	double calc_D();	// ������� � ������� ��������� �� ������� ���.��������
	double calc_dev();	// ������� ���������� � ��������� ����������, ��������� �� �����������
	void calc_all();	// ������� ���� ��������� ������������� �������� ��������
	void print_vp();	// ������ ������ ������������� 
	void print_nc();	// ������ ��������� �������������
};

discrete_distr::discrete_distr(int np, char nc)
{
	this->n = np;
	this->c = nc;
	this->v = (double *) calloc((size_t) this->n, sizeof (double));
	this->v2 = (double *) calloc((size_t) this->n, sizeof (double));
	this->dev = (double *) calloc((size_t) this->n, sizeof (double));
	this->dev2 = (double *) calloc((size_t) this->n, sizeof (double));
	this->p = (double *) calloc((size_t) this->n, sizeof (double));
	return;
}
discrete_distr::~discrete_distr()
{
	free(this->v);
	free(this->p);
	return;
}
double discrete_distr::getv(int np) 
{
	if (this->v != NULL)
		return this->v[np];
	else 
		return 0.0;
}
double discrete_distr::setv(int nv, double newv)
{
	if (nv < 0 || nv >= this->n)
		return 0.0;
	this->v[nv] = newv;
	this->v2[nv] = newv * newv;
	return this->v[nv];
}
double discrete_distr::getp(int np) 
{
	if (this->p != NULL)
		return this->p[np];
	else 
		return 0.0;
}
double discrete_distr::setp(int np, double newp)
{
	if (np < 0 || np >= this->n)
		return 0.0;
	this->sp -= this->p[np];
	this->p[np] = newp;
	this->sp += this->p[np];
	return this->p[np];
}
void discrete_distr::set_bin_distr(double prob)
{
	int i;
	double set_p;
	for (i = 0; i < n; i++)
	{
		setv(i, (double) i);	
		set_p = sBernulli(n-1, i, prob);
		setp(i, set_p);
	}
	return;
}
double discrete_distr::calc_sump()
{
	int k;
	sp = 0;
	for (k = 0; k < n; k++)
		sp += p[k];
	return sp;
}
double discrete_distr::calc_M()	// ������� � ������� �������������� �������� �������� � � ��������
{
	int k;
	M = 0;
	M2 = 0;
	for (k = 0; k < n; k++)
	{
		M += p[k]*v[k];
		M2 += p[k]*v2[k];
	}
	return M;
}
double discrete_distr::calc_D()	// ������� � ������� ��������� �� ������� ���.��������
{
	D = M2 - M*M;
	return D;
}
double discrete_distr::calc_dev()
{
	int j;
	D_2 = 0;
	for (j = 0; j < n; j++)
	{
		dev[j] = v[j] - M;
		dev2[j] = dev[j] * dev[j];
		D_2 += dev2[j] * p[j];
	}
	sigma = D_2 >= 0 ? sqrt(D_2) : 0;
	return D_2;
}
void discrete_distr::calc_all()
{
	calc_sump();
	calc_M();
	calc_D();
	calc_dev();
	return;
}
void discrete_distr::print_vp()
{
	int j;
	printf(" %c | ", c);
	for (j = 0; j < n; j++)
		printf("%.3f | ", v[j]);
	putchar('\n');
	printf(" p | ");
	for (j = 0; j < n; j++)
		printf("%.3f | ", p[j]);
	putchar('\n');
	return;
}
void discrete_distr::print_nc()
{
	printf("sum probability | %.3f\n", sp);
	printf("expectation | %6.3f\n", M);
	printf("expectation of square | %6.3f\n", M2);
	printf("dispersion by formula | %6.3f\n", D);
	printf("dispersion by definition | %6.3f\n", D_2);
	printf("sigma | %6.3f\n", sigma);
	return;
}

/* ������� ����� */
int main()
{
	/* ������� ��������� ������ */
	const int n = 9; // ����� ��������� ��������
	discrete_distr *X;
	X = new discrete_distr(n);

	/* ������ ��������� �������� �������� �������� */
	X->set_bin_distr(0.2);

	/* ������ ����������� */
	
	/* ��������� ��������� �������������� */
	X->calc_all();

	/* ������� ����� ������������� ���������� ��������� �������� � � ��������� �������������� */
	X->print_vp();
	X->print_nc();

	/* ������� ��������� ������ */
	delete X;

	/* ������� */
	return 0;
}