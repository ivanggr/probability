#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* возьмем немного констант из других источников */
#define M_E          2.7183
#define M_LOG2E      1.4427
#define M_LOG10E     0.4343
#define M_LN2        0.69315
#define M_LN10       2.30259
#define M_PI         3.14159
#define M_PI_2       1.57080
#define M_PI_4       0.78540
#define M_1_PI       0.31831
#define M_2_PI       0.63662
#define M_2_SQRTPI   1.12838
#define M_SQRT2      1.41421
#define M_SQRT1_2    0.70711

/* Вычисление факториала и числа сочетаний */
double fStirling(long n)
/* приближенное вычисление факториала по формуле Стирлинга n! ~= sqrt(2.0*pi*n)*pow(n/e,n) */
{
	double d;
	d = sqrt(2.0*M_PI*n) * pow(n/M_E, n);
	if (d == 0.0) 
		d = 1.0;
	return d;
}
double fact(int n)
/* вычислить факториал n */
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
/* вычислить число сочетаний n элементов по m */
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

/* Формулы Бернулли и Пуассона */
double fBernulli(int n, int k, double p)
/* Формула Бернулли. Вычисляет вероятность наступления случайного события k раз в 
серии из n испытаний, вероятность наступления каждого события постоянна и равна p */
{
	/* проверка аргументов */
	if (n <= 0 || k < 0 || k > n || p > 1 || p < 0)
		return -1.0;
	/* объявление переменных */
	double q, c, r, rp, rq;
	/* вычисление вероятности и возврат */
	q = 1.0 - p;
	rp = pow(p, k);			// rp = p^k
	rq = pow(q, n-k);		// rq = q^(n-k)
	c = combination(n, k);	// c = C(n,k)
	r = c*rp*rq;
	return r;
}
double fPuasson(long n, long k, double p, double lambda = 0)
/* Формула Пуассона. Вычисляет вероятность наступления события k раз в серии из n испытаний.
Вероятность p должна быть мала, а число испытаний n велико */
{
	if (n < 0 || p < 0.0 || p > 1.0)
		return 0.0;
	double r;
	if (lambda == 0)
		lambda = n*p;
	r = pow(lambda, k)/fact(k)*pow(M_E, -lambda);
	return r;
}

/* Вычисление приближенной вероятности количества появлений события в серии Бернулли */
double calc_x(long n, long k, const double p)
/*	Вычисляет аргумент x для локальной теоремы Муавра-Лапласа и интегральной теоремы Лапласа. 
	x = (k - n*p) / sqrt(n*p*q).
	Выводит результаты на экран.
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
/* вычисление функции фи от x */
{
	x = exp(- (x * x) / 2) * 1 / sqrt(2 * M_PI);
	return x;
}
double LocalLaplas(long n, long k, const double p)
/* Вычисляет приблизительную вероятность согласно локальной теореме Муавра-Лапласа */
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
/* Вычисляет вероятность наступления события от k1 до k2 раз в серии из n испытаний.
Внимание! Не вычисляет значение функции Лапласа. Требует задания вручную */
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
/* Вычисляет значение аргумента для функции Лапласа при нахождении вероятности отклонения 
относительной частоты события m/n от вероятности наступления события p за пределы не более eps */
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

/* Обобщенная формула для вычисления вероятности в серии Бернулли */
double sBernulli(int n, int k, double p)
{
	double r;
	r = n > 20 ? LocalLaplas(n, k, p) : fBernulli(n, k, p);
	return r;
}

/* Функции вывода таблиц */
int print_comb_table(int n)
/* печатает таблицу числа сочетаний n элементов по m начиная с 0 до n */  
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
/* напечатать биноминальный закон распределения дискретной случайной величины */
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
/* Печатает таблицу значений вероятности количества появлений события по формуле Пуассона
до количества событий не более maxn */
{
	int k;
	double fp;
	double sp = 0; // счетчик накопленной вероятности
	printf("k\tp\n");
	for (k = 0; k <= maxn; k++)
	{
		fp = fPuasson(n, k, p, lambda);
		printf("%d\t%.4f\n", k, fp);
		sp += fp;
	}
	/* вывод накопленной вероятности */
	printf("0..%d\t%.4f\n", k-1, sp); 
	/* вывод остатка вероятности */
	sp = 1.0 - sp;
	printf(">%d\t%.4f\n", k-1, sp);
	return k;
}

/* Случайные величины и их численные характеристики */
/* Класс для работы с дискретной случайной величиной */
class discrete_distr
{
	int n;			// количество возможных значений дискретной случайной величины
	char c;			// символ случайной величины
	double *v;		// массив значений дискретной величины
	double *v2;		// массив квадратов значений дискретной величины
	double *p;		// массив вероятностей появления значения дискретной величины
	double sp;		// суммарное значение всех вероятностей
	double M;		// математическое ожидание сл.величины
	double M2;		// математическое ожидание квадрата сл.величины
	double D;		// дисперсия дискретной случайной величины
	double D_2;		// дисперсия дискретной случайной величины, вычисленная по определению
	double sigma;	// среднее квадратическое отклонение
	double *dev;	// массив отклонений значений дискретной величины от матожидания
	double *dev2;	// массив квадратов отклонений значений дискретной величины от матожидания
public:
	discrete_distr(int np = 0, char nc = 'X');
	~discrete_distr();
	inline int getn() { return this->n; } // вернуть количество возможных значений		
	double getv(int np); // вернуть определенное значение 
	double setv(int nv, double newv); // установить определенное значение
	double getp(int np); // вернуть значение вероятности для определенного значения
	double setp(int np, double newp); // установить значение вероятности для определенного значения
	void set_bin_distr(double prob); // установить биноминальное распределение
	double calc_sump();	// подсчет и возврат суммы вероятностей всех возможных значений 
	double calc_M();	// подсчет и возврат математического ожидания и её квадрата
	double calc_D();	// подсчет и возврат дисперсии по формуле мат.ожиданий
	double calc_dev();	// подсчет отклонений и квадратов отклонений, дисперсии по определению
	void calc_all();	// подсчет всех численных характеристик случаной величины
	void print_vp();	// печать закона распределения 
	void print_nc();	// печать численных характеристик
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
double discrete_distr::calc_M()	// подсчет и возврат математических ожиданий величины и её квадраты
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
double discrete_distr::calc_D()	// подсчет и возврат дисперсии по формуле мат.ожиданий
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

/* Функция общая */
int main()
{
	/* создать экземпляр класса */
	const int n = 9; // число возможных значений
	discrete_distr *X;
	X = new discrete_distr(n);

	/* задать возможные значения случаной величины */
	X->set_bin_distr(0.2);

	/* задать вероятности */
	
	/* вычислить численные характеристики */
	X->calc_all();

	/* вывести закон распределения дискретной случайной величины и её численные характеристики */
	X->print_vp();
	X->print_nc();

	/* удалить экземпляр класса */
	delete X;

	/* возврат */
	return 0;
}
