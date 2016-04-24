#include <stdio.h>
#include <gmp.h>

struct qures {
	mpz_t a; /* number */
	mpz_t ainv; /* a inverse */
	mpz_t p; /* prime */
	mpz_t pm1; /* p-1 */
	mpz_t r; /* a^((s+1)/2) */
	mpz_t s; /* odd */
	mpz_t b; /* n^s */
	mpz_t n; /* odd */
	mpz_t alpha;
	mpz_t topal; /* 2^alpha */
	mpz_t j;
	mpz_t qur;
};

typedef struct qures QURES[1];

void
init_qures(QURES q)
{
	mpz_init(q->a);
	mpz_init(q->ainv);
	mpz_init(q->p);
	mpz_init(q->pm1);
	mpz_init(q->r);
	mpz_init(q->s);
	mpz_init(q->b);
	mpz_init(q->n);
	mpz_init(q->alpha);
	mpz_init(q->topal);
	mpz_init(q->j);
	mpz_init(q->qur);
}

void
clear_qures(QURES q)
{
	mpz_clear(q->a);
	mpz_clear(q->ainv);
	mpz_clear(q->p);
	mpz_clear(q->pm1);
	mpz_clear(q->r);
	mpz_clear(q->s);
	mpz_clear(q->b);
	mpz_clear(q->n);
	mpz_clear(q->alpha);
	mpz_clear(q->topal);
	mpz_clear(q->j);
	mpz_clear(q->qur);
}

void
set_qures_str(QURES q, const char *stra, const char *strn, const char *strp)
{
	mpz_t tmp;
	mpz_t count;
	mpz_t two;

	mpz_set_str(q->p, strp, 10);
	mpz_set_str(q->a, stra, 10);
	mpz_set_str(q->n, strn, 10);
	mpz_init_set_ui(count, 0);
	mpz_init_set_ui(two, 2);

	/* ainv */
	mpz_invert(q->ainv, q->a, q->p);

	mpz_sub_ui(q->pm1, q->p, 1); /* pm1 = p - 1 */
	mpz_init_set(tmp, q->pm1);

	/* calc alpha and s */
	while (mpz_even_p(tmp)) {
		mpz_divexact_ui(tmp, tmp, 2);
		mpz_add_ui(count, count, 1);
	}
	mpz_set(q->alpha, count);
	mpz_set(q->s, tmp);

	/* calc 2^alpha */
	mpz_powm(q->topal, two, q->alpha, q->p);

	/* calc r = a^((s+1)/2)*/
	mpz_add_ui(tmp, q->s, 1);
	mpz_divexact_ui(tmp, tmp, 2);
	mpz_powm(q->r, q->a, tmp, q->p);

	/* calc b */
	mpz_powm(q->b, q->n, q->s, q->p);

	mpz_set_ui(q->j, 0);

	mpz_clear(tmp);
	mpz_clear(count);
	mpz_clear(two);
}

void
calc_j(QURES q)
{
	mpz_t ra;
	mpz_t res;
	mpz_t exp;
	mpz_t k;
	mpz_t tmp;
	mpz_t two;
	mpz_init(ra);
	mpz_init(res);
	mpz_init(exp);
	mpz_init(tmp);
	mpz_init_set_ui(two, 2);
	mpz_init_set_ui(k, 1);

	/*
	 * TODO
	 * should I make mpz_mulm()?
	 */

	/* r^2 * a^-1 */
	mpz_powm_ui(ra, q->r, 2, q->p);
	mpz_mul(ra, ra, q->ainv);
	mpz_mod(ra, ra, q->p);

	/* exp = alpha - 2 */
	mpz_sub_ui(exp, q->alpha, 2);

	/* calc first j */
	mpz_powm_ui(res, q->r, 2, q->p);
	mpz_mul(res, res, q->ainv);
	mpz_mod(res, res, q->p);
	mpz_powm(tmp, two, exp, q->pm1);
	mpz_powm(res, res, tmp, q->p);

	mpz_sub_ui(exp, exp, 1);
	/* (r^2/a)^(2^(alpha - 2)) */
	if (mpz_cmp_ui(res, 1) != 0) {
		mpz_powm(tmp, two, k, q->pm1);
		mpz_add(q->j, q->j, tmp);
	}
	mpz_add_ui(k, k, 1);

	while (mpz_cmp_ui(exp, 0) >= 0) {
		/* ((b^2j)*(r^2/a))^(2^(alpha-k-1)) */
		mpz_powm(res, q->b, q->j, q->p);
		mpz_mul(res, res, ra);
		mpz_mod(res, res, q->p);
		mpz_powm(tmp, two, exp, q->pm1);
		mpz_powm(res, res, tmp, q->p);
		if (mpz_cmp_ui(res, 1) != 0) {
			mpz_powm(tmp, two, k, q->pm1);
			mpz_add(q->j, q->j, tmp);
		}
		mpz_sub_ui(exp, exp, 1);
		mpz_add_ui(k, k, 1);
	}

	mpz_divexact_ui(q->j, q->j, 2);

	mpz_clear(ra);
	mpz_clear(res);
	mpz_clear(exp);
	mpz_clear(k);
	mpz_clear(tmp);
	mpz_clear(two);
}

void
calc_qures(QURES q)
{
	mpz_powm(q->qur, q->b, q->j, q->p);
	mpz_mul(q->qur, q->qur, q->r);
	mpz_mod(q->qur, q->qur, q->p);
}

void
print_qures(QURES q)
{
	gmp_printf("%5s: %Zd\n", "p", q->p);
	gmp_printf("%5s: %Zd\n", "pm1", q->pm1);
	gmp_printf("%5s: %Zd\n", "a", q->a);
	gmp_printf("%5s: %Zd\n", "a^-1", q->ainv);
	gmp_printf("%5s: %Zd\n", "alpha", q->alpha);
	gmp_printf("%5s: %Zd\n", "n", q->n);
	gmp_printf("%5s: %Zd\n", "s", q->s);
	gmp_printf("%5s: %Zd\n", "b", q->b);
	gmp_printf("%5s: %Zd\n", "r", q->r);
	gmp_printf("%5s: %Zd\n", "j", q->j);
	gmp_printf("%5s: %Zd\n", "qur", q->qur);
}

int
main()
{
	QURES q;
	init_qures(q);

	//set_qures_str(q, "186", "3", "401");
	set_qures_str(q, "90", "3", "401");
	calc_j(q);
	calc_qures(q);
	print_qures(q);

	mpz_t tmp;
	mpz_init(tmp);
	mpz_powm_ui(tmp, q->qur, 2, q->p);
	gmp_printf("%Zd^2 = %Zd\n", q->qur, tmp);
	mpz_clear(tmp);

	clear_qures(q);
}
