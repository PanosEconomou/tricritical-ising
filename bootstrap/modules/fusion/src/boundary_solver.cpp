#include <glpk.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/arb.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/nf_elem.h>
#include <flint/qqbar.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "boundary_solver.h"

#define BOUND_PREC 256
#define NUM_PREC   128
#define LP_RELAX   1e-9   /* positivity rows of the float LPs are relaxed by this.
                             Affects speed only: every range and every pruning
                             decision is certified, see certified_lower */
#define INDEX_BASE 1      /* Ishibashi index files come from Mathematica: 1-based */

/* Reflection and transmission of the folded interface (Quella-Runkel-Watts), for two
   copies with equal central charge. Labels as in the input files (with INDEX_BASE):
   TRANS_DEN is the vacuum (1,1)+, TRANS_NUM the exchange-odd part of the vacuum
   module (1,1)-, whose lowest state is (L^(1)_{-2} - L^(2)_{-2})|0>. With unit
   Ishibashi states and r = b^TRANS_NUM / b^TRANS_DEN (signed),
       R = (1 + r)/2,   T = (1 - r)/2.
   Only computed after disentangling, where the signs are known (see the note
   printed with the boundaries). */
#define TRANS_NUM  35
// #define TRANS_NUM  17
#define TRANS_DEN  1

/* 1: every Ishibashi state in I must appear, |b^i| > 0 ("= 0 iff i not in I").
   0: |b^i| >= 0 for i in I, so boundaries using a subset of I are kept too
      (e.g. topological interfaces in the folded picture). */
#define STRICT_POSITIVITY 0

static inline const nf_elem_struct* entry(const gr_mat_struct* S, slong i, slong j,
                                          gr_ctx_struct* ctx)
{
    return (const nf_elem_struct*) gr_mat_entry_srcptr(S, i, j, ctx);
}

/* Real value of a at the real embedding given by gen. */
static void nf_elem_get_arb(arb_t v, const nf_elem_t a, const arb_t gen,
                            const nf_t field, slong prec)
{
    fmpq_poly_t p; fmpq_poly_init(p);
    nf_elem_get_fmpq_poly(p, a, field);
    _arb_fmpz_poly_evaluate_arb(v, fmpq_poly_numref(p), fmpq_poly_length(p), gen, prec);
    arb_div_fmpz(v, v, fmpq_poly_denref(p), prec);
    fmpq_poly_clear(p);
}

/* Exact sign of a (0, +1, -1); an enclosure of its value is left in v. */
static int nf_elem_sign(arb_t v, const nf_elem_t a, const qqbar_t generator,
                        const nf_t field)
{
    if (nf_elem_is_zero(a, field)) { arb_zero(v); return 0; }

    arb_t g; arb_init(g);
    for (slong prec = 64; ; prec *= 2) {      /* terminates since a != 0 */
        qqbar_get_arb(g, generator, prec);
        nf_elem_get_arb(v, a, g, field, prec);
        if (!arb_contains_zero(v)) break;
    }
    arb_clear(g);
    return arb_is_positive(v) ? 1 : -1;
}

/* Sign of sum_j S_ij x^j: numeric first, exact only when the ball contains zero. */
static int row_sign(arb_t out, const gr_mat_struct* S, slong i, const slong* x,
                    slong n, arb_srcptr Sa, nf_struct* field, qqbar_struct* gen,
                    gr_ctx_struct* ctx)
{
    arb_zero(out);
    for (slong j = 0; j < n; j++)
        if (x[j]) arb_addmul_si(out, Sa + i * n + j, x[j], NUM_PREC);
    if (arb_is_positive(out)) return 1;
    if (arb_is_negative(out)) return -1;

    nf_elem_t y, t; nf_elem_init(y, field); nf_elem_init(t, field);
    nf_elem_zero(y, field);
    for (slong j = 0; j < n; j++) {
        if (x[j] == 0) continue;
        nf_elem_scalar_mul_si(t, entry(S, i, j, ctx), x[j], field);
        nf_elem_add(y, y, t, field);
    }
    int sgn = nf_elem_sign(out, y, gen, field);
    nf_elem_clear(t, field); nf_elem_clear(y, field);
    return sgn;
}

/*
 * With x^0 = 1, x is a convex combination of v_j = (S_ij / S_0j)_i, j in I, so
 *     max(0, ceil(min_j S_ij/S_0j)) <= x^i <= floor(max_j S_ij/S_0j),
 * rounded outward so no true solution is cut off.
 * Returns 1 if the box is empty, -1 on error, 0 otherwise.
 */
static int build_bounds(slong* lb, slong* ub, const gr_mat_struct* S, const char* inI,
                        nf_struct* field, qqbar_struct* generator, gr_ctx_struct* ctx)
{
    slong n = gr_mat_ncols(S, ctx);
    int status = 0;

    nf_elem_t r;      nf_elem_init(r, field);
    arb_t g, v;       arb_init(g);  arb_init(v);
    arf_t t, mn, mx;  arf_init(t);  arf_init(mn); arf_init(mx);

    qqbar_get_arb(g, generator, BOUND_PREC);

    for (slong i = 0; i < n && status == 0; i++) {
        arf_pos_inf(mn);
        arf_neg_inf(mx);
        for (slong j = 0; j < n; j++) {
            if (!inI[j]) continue;
            if (nf_elem_is_zero(entry(S, 0, j, ctx), field)) {
                flint_fprintf(stderr, "ERROR: S_0%wd is zero\n", j);
                status = -1;
                break;
            }
            nf_elem_div(r, entry(S, i, j, ctx), entry(S, 0, j, ctx), field);
            nf_elem_get_arb(v, r, g, field, BOUND_PREC);
            arb_get_ubound_arf(t, v, BOUND_PREC); arf_max(mx, mx, t);
            arb_get_lbound_arf(t, v, BOUND_PREC); arf_min(mn, mn, t);
        }
        if (status) break;

        ub[i] = arf_get_si(mx, ARF_RND_FLOOR);
        lb[i] = FLINT_MAX(0, arf_get_si(mn, ARF_RND_CEIL));
        if (lb[i] > ub[i]) status = 1;
    }

    arf_clear(mx); arf_clear(mn); arf_clear(t);
    arb_clear(v);  arb_clear(g);
    nf_elem_clear(r, field);
    return status;
}

/* ------------------------------------------------------------------------ */
/* LP-guided depth-first search                                             */
/* ------------------------------------------------------------------------ */

typedef struct {
    glp_prob*          lp;
    glp_smcp           parm;
    slong              n;
    std::vector<slong> order;   /* free columns in branching order */
    std::vector<slong> lb, ub;  /* global box */
    std::vector<slong> cl, cu;  /* current bounds: the box, or lb = ub = x for fixed columns */
    std::vector<slong> x;       /* current assignment */

    /* certified bounds */
    glp_prob*          el;      /* elastic LP, used only to certify infeasibility */
    slong              rank, nI;
    std::vector<slong> Irow;    /* positivity row k <-> Ishibashi index Irow[k] */
    std::vector<double> y, z;   /* multipliers for the equality / positivity rows */

    /* dependents: x[piv[r]] = -(sum_{j != piv[r]} Z[r][j] x[j]) / Z[r][piv[r]] */
    fmpz_mat_struct*   Z;
    std::vector<slong> piv;
    fmpz_t             acc;

    /* exact check */
    const char*        inI;
    arb_srcptr         Sa;
    gr_mat_struct*     S;
    nf_struct*         field;
    qqbar_struct*      gen;
    gr_ctx_struct*     ctx;
    arb_t              v;

    std::vector<std::vector<slong>> sols;  /* candidates, in printed order */
    std::vector<int>   mult;               /* Ishibashi states per label: times listed in the input */
    slong              dimI;               /* sum of mult */

    slong              nodes, leaves, nsolutions;
    slong              fallbacks;          /* LP failed, box used for that end */
    slong              lp_prunes;          /* certified empty from an "optimal" LP */
    slong              elastic_calls, elastic_certified;
    double             max_gap;            /* max |GLPK optimum - certified bound| */
} search_t;

static void set_col_box(glp_prob* lp, int col, slong lo, slong hi)
{
    glp_set_col_bnds(lp, col, lo == hi ? GLP_FX : GLP_DB, (double) lo, (double) hi);
}

/*
 * Certified lower bound L on  min c^T x  over the TRUE feasible set
 *     { Z x = 0,  S_I x >= 0,  cl <= x <= cu }
 * (exact Z, ball enclosures of S), with c = sense * e_col, or c = 0 if col < 0.
 *
 * For any y and any z >= 0, with r = c - Z^T y - S_I^T z,
 *     c^T x = r^T x + y^T (Z x) + z^T (S_I x) >= r^T x >= sum_j min(r_j cl_j, r_j cu_j).
 * So L is valid for ANY multipliers; good ones (the LP duals) make it tight.
 * y is clipped to [-ymax, ymax] and z to [0, zmax] before use.
 */
static void certified_lower(arf_t L, search_t* e, slong col, int sense,
                            double ymax, double zmax)
{
    const slong n = e->n;
    arb_t r, t, m;  arb_init(r); arb_init(t); arb_init(m);
    arf_t a, b;     arf_init(a); arf_init(b);
    arf_zero(L);

    for (slong j = 0; j < n; j++) {
        arb_set_si(r, j == col ? sense : 0);

        for (slong k = 0; k < e->rank; k++) {            /* r -= Z^T y, Z exact */
            double yk = std::min(ymax, std::max(-ymax, e->y[k]));
            if (yk == 0.0 || fmpz_is_zero(fmpz_mat_entry(e->Z, k, j))) continue;
            arb_set_d(m, yk);
            arb_set_fmpz(t, fmpz_mat_entry(e->Z, k, j));
            arb_submul(r, t, m, NUM_PREC);
        }
        for (slong k = 0; k < e->nI; k++) {              /* r -= S_I^T z, S as balls */
            double zk = std::min(zmax, e->z[k]);
            if (!(zk > 0.0)) continue;                    /* z >= 0 (also drops NaN) */
            arb_set_d(m, zk);
            arb_submul(r, e->Sa + e->Irow[k] * n + j, m, NUM_PREC);
        }

        /* inf of r_j x_j over x_j in [cl_j, cu_j], taken at an endpoint */
        arb_mul_si(t, r, e->cl[j], NUM_PREC);
        arb_get_lbound_arf(a, t, NUM_PREC);
        if (e->cu[j] != e->cl[j]) {
            arb_mul_si(t, r, e->cu[j], NUM_PREC);
            arb_get_lbound_arf(b, t, NUM_PREC);
            arf_min(a, a, b);
        }
        arf_add(L, L, a, NUM_PREC, ARF_RND_FLOOR);       /* rounding down keeps it a lower bound */
    }

    arf_clear(b); arf_clear(a);
    arb_clear(m); arb_clear(t); arb_clear(r);
}

static int run_simplex(glp_prob* lp, glp_smcp* parm)
{
    int ret = glp_simplex(lp, parm);
    if (ret != 0) { glp_std_basis(lp); ret = glp_simplex(lp, parm); }
    return ret;
}

/*
 * Certify that the true feasible set (with the current bounds) is empty.
 * Elastic LP:  min sum(s+ + s- + s)  s.t.  Z x + s+ - s- = 0,  S_I x + s >= 0,
 * box on x, slacks >= 0. It is always feasible, and its optimum is 0 whenever
 * the true problem is feasible. With |y| <= 1 and 0 <= z <= 1 every slack has a
 * nonnegative coefficient in the Lagrangian, so certified_lower with c = 0
 * bounds its optimum from below; L > 0 proves infeasibility.
 * Returns 1 only if infeasibility is proven.
 */
static int certify_infeasible(search_t* e)
{
    e->elastic_calls++;
    for (slong j = 0; j < e->n; j++) set_col_box(e->el, (int) j + 1, e->cl[j], e->cu[j]);

    if (run_simplex(e->el, &e->parm) != 0 || glp_get_status(e->el) != GLP_OPT) return 0;

    for (slong k = 0; k < e->rank; k++) e->y[k] = glp_get_row_dual(e->el, (int) k + 1);
    for (slong k = 0; k < e->nI; k++)   e->z[k] = glp_get_row_dual(e->el, (int) (e->rank + k) + 1);

    arf_t L; arf_init(L);
    certified_lower(L, e, -1, 0, 1.0, 1.0);
    int proven = arf_sgn(L) > 0;
    arf_clear(L);

    e->elastic_certified += proven;
    return proven;
}

/*
 * Certified integer range [*lo, *hi] of column j over the true feasible set
 * with the current fixings. Returns 0 only if that set is PROVEN empty.
 * GLPK only proposes multipliers; when it fails, the box is used for that end.
 */
static int lp_range(search_t* e, slong j, slong* lo, slong* hi)
{
    const int col = (int) j + 1;
    int empty = 0;
    *lo = e->cl[j];
    *hi = e->cu[j];

    arf_t L; arf_init(L);
    glp_set_obj_coef(e->lp, col, 1.0);

    for (int dir = 0; dir < 2 && !empty; dir++) {
        const int sense = (dir == 0) ? 1 : -1;      /* minimize sense * x_j */
        glp_set_obj_dir(e->lp, dir == 0 ? GLP_MIN : GLP_MAX);

        if (run_simplex(e->lp, &e->parm) != 0) { e->fallbacks++; continue; }
        int st = glp_get_status(e->lp);
        if (st == GLP_NOFEAS) {                     /* not trusted: certify or keep the box */
            empty = certify_infeasible(e);
            break;
        }
        if (st != GLP_OPT) { e->fallbacks++; continue; }

        /* duals of "max x_j" are minus the duals of "min -x_j" */
        for (slong k = 0; k < e->rank; k++) e->y[k] = sense * glp_get_row_dual(e->lp, (int) k + 1);
        for (slong k = 0; k < e->nI; k++)   e->z[k] = sense * glp_get_row_dual(e->lp, (int) (e->rank + k) + 1);

        certified_lower(L, e, j, sense, INFINITY, INFINITY);
        if (!arf_is_finite(L)) { e->fallbacks++; continue; }

        double gap = std::fabs(sense * glp_get_obj_val(e->lp) - arf_get_d(L, ARF_RND_NEAR));
        if (gap > e->max_gap) e->max_gap = gap;

        if (dir == 0) {                             /* x_j >= L */
            if (arf_cmp_si(L, e->cu[j]) > 0)      { empty = 1; e->lp_prunes++; }
            else if (arf_cmp_si(L, e->cl[j]) > 0) *lo = arf_get_si(L, ARF_RND_CEIL);
        } else {                                    /* x_j <= -L */
            arf_neg(L, L);
            if (arf_cmp_si(L, *lo) < 0)           { empty = 1; e->lp_prunes++; }
            else if (arf_cmp_si(L, e->cu[j]) < 0) *hi = arf_get_si(L, ARF_RND_FLOOR);
        }
    }

    glp_set_obj_coef(e->lp, col, 0.0);
    arf_clear(L);
    return !empty;
}

static void leaf(search_t* e)
{
    slong n = e->n;
    e->leaves++;

    /* dependent coordinates, exactly */
    for (size_t r = 0; r < e->piv.size(); r++) {
        slong p = e->piv[r];
        fmpz_zero(e->acc);
        for (slong j = 0; j < n; j++)
            if (j != p && e->x[j] != 0)   /* other pivot columns are 0 in row r */
                fmpz_addmul_si(e->acc, fmpz_mat_entry(e->Z, r, j), e->x[j]);

        const fmpz* d = fmpz_mat_entry(e->Z, r, p);
        if (!fmpz_divisible(e->acc, d)) return;
        fmpz_divexact(e->acc, e->acc, d);
        fmpz_neg(e->acc, e->acc);
        if (!fmpz_fits_si(e->acc)) return;

        slong val = fmpz_get_si(e->acc);
        if (val < e->lb[p] || val > e->ub[p]) return;
        e->x[p] = val;
    }

    /* positivity, exactly (strict or not, see STRICT_POSITIVITY) */
    for (slong i = 0; i < n; i++) {
        if (!e->inI[i]) continue;
        int sg = row_sign(e->v, e->S, i, e->x.data(), n, e->Sa,
                          e->field, e->gen, e->ctx);
        if (sg < 0 || (STRICT_POSITIVITY && sg == 0)) return;
    }

    e->nsolutions++;
    e->sols.push_back(e->x);
    flint_printf("x = ");
    for (slong j = 0; j < n; j++) flint_printf("%wd ", e->x[j]);
    flint_printf("\n");
    for (slong i = 0; i < n; i++) {
        if (!e->inI[i]) continue;
        row_sign(e->v, e->S, i, e->x.data(), n, e->Sa, e->field, e->gen, e->ctx);
        arb_sqrt(e->v, e->v, NUM_PREC);
        flint_printf("  |b^%wd| = ", i + INDEX_BASE);
        arb_printn(e->v, 15, 0);
        flint_printf("\n");
    }

    /* |r| only: its sign, which decides R and T, is fixed after disentangling */
    slong inum = TRANS_NUM - INDEX_BASE, iden = TRANS_DEN - INDEX_BASE;
    if (inum >= 0 && inum < n && iden >= 0 && iden < n && e->inI[inum] && e->inI[iden]) {
        arb_t bnum, bden;
        arb_init(bnum); arb_init(bden);
        row_sign(bnum, e->S, inum, e->x.data(), n, e->Sa, e->field, e->gen, e->ctx);
        row_sign(bden, e->S, iden, e->x.data(), n, e->Sa, e->field, e->gen, e->ctx);
        if (!arb_is_zero(bden)) {
            arb_sqrt(bnum, bnum, NUM_PREC);
            arb_sqrt(bden, bden, NUM_PREC);
            arb_div(bnum, bnum, bden, NUM_PREC);
            flint_printf("  |b^%d|/|b^%d| = ", TRANS_NUM, TRANS_DEN);
            arb_printn(bnum, 15, 0);
            flint_printf("\n");
        }
        arb_clear(bden); arb_clear(bnum);
    }
    fflush(stdout);
}

static void dfs(search_t* e, size_t d)
{
    e->nodes++;
    if ((e->nodes & 0xFFFF) == 0) {
        flint_printf("  nodes %wd  leaves %wd  solutions %wd\n",
                     e->nodes, e->leaves, e->nsolutions);
        fflush(stdout);
    }

    if (d == e->order.size()) { leaf(e); return; }

    slong j   = e->order[d];
    int   col = (int) j + 1;
    slong a, b;
    if (!lp_range(e, j, &a, &b)) return;

    for (slong val = a; val <= b; val++) {
        glp_set_col_bnds(e->lp, col, GLP_FX, (double) val, (double) val);
        e->cl[j] = e->cu[j] = e->x[j] = val;
        dfs(e, d + 1);
    }
    set_col_box(e->lp, col, e->lb[j], e->ub[j]);
    e->cl[j] = e->lb[j];
    e->cu[j] = e->ub[j];
}

/* ------------------------------------------------------------------------ */
/* Disentangling: mutual Cardy conditions among the candidates              */
/* ------------------------------------------------------------------------ */
/*
 * The search only imposes the self-overlap of one boundary. Here the
 * candidates are combined into families of boundary states that satisfy
 * Cardy's condition pairwise:
 *     n_ab = S p_ab  is a nonnegative integer vector,  p_ab^i = <b_a^i, b_b^i>,
 * where b_a^i lives in R^{m_i} and m_i is the number of Ishibashi states with
 * label i: the number of times i is listed in the Ishibashi index input, i.e.
 * the diagonal mass-matrix entry Z_ii. (Real vectors suffice: n real and S real
 * force every p real, and a real Gram matrix of complex rank m is realised in R^m.)
 *
 * Facts used, all from S real and S^2 = 1:
 *   - |n_ab|^2 = |p_ab|^2 <= x_a . x_b;
 *   - distinct elementary boundaries have n_ab^0 = 0 (Cauchy-Schwarz with
 *     S_0i > 0 and x^0 = 1), so the vectors (sqrt(S_0i) b_a^i)_i of a family are
 *     orthonormal: a family has at most sum_i m_i members over the labels it
 *     uses, and a complete one is a basis.
 *
 * Gauge-free formulation. A family exists iff for every label i its Gram matrix
 * G^i_kl = (S n_kl)_i (with G^i_kk = (S x_k)_i) is positive semidefinite of rank
 * <= m_i. This is checked incrementally with an LDL^T (Gram-Schmidt without
 * square roots), exactly in the number field: per label, the members that raise
 * the rank are pivots with D-entries d_j > 0, every state carries its
 * coordinates y_j along the pivots and its residual r >= 0, and the inner
 * product of two states s, v is
 *     <s, v> = sum_j y_s,j y_v,j / d_j  +  (part in the orthogonal complement),
 * where the complement part is 0 if either residual vanishes or no complement
 * is left, +-sqrt(r_s r_v) if one dimension is left, and anything in between if
 * more are left. All acceptance tests are polynomial identities or exact sign
 * tests in the field; candidate n are only found numerically.
 *
 * Search. Families are searched per support pattern U (the labels their members
 * use): U = I for complete families, and the support of each candidate for
 * families complete on a subset. Members are added one at a time, starting from
 * none. The remaining states are kept as classes: states that are
 * indistinguishable given the members so far (same candidate, same overlaps).
 * Adding a member refines every class by its possible overlaps with it,
 * including the member's own class, which yields the states related to it by a
 * reflection (e.g. the three fixed boundaries of the Potts model). A class
 * whose states have residuals in a complement of total dimension k can supply
 * at most k + 1 members, and candidate B at most min_i m_i / (S_0i |b_B^i|^2)
 * (from the trace condition below).
 *
 * Trace condition. A complete family is an orthonormal basis, so for every
 * label sum_k S_0i b_k^i b_k^iT = 1_{m_i}; taking traces and applying S,
 *     sum_{k in F} x_k = S w,   w_i = m_i / S_0i on U.
 * This integer target is computed exactly; a partial family may only take
 * states whose x fits in what is left of it, and the remaining classes must be
 * able to cover what is left, component by component. The search branches on
 * the component still owed that the fewest classes can supply (as in exact
 * cover): every completion contains one of them, and the t-th branch drops the
 * classes before it, so each family is still found exactly once.
 *
 * Exactness. The search runs in floating point: candidate overlaps are found
 * numerically, and states are updated numerically (residuals below ZERO_TOL are
 * taken as 0). Integer conditions (Z n = 0 for labels outside I) are exact. Every
 * family the search returns is then replayed in the number field with its overlaps
 * fixed, which checks every Gram condition exactly; only families that pass are
 * reported. Soundness is therefore exact; completeness relies on the numerics,
 * as the candidate overlaps always did.
 *
 * Completion LP. At every node, with c_t the number of members still to be taken from
 * class t, three families of identities of a complete family are linear in c: the
 * trace identity, the Perron-Frobenius identity sum_l n^j_kl g_l = d_j g_k for each
 * member k (the vacuum column of Psi is an eigenvector of every n^j), and the NIM-rep
 * budget below as an equality. An infeasible LP prunes the node. The PF identity also
 * bounds every overlap: n^j_kl g_l <= d_j g_k - (what k already owes), so x^j <= d_j.
 *
 * NIM-rep budgets. In a family complete on its support, n^j = Psi diag(S_ij/S_0i)
 * Psi^T with Psi orthogonal, so (n^j)^2 = sum_r N_jj^r n^r (Verlinde) and, on the
 * diagonal, sum_{l != k} (n^j_kl)^2 = sum_r N_jj^r x_k^r - (x_k^j)^2 =: cap_k^j,
 * which depends on x_k alone. Every overlap search is boxed by
 * n^j_kl <= sqrt(what is left of cap_k^j and cap_l^j), and members' budgets are
 * tracked during the search. (Not used in the closure test, which is about
 * states outside the family.)
 */

#define DISENTANGLE        1
#define SUBSET_FAMILIES    0         /* 0: search only families complete for all of I (faster) */
#define MAX_FAMILIES_SHOWN 10
#define SEARCH_NODE_LIMIT  20000000  /* per run; a warning is printed if hit */
#define DIS_TOL            1e-7      /* numeric tolerance in the search; families are verified exactly */
#define ZERO_TOL           1e-9      /* residuals below this (relative) are taken as 0 in the search */
#define PROGRESS_SECONDS   2.0       /* minimum time between progress lines */

/* owning array of number-field elements */
struct nf_vec {
    nf_elem_struct* v = nullptr;
    slong           len = 0;
    nf_struct*      K = nullptr;
    nf_vec() {}
    nf_vec(slong l, nf_struct* k) : len(l), K(k)
    {
        v = (nf_elem_struct*) flint_malloc(FLINT_MAX(l, 1) * sizeof(nf_elem_struct));
        for (slong i = 0; i < l; i++) nf_elem_init(v + i, K);
    }
    nf_vec(nf_vec&& o) noexcept : v(o.v), len(o.len), K(o.K) { o.v = nullptr; o.len = 0; }
    nf_vec& operator=(nf_vec&& o) noexcept
    {
        if (this != &o) { release(); v = o.v; len = o.len; K = o.K; o.v = nullptr; o.len = 0; }
        return *this;
    }
    nf_vec(const nf_vec&) = delete;
    nf_vec& operator=(const nf_vec&) = delete;
    ~nf_vec() { release(); }
    void release()
    {
        if (!v) return;
        for (slong i = 0; i < len; i++) nf_elem_clear(v + i, K);
        flint_free(v);
        v = nullptr;
    }
    nf_vec clone() const
    {
        nf_vec c(len, K);
        for (slong i = 0; i < len; i++) nf_elem_set(c.v + i, v + i, K);
        return c;
    }
    nf_elem_struct* operator[](slong i) const { return v + i; }
};

/* Per-label storage: label i owns slots off[i] .. off[i] + m_i (m_i coordinates, then r). */

/* A state (or class of indistinguishable states) relative to the members so far. */
struct dstate {
    slong                           cand;
    std::vector<std::vector<slong>> nh;     /* n with each member, in order of addition */
    nf_vec                          ex;     /* exact y_j and r per label */
    std::vector<double>             nd;     /* the same, numerically */
    std::vector<char>               rz;     /* per label: r == 0 exactly */
    slong                           bound;  /* at most this many members from this class */
    std::vector<slong>              used;   /* sum over members of (n^j)^2 */
    std::vector<double>             pf;     /* sum over members l of n^j g_l */

    dstate clone() const
    {
        dstate c;
        c.cand = cand; c.nh = nh; c.ex = ex.clone(); c.nd = nd; c.rz = rz; c.bound = bound;
        c.used = used; c.pf = pf;
        return c;
    }
};

/* Pivot structure per label: count c_i and D-entries d_ij (exact, inverted, numeric).
   Inverses are kept because every projection divides by them. */
struct dframe {
    std::vector<slong>  c;
    nf_vec              d, dinv;
    std::vector<double> dd;

    dframe clone() const
    {
        dframe f;
        f.c = c; f.d = d.clone(); f.dinv = dinv.clone(); f.dd = dd;
        return f;
    }
};

/* a member of a found family */
struct member_t {
    slong                           cand;
    std::vector<std::vector<slong>> nh;
    nf_vec                          ex;     /* y, r relative to the frame when it was added */
    std::vector<slong>              c_at;   /* pivot counts at that time */
    member_t() {}
    member_t(member_t&&) = default;
    member_t& operator=(member_t&&) = default;
};

typedef std::vector<std::vector<std::vector<slong>>> nmat_t;   /* n[k][l] = n_{kl} */

struct family_t {
    std::vector<char>     U;                 /* support pattern */
    int                   kind;              /* 2: complete for I, 1: closed on a subset */
    std::vector<member_t> mem;
    nmat_t                N;
    dframe                frame;             /* final pivots, for coordinates */
};

/* Integer lattice {n : (S n)_i = 0 for i in a set of labels}, row-reduced over Q
   (each such label gives one rational equation per power-basis coefficient). */
struct dlattice {
    fmpz_mat_struct*    Z = nullptr;
    std::vector<slong>  piv, fcol;
    std::vector<std::vector<double>> W;    /* W[i][t]: (S n)_i = sum_t W[i][t] n_{fcol[t]} */
    std::vector<std::vector<double>> PR;   /* PR[r][t]: n_{piv[r]} = sum_t PR[r][t] n_{fcol[t]} */
    ~dlattice() { if (Z) { fmpz_mat_clear(Z); flint_free(Z); } }
};

struct dis_t {
    search_t*  e;
    slong      n, M, L;
    std::vector<slong>  off, mult;
    nf_elem_struct* q;                 /* q[a*n + i] = (S x_a)_i = |b_a^i|^2 exactly */
    std::vector<double> qd;            /* the same, numerically */
    int                 exact;         /* 1 while replaying a family in the number field */
    slong               rejected;      /* families that failed the exact replay */
    std::vector<double> Sd;            /* all of S, numerically */
    std::vector<double> Qs;            /* Qs[j*n + i] = sqrt(sum_{k >= j} S_ik^2) */
    std::vector<char>   supp;          /* supp[a*n + i] = (b_a^i != 0), exact */
    std::vector<char>   inA;           /* the current support pattern U */
    std::vector<slong>  cnt;           /* trace bound on members per candidate */
    std::map<std::pair<slong, slong>, std::vector<std::vector<slong>>> rawpair[2];  /* [use_caps] */
    arb_t               gball;         /* enclosure of the field generator */

    /* lattices of n with (S n)_i = 0 on a set of labels, by that set */
    std::map<std::vector<char>, std::unique_ptr<dlattice>> lattices;

    std::vector<slong>  cap;           /* cap[a*n + j]: NIM-rep budget of candidate a */
    std::vector<double> dj;            /* quantum dimensions S_j0 / S_00 */
    std::vector<double> g;             /* g = vacuum coefficient of each candidate */
    std::vector<char>   pf_ok;         /* x^j <= d_j for all j (else in no complete family) */
    slong               lp_prunes, lp_rows, lp_calls;
    int                 use_caps;

    std::vector<family_t> found;
    slong      checked, accepted, refinements, nodes, open_subset;
    int        truncated;
    slong      patterns_done, npatterns;
    slong      root_at, root_of;     /* first-member branch being explored, of how many */
    double     t0, t_last;
};

static double dis_now()
{
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

static void dis_progress(dis_t& D, int force)
{
    double t = dis_now();
    if (!force && t - D.t_last < PROGRESS_SECONDS) return;
    D.t_last = t;
    printf("  [%6.1fs] support patterns %ld/%ld (first member %ld/%ld), refinements %ld, "
           "overlap checks %ld (%ld accepted), search nodes %ld, families %ld\n",
           t - D.t0, (long) D.patterns_done, (long) D.npatterns, (long) D.root_at, (long) D.root_of,
           (long) D.refinements,
           (long) D.checked, (long) D.accepted, (long) D.nodes, (long) D.found.size());
    fflush(stdout);
}

/* out = S n, exactly */
static void dis_Sn(nf_vec& out, dis_t& D, const std::vector<slong>& nv)
{
    search_t* e = D.e;
    nf_elem_t t; nf_elem_init(t, e->field);
    for (slong i = 0; i < D.n; i++) {
        nf_elem_zero(out[i], e->field);
        for (slong j = 0; j < D.n; j++) {
            if (!nv[j]) continue;
            nf_elem_scalar_mul_si(t, entry(e->S, i, j, e->ctx), nv[j], e->field);
            nf_elem_add(out[i], out[i], t, e->field);
        }
    }
    nf_elem_clear(t, e->field);
}

static int dis_sign(const nf_elem_t a, dis_t& D)
{
    arb_t v; arb_init(v);
    int s = nf_elem_sign(v, a, D.e->gen, D.e->field);
    arb_clear(v);
    return s;
}

static double dis_d(const nf_elem_t a, dis_t& D)
{
    arb_t v; arb_init(v);
    nf_elem_get_arb(v, a, D.gball, D.e->field, NUM_PREC);
    double d = arf_get_d(arb_midref(v), ARF_RND_NEAR);
    arb_clear(v);
    return d;
}

/* Number of nonnegative integer vectors of length n with sum of squares <= c. */
static double count_sq(slong c, slong n)
{
    std::vector<double> w(c + 1, 0.0), w2(c + 1);
    w[0] = 1.0;
    for (slong j = 0; j < n; j++) {
        std::fill(w2.begin(), w2.end(), 0.0);
        for (slong s = 0; s <= c; s++)
            for (slong v = 0; v * v <= s; v++) w2[s] += w[s - v * v];
        w.swap(w2);
    }
    double t = 0.0;
    for (slong s = 0; s <= c; s++) t += w[s];
    return t;
}

/* Allowed values of p_i = <s^i, v^i>: kind 0 = cen, 1 = cen +- rad, 2 = [cen - rad, cen + rad]. */
struct target_t {
    std::vector<int>    kind;
    std::vector<double> cen, rad;
    std::vector<slong>  nmax;       /* n^j <= nmax[j] */
    std::vector<char>   zero;       /* p_i = 0 structurally (a state has no component there) */
};

static double dis_gap(const target_t& T, slong i, double p)
{
    if (T.kind[i] == 0) return std::fabs(p - T.cen[i]);
    if (T.kind[i] == 1) return std::min(std::fabs(p - T.cen[i] - T.rad[i]), std::fabs(p - T.cen[i] + T.rad[i]));
    double lo = T.cen[i] - T.rad[i], hi = T.cen[i] + T.rad[i];
    return p < lo ? lo - p : (p > hi ? p - hi : 0.0);
}

static dlattice* dis_lattice(dis_t& D, const std::vector<char>& zero)
{
    auto it = D.lattices.find(zero);
    if (it != D.lattices.end()) return it->second.get();
    search_t* e = D.e;
    const slong n = D.n, degree = fmpq_poly_degree(e->field->pol);
    slong extra = 0;
    for (slong i = 0; i < n; i++) if (zero[i] && e->inI[i]) extra++;
    const slong rz = (slong) e->piv.size();
    const slong rows = FLINT_MAX(rz + extra * degree, 1);

    fmpq_mat_t A, R; fmpq_mat_init(A, rows, n); fmpq_mat_init(R, rows, n);
    fmpz_t one; fmpz_init_set_ui(one, 1);
    for (slong r = 0; r < rz; r++)                     /* the search's equations (i not in I) */
        for (slong j = 0; j < n; j++)
            fmpq_set_fmpz_frac(fmpq_mat_entry(A, r, j), fmpz_mat_entry(e->Z, r, j), one);
    fmpq_poly_t p; fmpq_poly_init(p);
    for (slong i = 0, r = rz; i < n; i++) {
        if (!zero[i] || !e->inI[i]) continue;
        for (slong j = 0; j < n; j++) {
            nf_elem_get_fmpq_poly(p, entry(e->S, i, j, e->ctx), e->field);
            for (slong l = 0; l < degree; l++)
                fmpq_poly_get_coeff_fmpq(fmpq_mat_entry(A, r + l, j), p, l);
        }
        r += degree;
    }
    fmpq_poly_clear(p);
    fmpz_clear(one);
    slong rank = fmpq_mat_rref(R, A);

    std::unique_ptr<dlattice> L(new dlattice);
    L->Z = (fmpz_mat_struct*) flint_malloc(sizeof(fmpz_mat_struct));
    fmpz_mat_init(L->Z, FLINT_MAX(rank, 1), n);
    fmpq_mat_t Wn; fmpq_mat_window_init(Wn, R, 0, 0, FLINT_MAX(rank, 1), n);
    fmpq_mat_get_fmpz_mat_rowwise(L->Z, NULL, Wn);
    fmpq_mat_window_clear(Wn);
    fmpq_mat_clear(R); fmpq_mat_clear(A);

    std::vector<char> isp(n, 0);
    for (slong r = 0; r < rank; r++) {
        slong j = 0;
        while (j < n && fmpz_is_zero(fmpz_mat_entry(L->Z, r, j))) j++;
        L->piv.push_back(j);
        isp[j] = 1;
    }
    for (slong j = 0; j < n; j++) if (!isp[j]) L->fcol.push_back(j);
    const slong k = (slong) L->fcol.size();
    L->PR.assign(rank, std::vector<double>(k));
    for (slong r = 0; r < rank; r++) {
        double dz = fmpz_get_d(fmpz_mat_entry(L->Z, r, L->piv[r]));
        for (slong t = 0; t < k; t++)
            L->PR[r][t] = -fmpz_get_d(fmpz_mat_entry(L->Z, r, L->fcol[t])) / dz;
    }
    L->W.assign(n, std::vector<double>(k, 0.0));
    for (slong i = 0; i < n; i++)
        for (slong t = 0; t < k; t++) {
            double w = D.Sd[i * n + L->fcol[t]];
            for (slong r = 0; r < rank; r++) w += D.Sd[i * n + L->piv[r]] * L->PR[r][t];
            L->W[i][t] = w;
        }
    dlattice* out = L.get();
    D.lattices[zero] = std::move(L);
    return out;
}

/* Numeric search for n_sv: n >= 0, n^0 = 0, |n|^2 <= c, on the lattice Z n = 0
   (which encodes (S n)_i = 0 for i not in I exactly), with each (S n)_i inside
   the hull of its allowed values. Depth-first over the free coordinates; every
   linear row bounds the next coordinate given the rest's range (interval
   propagation). Pivot coordinates are recovered exactly at the leaves. */
struct drow { const double* w; double lo, hi; std::vector<double> rmin, rmax, p2, n2; std::vector<double> wo; };

static void dis_lat_dfs(dis_t& D, const dlattice& L, const target_t& T, slong c, std::vector<drow>& rows,
                        const std::vector<slong>& box, const std::vector<slong>& ord, slong d,
                        std::vector<double>& part, std::vector<slong>& f,
                        std::vector<std::vector<slong>>& out, slong sq = 0)
{
    const slong n = D.n, k = (slong) L.fcol.size();
    if (d == k) {
        std::vector<slong> nv(n, 0);
        for (slong t = 0; t < k; t++) nv[L.fcol[ord[t]]] = f[t];
        fmpz_t acc; fmpz_init(acc);
        int ok = 1;
        for (size_t r = 0; r < L.piv.size() && ok; r++) {         /* as in leaf() */
            slong p = L.piv[r];
            fmpz_zero(acc);
            for (slong j = 0; j < n; j++)
                if (j != p && nv[j] != 0) fmpz_addmul_si(acc, fmpz_mat_entry(L.Z, r, j), nv[j]);
            const fmpz* dz = fmpz_mat_entry(L.Z, r, p);
            if (!fmpz_divisible(acc, dz)) { ok = 0; break; }
            fmpz_divexact(acc, acc, dz);
            fmpz_neg(acc, acc);
            if (!fmpz_fits_si(acc) || fmpz_sgn(acc) < 0) { ok = 0; break; }
            nv[p] = fmpz_get_si(acc);
        }
        fmpz_clear(acc);
        if (!ok || nv[0] != 0) return;
        for (slong j = 0; j < n; j++) if (nv[j] > T.nmax[j]) return;
        slong nn = 0;
        for (slong j = 0; j < n; j++) nn += nv[j] * nv[j];
        if (nn > c) return;
        for (slong i = 0; i < n; i++) {                           /* allowed values, numerically */
            double p = 0.0;
            for (slong j = 0; j < n; j++) if (nv[j]) p += D.Sd[i * n + j] * nv[j];
            if (dis_gap(T, i, p) > DIS_TOL) return;
        }
        out.push_back(nv);
        return;
    }
    /* range of f_d allowed by every row, given the remaining coordinates: each in its
       box, all >= 0, and sum of squares <= c - sq (Cauchy-Schwarz on the row) */
    double lo = 0.0, hi = (double) box[d];
    const double budget = (double) (c - sq);
    for (size_t r = 0; r < rows.size(); r++) {
        const drow& R = rows[r];
        double w = R.w[d];
        double rmx = std::min(R.rmax[d + 1], std::sqrt(budget * R.p2[d + 1]) + 1e-12);
        double rmn = std::max(R.rmin[d + 1], -std::sqrt(budget * R.n2[d + 1]) - 1e-12);
        double a = R.lo - part[r] - rmx;
        double b = R.hi - part[r] - rmn;
        if (std::fabs(w) < 1e-12) {
            if (a > DIS_TOL || b < -DIS_TOL) return;
            continue;
        }
        double u = a / w, v = b / w;
        if (w < 0) std::swap(u, v);
        lo = std::max(lo, u); hi = std::min(hi, v);
        if (lo > hi + 1e-9) return;
    }
    slong a = (slong) std::ceil(lo - 1e-9), b = (slong) std::floor(hi + 1e-9);
    while (b >= a && sq + b * b > c) b--;                        /* |n|^2 <= c */
    for (slong val = a; val <= b; val++) {
        f[d] = val;
        for (size_t r = 0; r < rows.size(); r++) part[r] += rows[r].w[d] * val;
        dis_lat_dfs(D, L, T, c, rows, box, ord, d + 1, part, f, out, sq + val * val);
        for (size_t r = 0; r < rows.size(); r++) part[r] -= rows[r].w[d] * val;
    }
    f[d] = 0;
}

static void dis_find_lattice(dis_t& D, const dlattice& L, const target_t& T, slong c,
                             std::vector<std::vector<slong>>& out)
{
    search_t* e = D.e;
    const slong n = D.n, k = (slong) L.fcol.size();
    double pmax = 0.0;
    for (slong i = 0; i < n; i++) {
        double m = std::fabs(T.cen[i]) + (T.kind[i] ? T.rad[i] : 0.0);
        pmax += m * m;
    }
    c = FLINT_MIN(c, (slong) std::floor(pmax + 1e-6));        /* |n|^2 = |p|^2 */
    if (c < 0) return;
    const slong B = (slong) std::floor(std::sqrt((double) c) + 1e-9);

    std::vector<slong> box(k);
    for (slong t = 0; t < k; t++) box[t] = (L.fcol[t] == 0) ? 0 : FLINT_MIN(B, T.nmax[L.fcol[t]]);   /* n^0 = 0 */

    std::vector<drow> rows;
    for (slong i = 0; i < n; i++) {
        if (!e->inI[i] || T.zero[i]) continue;                    /* exact via the lattice */
        drow R; R.w = L.W[i].data();
        R.lo = T.cen[i] - (T.kind[i] ? T.rad[i] : 0.0) - DIS_TOL;
        R.hi = T.cen[i] + (T.kind[i] ? T.rad[i] : 0.0) + DIS_TOL;
        rows.push_back(R);
    }
    for (size_t r = 0; r < L.piv.size(); r++) {
        drow R; R.w = L.PR[r].data();
        R.lo = -DIS_TOL;
        R.hi = (L.piv[r] == 0) ? DIS_TOL : FLINT_MIN(B, T.nmax[L.piv[r]]) + DIS_TOL;
        rows.push_back(R);
    }
    /* fix the most constrained coordinates first: weight relative to row width */
    std::vector<double> score(k, 0.0);
    for (const drow& R : rows)
        for (slong t = 0; t < k; t++)
            score[t] += std::fabs(R.w[t]) * box[t] / (R.hi - R.lo + 1e-3);
    std::vector<slong> ord(k);
    for (slong t = 0; t < k; t++) ord[t] = t;
    std::stable_sort(ord.begin(), ord.end(), [&](slong a, slong b) { return score[a] > score[b]; });
    std::vector<slong> boxo(k);
    for (slong t = 0; t < k; t++) boxo[t] = box[ord[t]];
    for (drow& R : rows) {
        R.wo.resize(k);
        for (slong t = 0; t < k; t++) R.wo[t] = R.w[ord[t]];
        R.w = R.wo.data();
    }
    for (drow& R : rows) {                                        /* rest ranges over boxes */
        R.rmin.assign(k + 1, 0.0); R.rmax.assign(k + 1, 0.0);
        R.p2.assign(k + 1, 0.0); R.n2.assign(k + 1, 0.0);
        for (slong t = k - 1; t >= 0; t--) {
            double x = R.w[t] * boxo[t];
            R.rmin[t] = R.rmin[t + 1] + std::min(0.0, x);
            R.rmax[t] = R.rmax[t + 1] + std::max(0.0, x);
            double w = boxo[t] > 0 ? R.w[t] : 0.0;
            R.p2[t] = R.p2[t + 1] + (w > 0 ? w * w : 0.0);
            R.n2[t] = R.n2[t + 1] + (w < 0 ? w * w : 0.0);
        }
    }
    std::vector<double> part(rows.size(), 0.0);
    std::vector<slong> f(k, 0);
    dis_lat_dfs(D, L, T, c, rows, boxo, ord, 0, part, f, out);
}

/* Numeric search 1 (no interval labels): choose the point for every two-point label. */
static void dis_label_dfs(dis_t& D, const std::vector<slong>& Lb, const std::vector<double>& rr,
                          const std::vector<double>& tail, const std::vector<slong>& nmax,
                          size_t d, std::vector<double>& v, std::vector<std::vector<slong>>& out)
{
    const slong n = D.n;
    if (d == Lb.size()) {
        std::vector<slong> nv(n);
        for (slong j = 0; j < n; j++) {
            double k = std::nearbyint(v[j]);
            if (std::fabs(v[j] - k) > DIS_TOL || k < 0) return;
            nv[j] = (slong) k;
        }
        if (nv[0] != 0) return;                        /* distinct boundaries */
        for (slong j = 0; j < n; j++) if (nv[j] > nmax[j]) return;
        out.push_back(nv);
        return;
    }
    for (slong j = 0; j < n; j++) {
        double lo = v[j] - tail[d * n + j], hi = v[j] + tail[d * n + j];
        if (hi < -DIS_TOL) return;
        if (std::floor(hi + DIS_TOL) < std::ceil(lo - DIS_TOL)) return;
        if (j == 0 && (lo > DIS_TOL || hi < -DIS_TOL)) return;
        if (lo > nmax[j] + DIS_TOL) return;
    }
    const slong i = Lb[d];
    for (int eps = 1; eps >= -1; eps -= 2) {
        for (slong j = 0; j < n; j++) v[j] += eps * rr[d] * D.Sd[i * n + j];
        dis_label_dfs(D, Lb, rr, tail, nmax, d + 1, v, out);
        for (slong j = 0; j < n; j++) v[j] -= eps * rr[d] * D.Sd[i * n + j];
    }
}

/* Numeric search 2: n >= 0 with n^0 = 0 and |n|^2 <= c; coordinates not yet
   chosen move p_i by at most sqrt(budget) * Qs_i (Cauchy-Schwarz). */
static void dis_norm_dfs(dis_t& D, const target_t& T, slong j, slong budget,
                         std::vector<slong>& nv, std::vector<double>& p,
                         std::vector<std::vector<slong>>& out)
{
    const slong n = D.n;
    if (budget > 0 && j < n) {
        const double sb = std::sqrt((double) budget);
        for (slong i = 0; i < n; i++)
            if (dis_gap(T, i, p[i]) > sb * D.Qs[j * n + i] + DIS_TOL) return;
    }
    if (j == n || budget == 0) {
        for (slong i = 0; i < n; i++) if (dis_gap(T, i, p[i]) > DIS_TOL) return;
        std::vector<slong> full(nv);
        std::fill(full.begin() + j, full.end(), 0);
        out.push_back(full);
        return;
    }
    for (slong val = 0; val * val <= budget && val <= T.nmax[j]; val++) {
        nv[j] = val;
        if (val) for (slong i = 0; i < n; i++) p[i] += val * D.Sd[i * n + j];
        dis_norm_dfs(D, T, j + 1, budget - val * val, nv, p, out);
        if (val) for (slong i = 0; i < n; i++) p[i] -= val * D.Sd[i * n + j];
    }
    nv[j] = 0;
}

/* Numeric search 4 (few interval labels): n = S p with p fixed on pinned labels,
   +-rad on two-point labels and p_i = t_i free in an interval on the f interval
   labels. Since the columns S_.i are orthonormal, n runs over an f-dimensional
   affine set: enumerate the integers on f well-conditioned coordinates, solve
   for t, and check the rest. Returns 0 (doing nothing) if that would be too many. */
#define AFFINE_MAX_FREE   6
#define AFFINE_MAX_SIGNS 12
#define AFFINE_MAX_POINTS 2000000

static void dis_affine_rec(dis_t& D, const target_t& T, slong c, const std::vector<slong>& Fi,
                           const std::vector<slong>& J, const std::vector<double>& Minv,
                           const std::vector<double>& base, const std::vector<slong>& lo,
                           const std::vector<slong>& hi, size_t d, std::vector<slong>& nJ,
                           std::vector<std::vector<slong>>& out)
{
    const slong n = D.n, f = (slong) Fi.size();
    if (d < J.size()) {
        for (slong v = lo[d]; v <= hi[d]; v++) {
            nJ[d] = v;
            dis_affine_rec(D, T, c, Fi, J, Minv, base, lo, hi, d + 1, nJ, out);
        }
        return;
    }
    std::vector<double> t(f, 0.0);
    for (slong a = 0; a < f; a++)
        for (slong b = 0; b < f; b++) t[a] += Minv[a * f + b] * (nJ[b] - base[J[b]]);
    for (slong a = 0; a < f; a++) {
        slong i = Fi[a];
        if (t[a] < T.cen[i] - T.rad[i] - DIS_TOL || t[a] > T.cen[i] + T.rad[i] + DIS_TOL) return;
    }
    std::vector<slong> nv(n);
    slong nn = 0;
    for (slong j = 0; j < n; j++) {
        double x = base[j];
        for (slong a = 0; a < f; a++) x += D.Sd[Fi[a] * n + j] * t[a];
        double k = std::nearbyint(x);
        if (std::fabs(x - k) > DIS_TOL || k < 0 || k > T.nmax[j]) return;
        nv[j] = (slong) k;
        nn += nv[j] * nv[j];
    }
    if (nv[0] != 0 || nn > c) return;
    out.push_back(nv);
}

static int dis_find_affine(dis_t& D, const target_t& T, slong c, std::vector<std::vector<slong>>& out)
{
    const slong n = D.n;
    std::vector<slong> Fi, Tw;
    for (slong i = 0; i < n; i++) {
        if (T.kind[i] == 2) Fi.push_back(i);
        if (T.kind[i] == 1) Tw.push_back(i);
    }
    const slong f = (slong) Fi.size();
    if (f == 0 || f > AFFINE_MAX_FREE || (slong) Tw.size() > AFFINE_MAX_SIGNS) return 0;
    const slong B = (slong) std::floor(std::sqrt((double) c) + 1e-9);

    /* f coordinates of n by greedy elimination on A = (S_ji), j rows, i in Fi */
    std::vector<double> A(n * f);
    for (slong j = 0; j < n; j++) for (slong a = 0; a < f; a++) A[j * f + a] = D.Sd[Fi[a] * n + j];
    std::vector<slong> J;
    std::vector<char> used(n, 0);
    std::vector<double> W(A);
    for (slong a = 0; a < f; a++) {
        slong best = -1;
        for (slong j = 0; j < n; j++)
            if (!used[j] && (best < 0 || std::fabs(W[j * f + a]) > std::fabs(W[best * f + a]))) best = j;
        if (best < 0 || std::fabs(W[best * f + a]) < 1e-9) return 0;
        used[best] = 1; J.push_back(best);
        for (slong j = 0; j < n; j++) {
            if (j == best) continue;
            double r = W[j * f + a] / W[best * f + a];
            for (slong b = a; b < f; b++) W[j * f + b] -= r * W[best * f + b];
        }
    }
    /* Minv = (A restricted to rows J)^-1 */
    std::vector<double> G(f * 2 * f, 0.0);
    for (slong a = 0; a < f; a++) {
        for (slong b = 0; b < f; b++) G[a * 2 * f + b] = A[J[a] * f + b];
        G[a * 2 * f + f + a] = 1.0;
    }
    for (slong col = 0; col < f; col++) {
        slong pr = col;
        for (slong r = col + 1; r < f; r++)
            if (std::fabs(G[r * 2 * f + col]) > std::fabs(G[pr * 2 * f + col])) pr = r;
        if (std::fabs(G[pr * 2 * f + col]) < 1e-12) return 0;
        for (slong k = 0; k < 2 * f; k++) std::swap(G[col * 2 * f + k], G[pr * 2 * f + k]);
        double pv = G[col * 2 * f + col];
        for (slong k = 0; k < 2 * f; k++) G[col * 2 * f + k] /= pv;
        for (slong r = 0; r < f; r++) {
            if (r == col) continue;
            double m = G[r * 2 * f + col];
            for (slong k = 0; k < 2 * f; k++) G[r * 2 * f + k] -= m * G[col * 2 * f + k];
        }
    }
    std::vector<double> Minv(f * f);
    for (slong a = 0; a < f; a++) for (slong b = 0; b < f; b++) Minv[a * f + b] = G[a * 2 * f + f + b];

    for (ulong mask = 0; mask < (1ul << Tw.size()); mask++) {
        std::vector<double> base(n, 0.0);
        for (slong i = 0; i < n; i++) {
            if (T.kind[i] == 2) continue;
            double p = T.cen[i];
            if (T.kind[i] == 1) {
                size_t k = std::find(Tw.begin(), Tw.end(), i) - Tw.begin();
                p += ((mask >> k) & 1) ? -T.rad[i] : T.rad[i];
            }
            if (p != 0.0) for (slong j = 0; j < n; j++) base[j] += p * D.Sd[i * n + j];
        }
        std::vector<slong> lo(f), hi(f);
        double points = 1.0;
        int empty = 0;
        for (slong b = 0; b < f; b++) {
            slong j = J[b];
            double l = base[j], h = base[j];
            for (slong a = 0; a < f; a++) {
                double w = A[j * f + a], tl = T.cen[Fi[a]] - T.rad[Fi[a]], th = T.cen[Fi[a]] + T.rad[Fi[a]];
                l += std::min(w * tl, w * th);
                h += std::max(w * tl, w * th);
            }
            lo[b] = FLINT_MAX(0, (slong) std::ceil(l - DIS_TOL));
            hi[b] = FLINT_MIN(j == 0 ? 0 : FLINT_MIN(B, T.nmax[j]), (slong) std::floor(h + DIS_TOL));
            if (lo[b] > hi[b]) empty = 1;
            points *= (double) (hi[b] - lo[b] + 1);
        }
        if (empty) continue;
        if (points > AFFINE_MAX_POINTS) { out.clear(); return 0; }
        std::vector<slong> nJ(f);
        dis_affine_rec(D, T, c, Fi, J, Minv, base, lo, hi, 0, nJ, out);
    }
    return 1;
}

static void dis_find(dis_t& D, const target_t& T, slong c, std::vector<std::vector<slong>>& out)
{
    const slong n = D.n;
    int interval = 0;
    std::vector<slong> Lb;
    double pmax = 0.0;
    for (slong i = 0; i < n; i++) {
        if (T.kind[i] == 2) interval = 1;
        if (T.kind[i] == 1) Lb.push_back(i);
        double m = std::fabs(T.cen[i]) + (T.kind[i] ? T.rad[i] : 0.0);
        pmax += m * m;
    }
    c = FLINT_MIN(c, (slong) std::floor(pmax + 1e-6));    /* |n|^2 = |p|^2 */
    if (c < 0) return;
    /* With labels whose overlap is only bounded (several Ishibashi states), the
       lattice of the search's equations is the better space when it is much
       smaller than all of n; otherwise enumerate n by norm, or the two-point
       labels by sign, whichever is fewer. */
    if (interval && dis_find_affine(D, T, c, out)) return;
    if (interval) {
        dlattice* L = dis_lattice(D, T.zero);
        dis_find_lattice(D, *L, T, c, out); return;
    }
    if (interval || count_sq(c, n - 1) < std::ldexp(1.0, (int) FLINT_MIN((slong) Lb.size(), 1000))) {
        std::vector<slong> nv(n, 0);
        std::vector<double> p(n, 0.0);
        dis_norm_dfs(D, T, 1, c, nv, p, out);                 /* starts at j = 1: n^0 = 0 */
    } else {
        std::sort(Lb.begin(), Lb.end(), [&](slong u, slong w) { return T.rad[u] > T.rad[w]; });
        std::vector<double> rr(Lb.size()), tail((Lb.size() + 1) * n, 0.0), v(n, 0.0);
        for (slong i = 0; i < n; i++)
            if (T.cen[i] != 0.0) for (slong j = 0; j < n; j++) v[j] += T.cen[i] * D.Sd[i * n + j];
        for (size_t t = 0; t < Lb.size(); t++) rr[t] = T.rad[Lb[t]];
        for (slong d = (slong) Lb.size() - 1; d >= 0; d--)
            for (slong j = 0; j < n; j++)
                tail[d * n + j] = tail[(d + 1) * n + j] + std::fabs(D.Sd[Lb[d] * n + j]) * rr[d];
        dis_label_dfs(D, Lb, rr, tail, T.nmax, 0, v, out);
    }
}

/* the state of candidate B before any member: no coordinates, r = |b_B|^2 */
static dstate dis_raw(dis_t& D, slong B)
{
    const slong n = D.n;
    dstate s;
    s.cand = B;
    s.ex = nf_vec(D.exact ? D.L : 0, D.e->field);
    s.nd.assign(D.L, 0.0);
    s.rz.assign(n, 1);
    for (slong i = 0; i < n; i++) {
        slong ri = D.off[i] + D.mult[i];
        if (D.exact) nf_elem_set(s.ex[ri], D.q + B * n + i, D.e->field);
        s.rz[i] = !D.supp[B * n + i];
        if (!s.rz[i]) s.nd[ri] = D.qd[B * n + i];
    }
    s.bound = 0;
    s.used.assign(n, 0);
    s.pf.assign(n, 0.0);
    return s;
}

/* Does v raise the rank at label i (given frame F)? */
static int dis_pivot(dis_t& D, const dframe& F, const dstate& v, slong i)
{
    return D.inA[i] && F.c[i] < D.mult[i] && !v.rz[i];
}

/* Refine class s against the new member v (frame F before v, Fn after):
   every verified overlap n_sv gives one refined class, appended to out. */
static void dis_refine(dis_t& D, const dframe& F, const dframe& Fn, const dstate& s,
                       const dstate& v, std::vector<dstate>& out,
                       const std::vector<slong>* forced = nullptr)
{
    search_t* e = D.e;
    const slong n = D.n;
    D.refinements++;
    dis_progress(D, 0);

    /* targets */
    target_t T;
    T.kind.assign(n, 0); T.cen.assign(n, 0.0); T.rad.assign(n, 0.0);
    T.nmax.assign(n, WORD_MAX);
    T.zero.assign(n, 0);
    for (slong i = 0; i < n; i++)
        T.zero[i] = !e->inI[i] || !D.inA[i] || !D.supp[s.cand * n + i] || !D.supp[v.cand * n + i];
    if (D.use_caps) {
        const double gs = D.g[s.cand], gv = D.g[v.cand];
        for (slong j = 0; j < n; j++) {
            slong left = FLINT_MIN(D.cap[s.cand * n + j] - s.used[j], D.cap[v.cand * n + j] - v.used[j]);
            if (left < 0) return;
            slong m = (slong) std::floor(std::sqrt((double) left) + 1e-9);
            while (m * m > left) m--;
            while ((m + 1) * (m + 1) <= left) m++;
            /* Perron-Frobenius: sum_l n^j_kl g_l = d_j g_k, all terms >= 0 */
            double ss = (D.dj[j] - e->sols[s.cand][j]) * gs - s.pf[j];
            double sv = (D.dj[j] - e->sols[v.cand][j]) * gv - v.pf[j];
            double bnd = std::min(ss / gv, sv / gs);
            if (bnd < -1e-7) return;
            m = FLINT_MIN(m, (slong) std::floor(bnd + 1e-7));
            T.nmax[j] = m;
        }
    }
    for (slong i = 0; i < n; i++) {
        if (!D.inA[i]) continue;                                   /* p_i = 0 */
        const slong o = D.off[i], m = D.mult[i];
        double pr = 0.0;
        for (slong j = 0; j < F.c[i]; j++) pr += s.nd[o + j] * v.nd[o + j] / F.dd[o + j];
        T.cen[i] = pr;
        if (F.c[i] == m || s.rz[i] || v.rz[i]) continue;
        T.rad[i] = std::sqrt(std::max(0.0, s.nd[o + m] * v.nd[o + m]));
        T.kind[i] = (m - F.c[i] == 1) ? 1 : 2;
    }
    slong c = 0;
    for (slong j = 0; j < n; j++) c += e->sols[s.cand][j] * e->sols[v.cand][j];
    /* two candidates before any other member: the answer depends on them alone */
    const int raw = s.nh.empty() && v.nh.empty();
    const std::pair<slong, slong> key(FLINT_MIN(s.cand, v.cand), FLINT_MAX(s.cand, v.cand));
    std::vector<std::vector<slong>> cands;
    if (forced) cands.push_back(*forced);
    else if (raw && D.rawpair[D.use_caps].count(key)) cands = D.rawpair[D.use_caps][key];
    else {
        dis_find(D, T, c, cands);
        if (raw) D.rawpair[D.use_caps][key] = cands;
    }
    if (cands.empty()) return;

    if (!D.exact) {
        /* numeric acceptance and update; integer conditions exact */
        fmpz_t acc; fmpz_init(acc);
        std::vector<double> P(n), proj(n, 0.0), rr(n, 0.0);
        for (slong i = 0; i < n; i++) {
            if (!D.inA[i]) continue;
            const slong o = D.off[i], m = D.mult[i];
            for (slong j = 0; j < F.c[i]; j++) proj[i] += s.nd[o + j] * v.nd[o + j] / F.dd[o + j];
            rr[i] = s.nd[o + m] * v.nd[o + m];
        }
        for (const std::vector<slong>& nv : cands) {
            D.checked++;
            int ok = 1;
            for (slong r = 0; r < (slong) e->piv.size() && ok; r++) {
                fmpz_zero(acc);
                for (slong j = 0; j < n; j++)
                    if (nv[j]) fmpz_addmul_si(acc, fmpz_mat_entry(e->Z, r, j), nv[j]);
                ok = fmpz_is_zero(acc);
            }
            if (!ok) continue;
            for (slong i = 0; i < n && ok; i++) {
                if (!e->inI[i]) continue;
                double p = 0.0;
                for (slong j = 0; j < n; j++) if (nv[j]) p += D.Sd[i * n + j] * nv[j];
                P[i] = p;
                if (!D.inA[i]) { ok = std::fabs(p) < DIS_TOL; continue; }
                double dl = p - proj[i];
                if (T.kind[i] == 0) ok = std::fabs(dl) < DIS_TOL;
                else if (T.kind[i] == 1) ok = std::fabs(dl * dl - rr[i]) < DIS_TOL * (1.0 + rr[i]);
                else ok = dl * dl <= rr[i] + DIS_TOL * (1.0 + rr[i]);
            }
            if (!ok) continue;
            D.accepted++;

            dstate x = s.clone();
            x.nh.push_back(nv);
            for (slong j = 0; j < n; j++) { x.used[j] += nv[j] * nv[j]; x.pf[j] += nv[j] * D.g[v.cand]; }
            for (slong i = 0; i < n; i++) {
                if (!dis_pivot(D, F, v, i)) continue;
                const slong o = D.off[i], m = D.mult[i], k = F.c[i];
                double y = P[i] - proj[i];
                double r = x.nd[o + m] - y * y / v.nd[o + m];
                x.nd[o + k] = y;
                x.rz[i] = (r <= ZERO_TOL * (1.0 + D.qd[x.cand * n + i]));
                x.nd[o + m] = x.rz[i] ? 0.0 : r;
            }
            x.bound = 1;
            for (slong i = 0; i < n; i++)
                if (D.inA[i] && !x.rz[i]) x.bound += D.mult[i] - Fn.c[i];
            x.bound = FLINT_MIN(x.bound, D.cnt[x.cand]);
            out.push_back(std::move(x));
        }
        fmpz_clear(acc);
        return;
    }

    /* exact projections once per pair */
    nf_vec proj(n, e->field), rr(n, e->field), P(n, e->field);
    nf_elem_t t, u; nf_elem_init(t, e->field); nf_elem_init(u, e->field);
    for (slong i = 0; i < n; i++) {
        nf_elem_zero(proj[i], e->field);
        if (!D.inA[i]) continue;
        const slong o = D.off[i], m = D.mult[i];
        for (slong j = 0; j < F.c[i]; j++) {
            nf_elem_mul(t, s.ex[o + j], v.ex[o + j], e->field);
            nf_elem_mul(t, t, F.dinv[o + j], e->field);
            nf_elem_add(proj[i], proj[i], t, e->field);
        }
        nf_elem_mul(rr[i], s.ex[o + m], v.ex[o + m], e->field);
    }

    fmpz_t acc; fmpz_init(acc);
    for (const std::vector<slong>& nv : cands) {
        D.checked++;
        int ok = 1;
        /* (S n)_i = 0 for i not in I  <=>  Z n = 0, in integers */
        for (slong r = 0; r < (slong) e->piv.size() && ok; r++) {
            fmpz_zero(acc);
            for (slong j = 0; j < n; j++)
                if (nv[j]) fmpz_addmul_si(acc, fmpz_mat_entry(e->Z, r, j), nv[j]);
            ok = fmpz_is_zero(acc);
        }
        if (!ok) continue;
        for (slong i = 0; i < n; i++) {                           /* (S n)_i for i in I only */
            nf_elem_zero(P[i], e->field);
            if (!e->inI[i]) continue;
            for (slong j = 0; j < n; j++) {
                if (!nv[j]) continue;
                nf_elem_scalar_mul_si(t, entry(e->S, i, j, e->ctx), nv[j], e->field);
                nf_elem_add(P[i], P[i], t, e->field);
            }
        }
        for (slong i = 0; i < n && ok; i++) {
            if (!e->inI[i]) continue;
            if (!D.inA[i]) { ok = nf_elem_is_zero(P[i], e->field); continue; }
            nf_elem_sub(t, P[i], proj[i], e->field);               /* Delta = P - proj */
            if (T.kind[i] == 0) { ok = nf_elem_is_zero(t, e->field); continue; }
            nf_elem_mul(u, t, t, e->field);
            if (T.kind[i] == 1) ok = nf_elem_equal(u, rr[i], e->field);
            else { nf_elem_sub(u, rr[i], u, e->field); ok = (dis_sign(u, D) >= 0); }
        }
        if (!ok) continue;
        D.accepted++;

        dstate x = s.clone();
        x.nh.push_back(nv);
        for (slong j = 0; j < n; j++) x.used[j] += nv[j] * nv[j];
        for (slong i = 0; i < n; i++) {
            if (!dis_pivot(D, F, v, i)) continue;
            const slong o = D.off[i], m = D.mult[i], k = F.c[i];
            nf_elem_sub(x.ex[o + k], P[i], proj[i], e->field);           /* new coordinate */
            nf_elem_mul(t, x.ex[o + k], x.ex[o + k], e->field);
            nf_elem_mul(t, t, Fn.dinv[o + k], e->field);                 /* 1 / r_v */
            nf_elem_sub(x.ex[o + m], x.ex[o + m], t, e->field);          /* new residual */
            x.nd[o + k] = dis_d(x.ex[o + k], D);
            x.rz[i] = nf_elem_is_zero(x.ex[o + m], e->field);
            x.nd[o + m] = x.rz[i] ? 0.0 : dis_d(x.ex[o + m], D);
        }
        x.bound = 1;
        for (slong i = 0; i < n; i++)
            if (D.inA[i] && !x.rz[i]) x.bound += D.mult[i] - Fn.c[i];
        x.bound = FLINT_MIN(x.bound, D.cnt[x.cand]);
        out.push_back(std::move(x));
    }
    fmpz_clear(acc);
    nf_elem_clear(u, e->field); nf_elem_clear(t, e->field);
}

/* frame after adding v */
static dframe dis_empty_frame(dis_t& D)
{
    dframe F;
    F.c.assign(D.n, 0);
    F.d = nf_vec(D.exact ? D.L : 0, D.e->field);
    F.dinv = nf_vec(D.exact ? D.L : 0, D.e->field);
    F.dd.assign(D.L, 0.0);
    return F;
}

static dframe dis_extend_frame(dis_t& D, const dframe& F, const dstate& v)
{
    dframe G = F.clone();
    for (slong i = 0; i < D.n; i++) {
        if (!dis_pivot(D, F, v, i)) continue;
        const slong o = D.off[i], m = D.mult[i];
        if (D.exact) {
            nf_elem_set(G.d[o + G.c[i]], v.ex[o + m], D.e->field);
            nf_elem_inv(G.dinv[o + G.c[i]], v.ex[o + m], D.e->field);
        }
        G.dd[o + G.c[i]] = v.nd[o + m];
        G.c[i]++;
    }
    return G;
}

static nmat_t dis_nmatrix(dis_t& D, const std::vector<member_t>& F)
{
    nmat_t N(F.size(), std::vector<std::vector<slong>>(F.size()));
    for (size_t k = 0; k < F.size(); k++) {
        N[k][k] = D.e->sols[F[k].cand];
        for (size_t l = 0; l < k; l++) { N[k][l] = F[k].nh[l]; N[l][k] = F[k].nh[l]; }
    }
    return N;
}

/* Same family up to gauge: a candidate-preserving bijection with equal n_kl. */
static int dis_iso_rec(const std::vector<member_t>& F, const nmat_t& NF,
                       const std::vector<member_t>& G, const nmat_t& NG,
                       std::vector<slong>& pi, std::vector<char>& used, size_t k)
{
    if (k == F.size()) return 1;
    for (size_t g = 0; g < G.size(); g++) {
        if (used[g] || G[g].cand != F[k].cand) continue;
        int ok = 1;
        for (size_t l = 0; l < k && ok; l++) ok = (NF[k][l] == NG[g][pi[l]]);
        if (!ok) continue;
        used[g] = 1; pi[k] = (slong) g;
        if (dis_iso_rec(F, NF, G, NG, pi, used, k + 1)) return 1;
        used[g] = 0;
    }
    return 0;
}

typedef std::shared_ptr<dstate> sp_t;    /* pool entries are immutable once made */

struct dsearch {
    slong                 T;
    std::vector<dframe>   frames;     /* frames[k]: the frame before member k */
    std::vector<std::vector<slong>> mused;   /* members' budgets used so far */
    std::vector<slong>    Xrem;       /* S w minus the members' x so far */
    std::vector<dstate>   R;          /* members, as states relative to the frame when added */
    std::vector<std::vector<slong>> c_at;
};

/* closed: no state of any candidate (any support) fits with every member */
static int dis_closed(dis_t& D, dsearch& S)
{
    const slong M = D.M;
    std::vector<char> saveA = D.inA;
    int saved_caps = D.use_caps;
    D.use_caps = 0;
    int closed = 1;
    for (slong X = 0; X < M && closed; X++) {
        D.inA = saveA;
        dframe F = dis_empty_frame(D);
        std::vector<dstate> cur;
        cur.push_back(dis_raw(D, X));
        for (size_t k = 0; k < S.R.size() && !cur.empty(); k++) {
            dframe Fn = dis_extend_frame(D, F, S.R[k]);
            std::vector<dstate> nxt;
            for (const dstate& s : cur) {
                /* X's labels outside U carry no overlap with the members */
                dis_refine(D, F, Fn, s, S.R[k], nxt);
            }
            cur.swap(nxt);
            F = std::move(Fn);
        }
        if (!cur.empty()) closed = 0;
    }
    D.inA = saveA;
    D.use_caps = saved_caps;
    return closed;
}

static void dis_record(dis_t& D, dsearch& S)
{
    family_t fam;
    fam.U = D.inA;
    slong N = 0;
    for (slong i = 0; i < D.n; i++) if (D.inA[i]) N += D.mult[i];
    fam.kind = (N == D.e->dimI) ? 2 : 1;
    for (size_t k = 0; k < S.R.size(); k++) {
        member_t m;
        m.cand = S.R[k].cand; m.nh = S.R[k].nh; m.c_at = S.c_at[k];
        fam.mem.push_back(std::move(m));
    }
    fam.N = dis_nmatrix(D, fam.mem);
    for (const family_t& g : D.found) {
        if (g.mem.size() != fam.mem.size()) continue;
        std::vector<slong> pi(fam.mem.size());
        std::vector<char> used(fam.mem.size(), 0);
        if (dis_iso_rec(fam.mem, fam.N, g.mem, g.N, pi, used, 0)) return;
    }
    if (fam.kind == 1 && !dis_closed(D, S)) { D.open_subset++; return; }

    /* exact replay: rebuild every member in the number field with its overlaps
       fixed; this checks all Gram conditions exactly */
    {
        int saved_caps = D.use_caps;
        D.exact = 1; D.use_caps = 0;
        std::vector<dframe> fr;
        fr.push_back(dis_empty_frame(D));
        std::vector<dstate> Rx;
        int ok = 1;
        for (size_t k = 0; k < fam.mem.size() && ok; k++) {
            std::vector<dstate> cur;
            cur.push_back(dis_raw(D, fam.mem[k].cand));
            for (size_t l = 0; l < k && ok; l++) {
                std::vector<dstate> nxt;
                dis_refine(D, fr[l], fr[l + 1], cur[0], Rx[l], nxt, &fam.mem[k].nh[l]);
                if (nxt.size() != 1) ok = 0;
                else cur.swap(nxt);
            }
            if (!ok) break;
            fam.mem[k].ex = cur[0].ex.clone();
            fam.mem[k].c_at = fr[k].c;
            fr.push_back(dis_extend_frame(D, fr[k], cur[0]));
            Rx.push_back(std::move(cur[0]));
        }
        D.exact = 0; D.use_caps = saved_caps;
        if (!ok) {
            D.rejected++;
            printf("  [%6.1fs] a %ld-member family failed the exact check and was discarded\n",
                   dis_now() - D.t0, (long) fam.mem.size());
            fflush(stdout);
            return;
        }
        fam.frame = fr.back().clone();
    }
    printf("  [%6.1fs] family %ld: %ld boundaries, %s\n", dis_now() - D.t0,
           (long) D.found.size() + 1, (long) fam.mem.size(),
           fam.kind == 2 ? "complete for I" : "complete and closed on a subset of I");
    fflush(stdout);
    D.found.push_back(std::move(fam));
}

static int dis_fits(dis_t& D, const std::vector<slong>& Xrem, slong cand)
{
    const std::vector<slong>& x = D.e->sols[cand];
    for (slong j = 0; j < D.n; j++) if (x[j] > Xrem[j]) return 0;
    return 1;
}

/* Knapsack: can the classes not dropped, each supplying at most its bound, sum
   exactly to Xrem? LP relaxation, and optionally the integer problem (with a time
   limit). Returns 0 only if infeasibility is established. */
static int dis_knapsack(dis_t& D, const std::vector<sp_t>& P, const std::vector<char>& dropped,
                        const std::vector<slong>& Xrem, int integer = 0)
{
    const slong n = D.n;
    std::map<slong, slong> cap;                    /* candidate -> total bound */
    for (size_t t = 0; t < P.size(); t++) if (!dropped[t]) cap[P[t]->cand] += P[t]->bound;
    for (slong j = 0; j < n; j++) if (Xrem[j] < 0) return 0;
    if (cap.empty()) { for (slong j = 0; j < n; j++) if (Xrem[j]) return 0; return 1; }
    glp_prob* lp = glp_create_prob();
    glp_add_rows(lp, (int) n);
    for (slong j = 0; j < n; j++) glp_set_row_bnds(lp, (int) j + 1, GLP_FX, (double) Xrem[j], (double) Xrem[j]);
    glp_add_cols(lp, (int) cap.size());
    std::vector<int> ia(1, 0), ja(1, 0);
    std::vector<double> ar(1, 0.0);
    int col = 0;
    for (const auto& kv : cap) {
        col++;
        glp_set_col_bnds(lp, col, GLP_DB, 0.0, (double) kv.second);
        const std::vector<slong>& x = D.e->sols[kv.first];
        for (slong j = 0; j < n; j++)
            if (x[j]) { ia.push_back((int) j + 1); ja.push_back(col); ar.push_back((double) x[j]); }
    }
    glp_load_matrix(lp, (int) ia.size() - 1, ia.data(), ja.data(), ar.data());
    glp_smcp parm; glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;
    parm.presolve = GLP_ON;
    int ret = glp_simplex(lp, &parm);
    int feasible = 1;
    if (ret == GLP_ENOPFS || (ret == 0 && glp_get_status(lp) == GLP_NOFEAS)) feasible = 0;
    if (feasible && integer) {                     /* integer knapsack, with a time limit */
        for (int c = 1; c <= col; c++) glp_set_col_kind(lp, c, GLP_IV);
        glp_iocp ip; glp_init_iocp(&ip);
        ip.msg_lev = GLP_MSG_OFF;
        ip.presolve = GLP_ON;
        ip.tm_lim = 10000;                         /* ms; inconclusive counts as feasible */
        int r2 = glp_intopt(lp, &ip);
        if (r2 == GLP_ENOPFS || (r2 == 0 && glp_mip_status(lp) == GLP_NOFEAS)) feasible = 0;
    }
    glp_delete_prob(lp);
    return feasible;
}

/* Bring a pool state up to date: refine it against the members it has not seen
   (it has seen s.nh.size() of them). */
static void dis_catchup(dis_t& D, dsearch& S, const dframe& Fcur, const sp_t& s, std::vector<sp_t>& out)
{
    if (s->nh.size() == S.R.size()) { out.push_back(s); return; }
    std::vector<dstate> cur;
    cur.push_back(s->clone());
    for (size_t k = s->nh.size(); k < S.R.size() && !cur.empty(); k++) {
        const dframe& Fk = S.frames[k];
        const dframe& Fk1 = (k + 1 < S.R.size()) ? S.frames[k + 1] : Fcur;
        std::vector<dstate> nxt;
        for (const dstate& st : cur) dis_refine(D, Fk, Fk1, st, S.R[k], nxt);
        cur.swap(nxt);
    }
    for (dstate& x : cur)
        if (dis_fits(D, S.Xrem, x.cand)) out.push_back(std::make_shared<dstate>(std::move(x)));
}

/* Completion LP. For a complete family F the matrices n^j = Psi diag(S_ji/S_0i) Psi^T
   satisfy, for every member k and label j,
     (PF)   sum_{l in F} n^j_kl g_l    = d_j g_k             (g = vacuum coefficient,
                                                               d_j = S_j0/S_00: the vacuum
                                                               column of Psi is an eigenvector)
     (NIM)  sum_{l in F} (n^j_kl)^2    = sum_r N_jj^r x_k^r  (diagonal of (n^j)^2)
   and the trace identity sum_{l in F} x_l = S w. With c_t = number of members still to be
   taken from class t (0 <= c_t <= bound_t), the parts not yet fixed are linear in c, with
   coefficients from each class's overlaps with member k. Infeasible LP: no completion.
   Needs every class refined against every member. */
static int dis_completion_lp(dis_t& D, dsearch& S, const std::vector<sp_t>& P)
{
    const slong n = D.n, R = (slong) S.R.size(), np = (slong) P.size();
    if (np == 0) {
        for (slong j = 0; j < n; j++) if (S.Xrem[j]) return 0;
        return 1;
    }
    auto gof = [&](slong cand) { return D.g[cand]; };
    auto nkl = [&](slong k, slong l) -> const std::vector<slong>& {
        if (k == l) return D.e->sols[S.R[k].cand];
        return k > l ? S.R[k].nh[l] : S.R[l].nh[k];
    };
    struct lrow { double lo, hi; std::vector<std::pair<int, double>> c; };
    std::vector<lrow> rows;
    const double tol = 1e-6;
    int bad = 0;
    auto push = [&](double lo, double hi, std::vector<std::pair<int, double>>& c) {
        if (c.empty()) { if (lo > 0 || hi < 0) bad = 1; return; }      /* constant row */
        double mx = 0.0;                                                 /* range over the box */
        double mn = 0.0;
        for (auto& pr : c) {
            double b = pr.second * P[pr.first]->bound;
            if (b > 0) mx += b; else mn += b;
        }
        if (mx < lo || mn > hi) { bad = 1; return; }
        if (mn >= lo && mx <= hi) return;                                /* never binding */
        rows.push_back({lo, hi, std::move(c)});
    };
    std::vector<std::pair<int, double>> c;
    /* trace */
    for (slong j = 0; j < n && !bad; j++) {
        c.clear();
        for (slong t = 0; t < np; t++) {
            slong a = D.e->sols[P[t]->cand][j];
            if (a) c.push_back({(int) t, (double) a});
        }
        push((double) S.Xrem[j], (double) S.Xrem[j], c);
    }
    for (slong k = 0; k < R && !bad; k++) {
        const double gk = gof(S.R[k].cand);
        for (slong j = 0; j < n && !bad; j++) {
            /* PF */
            double rhs = D.dj[j] * gk;
            for (slong l = 0; l < R; l++) rhs -= nkl(k, l)[j] * gof(S.R[l].cand);
            c.clear();
            for (slong t = 0; t < np; t++) {
                slong a = P[t]->nh[k][j];
                if (a) c.push_back({(int) t, a * gof(P[t]->cand)});
            }
            push(rhs - tol * (1 + std::fabs(rhs)), rhs + tol * (1 + std::fabs(rhs)), c);
            /* NIM */
            if (D.use_caps && !bad) {
                slong left = D.cap[S.R[k].cand * n + j] - S.mused[k][j];
                c.clear();
                for (slong t = 0; t < np; t++) {
                    slong a = P[t]->nh[k][j];
                    if (a) c.push_back({(int) t, (double) (a * a)});
                }
                push((double) left, (double) left, c);
            }
        }
    }
    if (bad) return 0;
    if (rows.empty()) return 1;
    /* drop exact duplicate rows */
    std::sort(rows.begin(), rows.end(), [](const lrow& a, const lrow& b) {
        if (a.c != b.c) return a.c < b.c;
        return a.lo < b.lo;
    });
    std::vector<lrow> uniq;
    for (auto& r : rows) {
        if (!uniq.empty() && uniq.back().c == r.c) {
            uniq.back().lo = std::max(uniq.back().lo, r.lo);
            uniq.back().hi = std::min(uniq.back().hi, r.hi);
            if (uniq.back().lo > uniq.back().hi) return 0;
            continue;
        }
        uniq.push_back(std::move(r));
    }
    D.lp_rows += (slong) uniq.size(); D.lp_calls++;

    glp_prob* lp = glp_create_prob();
    glp_add_cols(lp, (int) np);
    for (slong t = 0; t < np; t++) glp_set_col_bnds(lp, (int) t + 1, GLP_DB, 0.0, (double) P[t]->bound);
    glp_add_rows(lp, (int) uniq.size());
    std::vector<int> ia(1, 0), ja(1, 0);
    std::vector<double> ar(1, 0.0);
    for (size_t r = 0; r < uniq.size(); r++) {
        const lrow& L = uniq[r];
        if (L.hi - L.lo < 1e-300) glp_set_row_bnds(lp, (int) r + 1, GLP_FX, L.lo, L.lo);
        else glp_set_row_bnds(lp, (int) r + 1, GLP_DB, L.lo, L.hi);
        for (auto& pr : L.c) { ia.push_back((int) r + 1); ja.push_back(pr.first + 1); ar.push_back(pr.second); }
    }
    glp_load_matrix(lp, (int) ia.size() - 1, ia.data(), ja.data(), ar.data());
    glp_smcp parm; glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;
    parm.presolve = GLP_ON;
    int ret = glp_simplex(lp, &parm);
    int feasible = 1;
    if (ret == GLP_ENOPFS || (ret == 0 && glp_get_status(lp) == GLP_NOFEAS)) feasible = 0;
    glp_delete_prob(lp);
    return feasible;
}

/* Exact-cover search over classes. Pool states may lag behind the members (their
   bounds stay valid upper bounds, since refinement only splits or removes); a state
   is brought up to date only when the branching column selects it. */
static void dis_rec(dis_t& D, dsearch& S, const dframe& F, const std::vector<sp_t>& P0)
{
    const slong n = D.n;
    if (D.truncated) return;
    if (++D.nodes > SEARCH_NODE_LIMIT) { D.truncated = 1; return; }
    if ((D.nodes & 0xFF) == 0) dis_progress(D, 0);
    if ((slong) S.R.size() == S.T) { dis_record(D, S); return; }

    /* coverage: what the classes can still supply, per component of x
       (component 0 is the size bound, since x^0 = 1) */
    std::vector<slong> cov(n, 0), opts(n, 0);
    for (const sp_t& s : P0) {
        const std::vector<slong>& x = D.e->sols[s->cand];
        for (slong j = 0; j < n; j++) if (x[j]) { cov[j] += s->bound * x[j]; opts[j]++; }
    }
    for (slong j = 0; j < n; j++) if (cov[j] < S.Xrem[j]) return;
    if (!dis_knapsack(D, P0, std::vector<char>(P0.size(), 0), S.Xrem)) return;

    /* bring every class up to date and test the completion LP */
    std::vector<sp_t> Pc;
    if (!S.R.empty()) {
        for (const sp_t& s : P0) dis_catchup(D, S, F, s, Pc);
        if (!dis_completion_lp(D, S, Pc)) { D.lp_prunes++; return; }
    }
    const std::vector<sp_t>& Pu = S.R.empty() ? P0 : Pc;

    /* Branch on the component still owed with the fewest classes able to supply it:
       every completion contains one of them (exact-cover column choice). With no
       members yet, branch on component 0 (every class) instead: first members with
       small x make every later comparison cheap, and each branch drops the earlier
       ones. Within a branch list, classes with small x come first. */
    slong jb = 0;
    if (!S.R.empty())
        for (slong j = 1; j < n; j++)
            if (S.Xrem[j] > 0 && (jb == 0 || opts[j] < opts[jb])) jb = j;

    std::vector<sp_t> P;                 /* the column's states brought up to date */
    for (const sp_t& s : Pu) {
        if (D.e->sols[s->cand][jb]) dis_catchup(D, S, F, s, P);
        else P.push_back(s);
    }
    std::vector<size_t> list;
    for (size_t k = 0; k < P.size(); k++) if (D.e->sols[P[k]->cand][jb]) list.push_back(k);
    auto norm2 = [&](size_t k) {
        slong s = 0;
        for (slong j = 0; j < n; j++) s += D.e->sols[P[k]->cand][j] * D.e->sols[P[k]->cand][j];
        return s;
    };
    std::stable_sort(list.begin(), list.end(), [&](size_t a, size_t b) { return norm2(a) < norm2(b); });

    std::vector<char> dropped(P.size(), 0);
    if (S.R.empty()) D.root_of = (slong) list.size();
    for (size_t li = 0; li < list.size() && !D.truncated; li++) {
        const size_t k = list[li];
        if (S.R.empty()) D.root_at = (slong) li + 1;
        if (li > 0) {                   /* what is left must still cover the remainder */
            std::vector<slong> c2(n, 0);
            for (size_t t = 0; t < P.size(); t++) {
                if (dropped[t]) continue;
                const std::vector<slong>& x = D.e->sols[P[t]->cand];
                for (slong j = 0; j < n; j++) c2[j] += P[t]->bound * x[j];
            }
            int enough = 1;
            for (slong j = 0; j < n && enough; j++) enough = (c2[j] >= S.Xrem[j]);
            if (!enough) break;
            if (!dis_knapsack(D, P, dropped, S.Xrem, 0)) break;   /* nothing later completes */
        }
        const dstate& v = *P[k];
        const std::vector<slong>& xv = D.e->sols[v.cand];
        std::vector<slong> Xn(S.Xrem);
        for (slong j = 0; j < n; j++) Xn[j] -= xv[j];

        /* the members' NIM-rep budgets with v added */
        std::vector<std::vector<slong>> saved(S.mused);
        int within = 1;
        for (size_t m = 0; m < S.R.size() && within; m++)
            for (slong j = 0; j < n; j++) {
                S.mused[m][j] += v.nh[m][j] * v.nh[m][j];
                if (S.mused[m][j] > D.cap[S.R[m].cand * n + j]) within = 0;
            }
        if (within) {
            std::vector<sp_t> P2;       /* everything not dropped that still fits; v's own
                                           entry stays, and yields its class-mates later */
            for (size_t t = 0; t < P.size(); t++)
                if (!dropped[t] && dis_fits(D, Xn, P[t]->cand)) P2.push_back(P[t]);
            dframe Fn = dis_extend_frame(D, F, v);
            S.mused.push_back(v.used);
            S.R.push_back(v.clone());
            S.c_at.push_back(F.c);
            S.frames.push_back(F.clone());
            std::swap(S.Xrem, Xn);
            dis_rec(D, S, Fn, P2);
            std::swap(S.Xrem, Xn);
            S.frames.pop_back();
            S.c_at.pop_back();
            S.R.pop_back();
        }
        S.mused.swap(saved);
        dropped[k] = 1;                 /* later branches contain no state of this class */
    }
}

/* S w with w_i = m_i / S_0i on the support pattern U; 0 if it is not a vector of
   nonnegative integers (then no complete family has this support). */
static int dis_trace_target(dis_t& D, std::vector<slong>& X)
{
    search_t* e = D.e;
    const slong n = D.n;
    nf_vec w(n, e->field);
    for (slong i = 0; i < n; i++) {
        if (!D.inA[i]) continue;
        nf_elem_set_si(w[i], D.mult[i], e->field);
        nf_elem_div(w[i], w[i], entry(e->S, 0, i, e->ctx), e->field);
    }
    nf_elem_t s, t; nf_elem_init(s, e->field); nf_elem_init(t, e->field);
    fmpq_poly_t p; fmpq_poly_init(p);
    int ok = 1;
    X.assign(n, 0);
    for (slong j = 0; j < n && ok; j++) {
        nf_elem_zero(s, e->field);
        for (slong i = 0; i < n; i++) {
            if (!D.inA[i]) continue;
            nf_elem_mul(t, entry(e->S, j, i, e->ctx), w[i], e->field);
            nf_elem_add(s, s, t, e->field);
        }
        if (!nf_elem_is_integer(s, e->field)) { ok = 0; break; }
        nf_elem_get_fmpq_poly(p, s, e->field);
        fmpz_t z; fmpz_init(z);
        fmpq_poly_get_coeff_fmpz(z, p, 0);
        if (fmpz_sgn(z) < 0 || !fmpz_fits_si(z)) ok = 0;
        else X[j] = fmpz_get_si(z);
        fmpz_clear(z);
    }
    fmpq_poly_clear(p);
    nf_elem_clear(t, e->field); nf_elem_clear(s, e->field);
    return ok;
}

/* Families complete on the support pattern U, from no members. */
static void dis_pattern(dis_t& D, const std::vector<char>& U)
{
    const slong n = D.n, M = D.M;
    D.inA = U;
    slong N = 0;
    for (slong i = 0; i < n; i++) if (U[i]) N += D.mult[i];

    std::vector<slong> Xtot;
    if (!dis_trace_target(D, Xtot) || Xtot[0] != N) return;

    std::vector<sp_t> P;
    for (slong B = 0; B < M; B++) {
        int inside = 1;
        for (slong i = 0; i < n; i++) if (D.supp[B * n + i] && !U[i]) inside = 0;
        if (!inside || !D.pf_ok[B] || !dis_fits(D, Xtot, B)) continue;
        dstate b = dis_raw(D, B);
        b.bound = D.cnt[B];
        P.push_back(std::make_shared<dstate>(std::move(b)));
    }
    slong total = 0;
    for (const sp_t& s : P) total += s->bound;
    if (total < N) return;

    dframe F0 = dis_empty_frame(D);
    dsearch S;
    S.T = N;
    S.Xrem = Xtot;
    dis_rec(D, S, F0, P);
}

static void disentangle(search_t* e)
{
    const slong n = e->n, M = (slong) e->sols.size();
    flint_printf("\n--- disentangling %wd candidates with the mutual Cardy condition ---\n", M);
    if (M == 0) return;

    dis_t D;
    D.e = e; D.n = n; D.M = M;
    D.checked = D.accepted = D.refinements = D.nodes = D.open_subset = 0;
    D.exact = 0; D.rejected = 0;
    D.truncated = 0;
    D.patterns_done = D.npatterns = 0;
    D.root_at = D.root_of = 0;
    D.t0 = D.t_last = dis_now();
    arb_init(D.gball);
    qqbar_get_arb(D.gball, e->gen, NUM_PREC);
    D.mult.assign(n, 0); D.off.assign(n, 0);
    D.L = 0;
    for (slong i = 0; i < n; i++) { D.mult[i] = e->mult[i]; D.off[i] = D.L; D.L += e->mult[i] + 1; }

    /* exact |b_a^i|^2 for every i (it is exactly 0 for i not in I) */
    D.q = (nf_elem_struct*) flint_malloc(M * n * sizeof(nf_elem_struct));
    for (slong k = 0; k < M * n; k++) nf_elem_init(D.q + k, e->field);
    D.supp.assign(M * n, 0);
    D.qd.assign(M * n, 0.0);
    {
        nf_vec P(n, e->field);
        for (slong a = 0; a < M; a++) {
            dis_Sn(P, D, e->sols[a]);
            for (slong i = 0; i < n; i++) {
                nf_elem_set(D.q + a * n + i, P[i], e->field);
                D.supp[a * n + i] = !nf_elem_is_zero(P[i], e->field);
                if (D.supp[a * n + i]) D.qd[a * n + i] = dis_d(P[i], D);
            }
        }
    }
    arb_t g, v; arb_init(g); arb_init(v);
    qqbar_get_arb(g, e->gen, NUM_PREC);
    D.Sd.assign(n * n, 0.0);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) {
            nf_elem_get_arb(v, entry(e->S, i, j, e->ctx), g, e->field, NUM_PREC);
            D.Sd[i * n + j] = arf_get_d(arb_midref(v), ARF_RND_NEAR);
        }
    arb_clear(v); arb_clear(g);
    D.Qs.assign((n + 1) * n, 0.0);
    for (slong j = n - 1; j >= 0; j--)
        for (slong i = 0; i < n; i++)
            D.Qs[j * n + i] = std::sqrt(D.Qs[(j + 1) * n + i] * D.Qs[(j + 1) * n + i]
                                        + D.Sd[i * n + j] * D.Sd[i * n + j]);
    D.dj.assign(n, 0.0);
    for (slong j = 0; j < n; j++) D.dj[j] = D.Sd[j * n] / D.Sd[0];
    D.g.assign(M, 0.0);
    D.pf_ok.assign(M, 1);
    slong npf = 0;
    for (slong a = 0; a < M; a++) {
        D.g[a] = std::sqrt(std::max(0.0, D.qd[a * n]));
        for (slong j = 0; j < n; j++)
            if (e->sols[a][j] > D.dj[j] + 1e-7) { D.pf_ok[a] = 0; break; }
        npf += D.pf_ok[a];
    }
    printf("  candidates with x^j <= d_j for all j (others are in no complete family): %ld of %ld\n",
           (long) npf, (long) M);
    D.lp_prunes = D.lp_rows = D.lp_calls = 0;

    /* NIM-rep budgets: N_jj^r by Verlinde (numerically, they are integers) */
    D.use_caps = 1;
    D.cap.assign(M * n, 0);
    {
        std::vector<slong> Njj(n * n, 0);
        for (slong j = 0; j < n && D.use_caps; j++)
            for (slong r = 0; r < n; r++) {
                double s = 0.0;
                for (slong i = 0; i < n; i++)
                    s += D.Sd[j * n + i] * D.Sd[j * n + i] * D.Sd[r * n + i] / D.Sd[i];
                double k = std::nearbyint(s);
                if (std::fabs(s - k) > 1e-6 || k < 0) { D.use_caps = 0; break; }
                Njj[j * n + r] = (slong) k;
            }
        if (!D.use_caps) printf("  fusion rules from S are not integral; NIM-rep budgets not used\n");
        else
            for (slong a = 0; a < M; a++)
                for (slong j = 0; j < n; j++) {
                    slong y = 0;
                    for (slong r = 0; r < n; r++) y += Njj[j * n + r] * e->sols[a][r];
                    D.cap[a * n + j] = y - e->sols[a][j] * e->sols[a][j];
                }
    }

    slong mmax = 0;
    for (slong i = 0; i < n; i++) mmax = FLINT_MAX(mmax, D.mult[i]);
    printf("  [%6.1fs] |b|^2 computed exactly for %ld candidates, %ld Ishibashi states "
           "over %ld labels (at most %ld per label), searching families\n",
           dis_now() - D.t0, (long) M, (long) e->dimI, (long) e->nI, (long) mmax);
    fflush(stdout);
    /* trace bound per candidate: members of candidate B in a complete family */
    D.cnt.assign(M, 0);
    for (slong a = 0; a < M; a++) {
        double best = 1e300;
        for (slong i = 0; i < n; i++) {
            if (!D.supp[a * n + i]) continue;
            double w = D.Sd[i] * dis_d(D.q + a * n + i, D);
            best = std::min(best, D.mult[i] / w);
        }
        D.cnt[a] = FLINT_MAX(1, (slong) std::floor(best + 1e-6));
    }

    /* support patterns: all of I first, then each candidate's support */
    std::vector<std::vector<char>> pats;
    pats.push_back(std::vector<char>(e->inI, e->inI + n));
    for (slong a = 0; a < M && SUBSET_FAMILIES; a++) {
        std::vector<char> u(D.supp.begin() + a * n, D.supp.begin() + (a + 1) * n);
        if (std::find(pats.begin(), pats.end(), u) == pats.end()) pats.push_back(u);
    }
    D.npatterns = (slong) pats.size();
    for (size_t k = 0; k < pats.size() && !D.truncated; k++) {
        dis_pattern(D, pats[k]);
        D.patterns_done++;
        dis_progress(D, 0);
    }
    dis_progress(D, 1);

    /* report: complete for I first, then closed subset families, larger first */
    std::vector<slong> ncov(D.found.size(), 0), in_complete(M, 0);
    for (size_t f = 0; f < D.found.size(); f++) {
        for (slong i = 0; i < n; i++)
            if (D.found[f].U[i]) ncov[f] += D.mult[i];
        if (D.found[f].kind == 2) for (const member_t& m : D.found[f].mem) in_complete[m.cand]++;
    }
    std::vector<size_t> ord(D.found.size());
    for (size_t f = 0; f < ord.size(); f++) ord[f] = f;
    std::stable_sort(ord.begin(), ord.end(), [&](size_t u, size_t w) {
        if (D.found[u].kind != D.found[w].kind) return D.found[u].kind > D.found[w].kind;
        return D.found[u].mem.size() > D.found[w].mem.size();
    });
    slong nk[3] = {0, 0, 0};
    for (const family_t& f : D.found) nk[f.kind]++;

    flint_printf("refinements: %wd; overlap checks: %wd, %wd accepted; search nodes: %wd; "
                 "families failing the exact check: %wd%s\n",
                 D.refinements, D.checked, D.accepted, D.nodes, D.rejected,
                 D.truncated ? "  (SEARCH NODE LIMIT HIT: list below may be incomplete)" : "");
    flint_printf("complete families (up to gauge): %wd for I, %wd closed on a subset of I"
                 " (%wd more on a subset but extendable, not listed)\n",
                 nk[2], nk[1], D.open_subset);

    for (size_t c = 0; c < ord.size() && c < (size_t) MAX_FAMILIES_SHOWN; c++) {
        const family_t& F = D.found[ord[c]];
        if (F.kind == 2)
            flint_printf("\nfamily %wd: %wd boundaries, complete for I\n", (slong) c + 1, (slong) F.mem.size());
        else
            flint_printf("\nfamily %wd: %wd boundaries, complete on %wd of the %wd Ishibashi states\n",
                         (slong) c + 1, (slong) F.mem.size(), ncov[ord[c]], e->dimI);
        for (const member_t& m : F.mem) {
            flint_printf("  cand %wd  x = ", m.cand + 1);
            for (slong j = 0; j < n; j++) flint_printf("%wd ", e->sols[m.cand][j]);
            flint_printf("\n");
        }
    }
    if (ord.size() > (size_t) MAX_FAMILIES_SHOWN)
        flint_printf("  ... %wd more families not shown\n", (slong) ord.size() - MAX_FAMILIES_SHOWN);

    flint_printf("\ncandidates in no family complete for I:");
    int any = 0;
    for (slong a = 0; a < M; a++) if (!in_complete[a]) { flint_printf(" %wd", a + 1); any = 1; }
    flint_printf(any ? "\n" : " none\n");

    /* Full list of boundaries in the notation of the search output: b^i for labels
       with one Ishibashi state, (b^{i,1}, ..., b^{i,m}) for labels with m, in the
       Gram-Schmidt basis of the family's members in the order listed. */
    slong mx = 0;
    for (slong i = 0; i < n; i++) mx = FLINT_MAX(mx, D.mult[i]);
    arb_t gb, c1, t1, bnum, bden;
    arb_init(gb); arb_init(c1); arb_init(t1); arb_init(bnum); arb_init(bden);
    qqbar_get_arb(gb, e->gen, NUM_PREC);
    const slong inum = TRANS_NUM - INDEX_BASE, iden = TRANS_DEN - INDEX_BASE;
    const int have_T = inum >= 0 && inum < n && iden >= 0 && iden < n
                       && e->inI[inum] && e->inI[iden];
    const int have_RT = have_T && D.mult[inum] == 1 && D.mult[iden] == 1;

    /* coordinate j of member m at label i, in the printed basis: y_j / sqrt(d_j) along
       earlier pivots, sqrt(r) on its own new axis if it raised the rank, 0 beyond */
    auto coord = [&](arb_t out, const family_t& F, const member_t& m, slong i, slong j) {
        const slong o = D.off[i], mi = D.mult[i], ck = m.c_at[i];
        if (j < ck) {
            nf_elem_get_arb(out, m.ex[o + j], gb, e->field, NUM_PREC);
            nf_elem_get_arb(t1, F.frame.d[o + j], gb, e->field, NUM_PREC);
            arb_sqrt(t1, t1, NUM_PREC);
            arb_div(out, out, t1, NUM_PREC);
        } else if (j == ck && F.U[i] && ck < mi && !nf_elem_is_zero(m.ex[o + mi], e->field)) {
            nf_elem_get_arb(out, m.ex[o + mi], gb, e->field, NUM_PREC);
            arb_sqrt(out, out, NUM_PREC);
        } else arb_zero(out);
    };

    flint_printf("\n=== boundaries ===\n");
    if (mx > 1)
        flint_printf("labels with m > 1 Ishibashi states are printed as (b^{i,1}, ..., b^{i,m}), in the\n"
                     "Gram-Schmidt basis of the family's boundaries in the order listed\n");
    if (have_RT)
        flint_printf("R = (1 + r)/2, T = (1 - r)/2 with r = b^%d/b^%d. Flipping the sign of the\n"
                     "Ishibashi state of label %d for every boundary leaves all Cardy overlaps\n"
                     "unchanged and exchanges R and T, so the data cannot fix that sign. Convention\n"
                     "used: in each family, the boundary with the smallest open spectrum (smallest\n"
                     "|x|^2, ties to the first listed) has r >= 0; it is marked [reference]\n",
                     TRANS_NUM, TRANS_DEN, TRANS_NUM);
    if (nk[2] + nk[1] == 0) flint_printf("no complete family found\n");
    if (nk[2] == 0 && e->dimI == e->nI)
        flint_printf("note: every label was given one Ishibashi state. For a non-diagonal modular\n"
                     "invariant, list each label as many times as it has Ishibashi states in the\n"
                     "Ishibashi index input (its diagonal mass-matrix entry Z_ii)\n");
    slong shown = 0;
    for (size_t c = 0; c < ord.size(); c++) {
        if (shown++ >= MAX_FAMILIES_SHOWN) {
            flint_printf("... further complete families not shown\n");
            break;
        }
        const family_t& F = D.found[ord[c]];
        if (F.kind == 2)
            flint_printf("\nfamily %wd (complete for I): %wd boundaries\n", (slong) c + 1, (slong) F.mem.size());
        else
            flint_printf("\nfamily %wd (complete on %wd of the %wd Ishibashi states): %wd boundaries\n",
                         (slong) c + 1, ncov[ord[c]], e->dimI, (slong) F.mem.size());
        /* reference for the sign of r */
        slong ref = -1, refn = WORD_MAX;
        int sref = 1;
        if (have_RT) {
            for (size_t k = 0; k < F.mem.size(); k++) {
                if (!D.supp[F.mem[k].cand * n + iden]) continue;
                slong nn = 0;
                for (slong j = 0; j < n; j++) nn += e->sols[F.mem[k].cand][j] * e->sols[F.mem[k].cand][j];
                if (nn < refn) { refn = nn; ref = (slong) k; }
            }
            if (ref >= 0) {
                coord(bnum, F, F.mem[ref], inum, 0);
                coord(bden, F, F.mem[ref], iden, 0);
                arb_div(bnum, bnum, bden, NUM_PREC);
                if (arb_is_negative(bnum)) sref = -1;
            }
        }
        for (size_t k = 0; k < F.mem.size(); k++) {
            const member_t& m = F.mem[k];
            flint_printf("boundary %wd.%wd (candidate %wd)%s\nx = ", (slong) c + 1, (slong) k + 1, m.cand + 1,
                         (slong) k == ref ? " [reference]" : "");
            for (slong j = 0; j < n; j++) flint_printf("%wd ", e->sols[m.cand][j]);
            flint_printf("\n");
            for (slong i = 0; i < n; i++) {
                if (!e->inI[i]) continue;
                const slong mi = D.mult[i];
                flint_printf("  b^%wd = ", i + INDEX_BASE);
                if (mi > 1) flint_printf("(");
                for (slong j = 0; j < mi; j++) {
                    coord(c1, F, m, i, j);
                    if (j) flint_printf(", ");
                    arb_printn(c1, 15, 0);
                }
                if (mi > 1) flint_printf(")");
                flint_printf("\n");
            }
            if (have_RT) {
                if (!D.supp[m.cand * n + iden]) flint_printf("  R, T undefined (b^%d = 0)\n", TRANS_DEN);
                else {
                    coord(bnum, F, m, inum, 0);
                    coord(bden, F, m, iden, 0);
                    arb_div(bnum, bnum, bden, NUM_PREC);           /* r */
                    if (sref < 0) arb_neg(bnum, bnum);
                    arb_add_ui(bden, bnum, 1, NUM_PREC);           /* R = (1 + r)/2 */
                    arb_mul_2exp_si(bden, bden, -1);
                    arb_sub_ui(bnum, bnum, 1, NUM_PREC);           /* T = (1 - r)/2 */
                    arb_neg(bnum, bnum);
                    arb_mul_2exp_si(bnum, bnum, -1);
                    flint_printf("  R = ");
                    arb_printn(bden, 15, 0);
                    flint_printf(",  T = ");
                    arb_printn(bnum, 15, 0);
                    flint_printf("\n");
                }
            } else if (have_T) {
                flint_printf("  R, T not computed: label %d or %d has several Ishibashi states\n",
                             TRANS_NUM, TRANS_DEN);
            }
        }
    }
    arb_clear(bden); arb_clear(bnum); arb_clear(t1); arb_clear(c1); arb_clear(gb);
    fflush(stdout);

    for (slong k = 0; k < M * n; k++) nf_elem_clear(D.q + k, e->field);
    flint_free(D.q);
    arb_clear(D.gball);
}

/*
 * Find all nonnegative integer x with
 *     sum_j S*_ij x^j  = 0   for i not in I
 *     sum_j S*_ij x^j  > 0   for i in I
 *     x^i <= max_{j in I} S_ij / S_0j
 * and print |b^i| = sqrt(sum_j S*_ij x^j) for i in I.
 *
 * All values in the field are real, so S* = S and S^2 = 1. Then
 * x^0 = sum_{j in I} S_0j (Sx)_j > 0, while the bound gives x^0 <= 1, so x^0 = 1
 * (it is enforced through the box, since lb_0 = ub_0 = 1).
 *
 * The equalities are row-reduced; the search branches only on the free
 * columns and recovers the pivot columns exactly at the leaves, where
 * positivity is checked exactly. Each free column's range comes from a
 * floating-point LP, but only its dual multipliers are used: the range itself
 * is a certified bound over the exact constraints (certified_lower), and a
 * subtree is pruned only when its emptiness is proven. The search is therefore
 * complete as well as sound.
 * Returns the number of solutions, or -1 on error.
 */
int solve_boundary_bootstrap(const slong* indices, slong nindices,
                             gr_mat_struct* matrix, nf_struct* field,
                             qqbar_struct* field_generator, gr_ctx_struct* context)
{
    const slong n      = gr_mat_nrows(matrix, context);
    const slong degree = fmpq_poly_degree(field->pol);

    if (gr_mat_ncols(matrix, context) != n || nindices <= 0) return -1;

    /* I as a mask; a label listed k times carries k Ishibashi states */
    std::vector<char> inI(n, 0);
    std::vector<int>  mult(n, 0);
    for (slong k = 0; k < nindices; k++) {
        slong idx = indices[k] - INDEX_BASE;
        if (idx < 0 || idx >= n) {
            flint_fprintf(stderr, "ERROR: Ishibashi index %wd out of range\n", indices[k]);
            return -1;
        }
        inI[idx] = 1;
        mult[idx]++;
    }
    slong nout = 0;
    for (slong i = 0; i < n; i++) nout += !inI[i];

    /* 1. Equalities for i not in I, one rational row per coordinate,
          row-reduced over Q and cleared to integers row by row. */
    slong rank;
    fmpz_mat_t Z;
    {
        fmpq_mat_t  A;  fmpq_mat_init(A, FLINT_MAX(nout, 1) * degree, n);
        fmpq_mat_t  R;  fmpq_mat_init(R, FLINT_MAX(nout, 1) * degree, n);
        fmpq_poly_t p;  fmpq_poly_init(p);

        for (slong i = 0, r = 0; i < n; i++) {
            if (inI[i]) continue;
            for (slong j = 0; j < n; j++) {
                nf_elem_get_fmpq_poly(p, entry(matrix, i, j, context), field);
                for (slong l = 0; l < degree; l++)
                    fmpq_poly_get_coeff_fmpq(fmpq_mat_entry(A, r * degree + l, j), p, l);
            }
            r++;
        }

        rank = fmpq_mat_rref(R, A);
        fmpz_mat_init(Z, FLINT_MAX(rank, 1), n);

        fmpq_mat_t W; fmpq_mat_window_init(W, R, 0, 0, FLINT_MAX(rank, 1), n);
        fmpq_mat_get_fmpz_mat_rowwise(Z, NULL, W);
        fmpq_mat_window_clear(W);

        fmpq_poly_clear(p);
        fmpq_mat_clear(R);
        fmpq_mat_clear(A);
    }

    int     status = -1;
    arb_ptr Sa = _arb_vec_init(n * n);
    arb_t   g;  arb_init(g);

    search_t e;
    e.lp = NULL;
    e.el = NULL;
    e.mult = mult;
    e.dimI = nindices;
    e.n  = n;
    e.Z  = Z;
    e.rank = rank;
    e.fallbacks = e.lp_prunes = e.elastic_calls = e.elastic_certified = 0;
    e.max_gap = 0.0;
    fmpz_init(e.acc);
    arb_init(e.v);
    e.nodes = e.leaves = e.nsolutions = 0;

    do {
        /* 2. Bounds */
        e.lb.resize(n); e.ub.resize(n); e.x.assign(n, 0);
        int bstatus = build_bounds(e.lb.data(), e.ub.data(), matrix, inI.data(),
                                   field, field_generator, context);
        if (bstatus < 0) break;
        if (bstatus > 0) { flint_printf("empty box, no solutions\n"); status = 0; break; }
        e.cl = e.lb;
        e.cu = e.ub;

        /* 3. Pivot and free columns of the reduced equalities */
        std::vector<char> is_pivot(n, 0);
        for (slong r = 0; r < rank; r++) {
            slong j = 0;
            while (j < n && fmpz_is_zero(fmpz_mat_entry(Z, r, j))) j++;
            e.piv.push_back(j);
            is_pivot[j] = 1;
        }
        for (slong j = 0; j < n; j++)
            if (!is_pivot[j]) e.order.push_back(j);

        /* branch on the tightest boxes first */
        std::stable_sort(e.order.begin(), e.order.end(), [&](slong a, slong b) {
            return e.ub[a] - e.lb[a] < e.ub[b] - e.lb[b];
        });

        /* 4. Numeric S, used by the LP and by the exact check */
        qqbar_get_arb(g, field_generator, NUM_PREC);
        for (slong i = 0; i < n; i++)
            if (inI[i])
                for (slong j = 0; j < n; j++)
                    nf_elem_get_arb(Sa + i * n + j, entry(matrix, i, j, context),
                                    g, field, NUM_PREC);

        /* 5. LPs. Main:   Z x = 0,  S_I x >= -LP_RELAX,  box.
                 Elastic: Z x + s+ - s- = 0,  S_I x + s >= -LP_RELAX,  box,
                          slacks >= 0,  min sum of slacks.
              Rows 1..rank are the equalities, rank+1..rank+nI the positivity
              rows, in both. The float coefficients only steer the multipliers. */
        slong nI = n - nout;
        e.nI = nI;
        for (slong i = 0; i < n; i++) if (inI[i]) e.Irow.push_back(i);
        e.y.assign(FLINT_MAX(rank, 1), 0.0);
        e.z.assign(nI, 0.0);

        std::vector<int>    ia(1, 0), ja(1, 0);
        std::vector<double> ar(1, 0.0);
        for (slong r = 0; r < rank; r++)
            for (slong j = 0; j < n; j++) {
                if (fmpz_is_zero(fmpz_mat_entry(Z, r, j))) continue;
                ia.push_back((int) r + 1); ja.push_back((int) j + 1);
                ar.push_back(fmpz_get_d(fmpz_mat_entry(Z, r, j)));
            }
        for (slong k = 0; k < nI; k++)
            for (slong j = 0; j < n; j++) {
                ia.push_back((int) (rank + k) + 1); ja.push_back((int) j + 1);
                ar.push_back(arf_get_d(arb_midref(Sa + e.Irow[k] * n + j), ARF_RND_NEAR));
            }
        const int nnz = (int) ia.size() - 1;

        glp_term_out(GLP_OFF);
        e.lp = glp_create_prob();
        e.el = glp_create_prob();
        for (glp_prob* lp : {e.lp, e.el}) {
            glp_add_cols(lp, (int) n);
            for (slong j = 0; j < n; j++) set_col_box(lp, (int) j + 1, e.lb[j], e.ub[j]);
            glp_add_rows(lp, (int) (rank + nI));
            for (slong r = 0; r < rank; r++)
                glp_set_row_bnds(lp, (int) r + 1, GLP_FX, 0.0, 0.0);
            for (slong k = 0; k < nI; k++)
                glp_set_row_bnds(lp, (int) (rank + k) + 1, GLP_LO, -LP_RELAX, 0.0);
        }

        glp_load_matrix(e.lp, nnz, ia.data(), ja.data(), ar.data());

        /* slack columns of the elastic LP: s+_r, s-_r (r < rank), then s_k */
        {
            const int first = (int) n + 1;
            const int nslack = (int) (2 * rank + nI);
            glp_add_cols(e.el, nslack);
            for (int c = first; c < first + nslack; c++) {
                glp_set_col_bnds(e.el, c, GLP_LO, 0.0, 0.0);
                glp_set_obj_coef(e.el, c, 1.0);
            }
            for (slong r = 0; r < rank; r++) {
                ia.push_back((int) r + 1); ja.push_back(first + (int) (2 * r));     ar.push_back( 1.0);
                ia.push_back((int) r + 1); ja.push_back(first + (int) (2 * r) + 1); ar.push_back(-1.0);
            }
            for (slong k = 0; k < nI; k++) {
                ia.push_back((int) (rank + k) + 1); ja.push_back(first + (int) (2 * rank + k));
                ar.push_back(1.0);
            }
            glp_load_matrix(e.el, (int) ia.size() - 1, ia.data(), ja.data(), ar.data());
            glp_set_obj_dir(e.el, GLP_MIN);
        }

        glp_scale_prob(e.lp, GLP_SF_AUTO);
        glp_scale_prob(e.el, GLP_SF_AUTO);

        glp_init_smcp(&e.parm);
        e.parm.msg_lev  = GLP_MSG_OFF;
        e.parm.presolve = GLP_OFF;

        /* 6. Search */
        e.inI = inI.data(); e.Sa = Sa; e.S = matrix;
        e.field = field; e.gen = field_generator; e.ctx = context;

        flint_printf("equations: %wd, free variables: %wd, positivity rows: %wd\n",
                     rank, (slong) e.order.size(), nI);
        fflush(stdout);

        dfs(&e, 0);

        flint_printf("nodes: %wd  leaves: %wd  solutions: %wd\n",
                     e.nodes, e.leaves, e.nsolutions);
        flint_printf("certified pruning: LP prunes %wd, elastic %wd/%wd proven, "
                     "box fallbacks %wd, max LP-vs-certified gap %.3g\n",
                     e.lp_prunes, e.elastic_certified, e.elastic_calls,
                     e.fallbacks, e.max_gap);
        status = (int) e.nsolutions;

        if (DISENTANGLE) disentangle(&e);
    } while (0);

    if (e.lp) glp_delete_prob(e.lp);
    if (e.el) glp_delete_prob(e.el);
    arb_clear(e.v);
    fmpz_clear(e.acc);
    arb_clear(g);
    _arb_vec_clear(Sa, n * n);
    fmpz_mat_clear(Z);
    return status;
}
