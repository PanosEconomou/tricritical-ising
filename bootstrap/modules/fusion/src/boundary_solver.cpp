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

/* Transmission coefficient T = 1 - |b^TRANS_NUM| / |b^TRANS_DEN| (labels as in
   the input files, i.e. with INDEX_BASE) */
#define TRANS_NUM  17
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

    /* transmission coefficient */
    slong inum = TRANS_NUM - INDEX_BASE, iden = TRANS_DEN - INDEX_BASE;
    if (inum >= 0 && inum < n && iden >= 0 && iden < n && e->inI[inum] && e->inI[iden]) {
        arb_t bnum, bden;
        arb_init(bnum); arb_init(bden);

        row_sign(bnum, e->S, inum, e->x.data(), n, e->Sa, e->field, e->gen, e->ctx);
        row_sign(bden, e->S, iden, e->x.data(), n, e->Sa, e->field, e->gen, e->ctx);
        arb_sqrt(bnum, bnum, NUM_PREC);
        arb_sqrt(bden, bden, NUM_PREC);

        arb_div(bnum, bnum, bden, NUM_PREC);   /* |b^num| / |b^den| */
        arb_neg(bnum, bnum);
        arb_add_ui(bnum, bnum, 1, NUM_PREC);   /* 1 - ratio */

        flint_printf("  T = 1 - |b^%d|/|b^%d| = ", TRANS_NUM, TRANS_DEN);
        arb_printn(bnum, 15, 0);
        flint_printf("\n");

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
 * Virasoro (chiral-algebra) label i: the number of times i is listed in the
 * Ishibashi index input. For m_i = 1, p^i = +-|b_a^i||b_b^i|; for m_i = 2 it
 * can be anything with |p^i| <= |b_a^i||b_b^i|. (Real vectors suffice: n real
 * and S real force every p real, and a real Gram matrix of complex rank m is
 * realised in R^m.) m_i <= 2 is supported.
 *
 * Facts used, all from S real and S^2 = 1:
 *   - |n_ab|^2 = |p_ab|^2 <= x_a . x_b, with equality when all m_i = 1;
 *   - distinct elementary boundaries have n_ab^0 = 0 (Cauchy-Schwarz with
 *     S_0i > 0 and x^0 = 1), so the vectors (sqrt(S_0i) b_a^i)_i of a family are
 *     orthonormal: a family has at most sum_i m_i members over the labels it
 *     uses, and a complete one is a basis.
 *
 * Family search. Each family is searched from an anchor A whose support
 * contains every member's. The gauge (O(m_i) acting on the Ishibashi states of
 * each label) is fixed by putting b_A^i along the first axis, b_A^i = (|b_A^i|, 0).
 * Then for every other state B, n_AB fixes its first component p_AB^i/|b_A^i|
 * and, for m_i = 2, the second up to a sign:  sqrt(q_A q_B - p_AB^2) / |b_A^i|.
 * With those signs chosen, every state is an explicit vector, the family is a
 * clique of size sum_{i in supp A} m_i in a fixed compatibility graph, and two
 * states B, C are compatible iff for every label
 *     (q_A P - p_AB p_AC)^2 = (q_A q_B - p_AB^2)(q_A q_C - p_AC^2),
 *     sign(q_A P - p_AB p_AC) = s_B s_C  (when the right side is nonzero),
 * with P = (S n_BC)^i: exact identities in the number field, no square roots.
 * Families that differ by the leftover reflections are merged (they have the
 * same n_kl). Candidate n are found numerically and always verified exactly.
 * Families in which no member's support contains all the others' are not
 * searched (none occur in the Cardy case: the vacuum boundary has full support).
 */

#define DISENTANGLE        1
#define MAX_FAMILIES_SHOWN 10
#define CLIQUE_NODE_LIMIT  50000000  /* per run; a warning is printed if hit */
#define DIS_TOL            1e-7      /* numeric pre-filter only; acceptance is exact */
#define PROGRESS_SECONDS   2.0       /* minimum time between progress lines */

typedef std::vector<signed char> signs_t;

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
    nf_elem_struct* operator[](slong i) const { return v + i; }
};

struct pair_opt_t {
    std::vector<slong> n;      /* open-string multiplicities n_ab^j, with n^0 = 0 */
    signs_t            sg;     /* sign of p_ab^i = (S n)_i, exact */
};

/* a boundary state in the gauge of an anchor A */
struct vtx_t {
    slong               cand;
    std::vector<slong>  nA;    /* n_{A, this}; x_A for the anchor itself */
    signs_t             s2;    /* sign of the second component (m_i = 2, residual > 0), else 0 */
    std::vector<double> w1, w2;
    nf_vec              p;     /* p_{A, this}^i, exact */
    nf_vec              res;   /* q_A^i q_this^i - (p^i)^2, exact */
};

struct member_t {
    slong               cand, anchor;
    std::vector<slong>  nA;
    signs_t             s2;
    std::vector<double> w1, w2;
};

typedef std::vector<std::vector<std::vector<slong>>> nmat_t;   /* n[k][l] = n_{kl} */

struct dis_t {
    search_t*  e;
    slong      n, M;
    nf_elem_struct* q;                 /* q[a*n + i] = (S x_a)_i = |b_a^i|^2 exactly */
    std::vector<double> Sd;            /* all of S, numerically */
    std::vector<double> Qs;            /* Qs[j*n + i] = sqrt(sum_{k >= j} S_ik^2) */
    std::vector<double> r;             /* r[a*n + i] = |b_a^i| */
    std::vector<char>   supp;          /* supp[a*n + i] = (b_a^i != 0), exact */
    std::vector<std::vector<pair_opt_t>> opts;
    std::vector<char>   have;
    slong      checked, accepted;

    std::vector<std::vector<member_t>> found;      /* complete families, up to gauge */
    std::vector<nmat_t>    found_n;
    std::vector<int>       found_kind;             /* 2: complete for I, 1: closed on a subset */
    slong      nodes, edges, open_subset;
    int        truncated;

    slong      pairs_done, anchors_done;
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
    printf("  [%6.1fs] anchors %ld/%ld, pairs %ld/%ld, exact checks %ld (%ld accepted), "
           "edges %ld, clique nodes %ld, families %ld\n",
           t - D.t0, (long) D.anchors_done, (long) D.M,
           (long) D.pairs_done, (long) (D.M * (D.M + 1) / 2),
           (long) D.checked, (long) D.accepted, (long) D.edges, (long) D.nodes,
           (long) D.found.size());
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

/* Exact test of a numerically found n for the pair (a, b):
   p = S n with p_i^2 = q_a q_b (m_i <= 1) or p_i^2 <= q_a q_b (m_i = 2). */
static int dis_verify(dis_t& D, slong a, slong b, const std::vector<slong>& nv, pair_opt_t& out)
{
    search_t* e = D.e;
    const slong n = D.n;
    D.checked++;
    nf_vec P(n, e->field);
    dis_Sn(P, D, nv);
    nf_elem_t u, w; nf_elem_init(u, e->field); nf_elem_init(w, e->field);
    out.n = nv;
    out.sg.assign(n, 0);
    int ok = 1;
    for (slong i = 0; i < n && ok; i++) {
        nf_elem_mul(u, P[i], P[i], e->field);
        nf_elem_mul(w, D.q + a * n + i, D.q + b * n + i, e->field);
        if (e->mult[i] == 2) {
            nf_elem_sub(w, w, u, e->field);
            ok = (dis_sign(w, D) >= 0);
        } else {
            ok = nf_elem_equal(u, w, e->field);
        }
        if (ok) out.sg[i] = (signed char) dis_sign(P[i], D);
    }
    nf_elem_clear(w, e->field); nf_elem_clear(u, e->field);
    D.accepted += ok;
    return ok;
}

/* Number of nonnegative integer vectors of length n with sum of squares c (as a double). */
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
    return w[c];
}

/* Strategy 1 (all m_i = 1 on the common support): choose the signs of
   p_i = +-|b_a^i||b_b^i|. */
static void dis_sign_dfs(dis_t& D, slong a, slong b, const std::vector<slong>& C,
                         const std::vector<double>& rr, const std::vector<double>& tail,
                         size_t d, std::vector<double>& v, std::vector<pair_opt_t>& out)
{
    const slong n = D.n;
    if (d == C.size()) {
        std::vector<slong> nv(n);
        for (slong j = 0; j < n; j++) {
            double k = std::nearbyint(v[j]);
            if (std::fabs(v[j] - k) > DIS_TOL || k < 0) return;
            nv[j] = (slong) k;
        }
        if (nv[0] != 0) return;                        /* distinct boundaries */
        pair_opt_t o;
        if (dis_verify(D, a, b, nv, o)) out.push_back(o);
        return;
    }
    for (slong j = 0; j < n; j++) {                    /* can column j still land on an integer >= 0? */
        double lo = v[j] - tail[d * n + j], hi = v[j] + tail[d * n + j];
        if (hi < -DIS_TOL) return;
        if (std::floor(hi + DIS_TOL) < std::ceil(lo - DIS_TOL)) return;
        if (j == 0 && (lo > DIS_TOL || hi < -DIS_TOL)) return;   /* n^0 = 0 */
    }
    const slong i = C[d];
    for (int eps = 1; eps >= -1; eps -= 2) {
        for (slong j = 0; j < n; j++) v[j] += eps * rr[d] * D.Sd[i * n + j];
        dis_sign_dfs(D, a, b, C, rr, tail, d + 1, v, out);
        for (slong j = 0; j < n; j++) v[j] -= eps * rr[d] * D.Sd[i * n + j];
    }
}

/* Strategy 2: n >= 0 with n^0 = 0 and |n|^2 = x_a.x_b (all m_i = 1) or <= x_a.x_b.
   p = S n must reach +-|b_a^i||b_b^i| (m_i <= 1) or stay within it (m_i = 2);
   coordinates not yet chosen move p_i by at most sqrt(budget) * Qs_i. */
static void dis_norm_dfs(dis_t& D, slong a, slong b, int exact_norm, slong j, slong budget,
                         std::vector<slong>& nv, std::vector<double>& p,
                         std::vector<pair_opt_t>& out)
{
    const slong n = D.n;
    const std::vector<int>& mult = D.e->mult;
    if (budget > 0 && j < n) {
        const double sb = std::sqrt((double) budget);
        for (slong i = 0; i < n; i++) {
            double target = D.r[a * n + i] * D.r[b * n + i];
            double rad = sb * D.Qs[j * n + i] + DIS_TOL;
            if (mult[i] == 2) { if (std::fabs(p[i]) - rad > target) return; }
            else if (std::fabs(p[i] - target) > rad && std::fabs(p[i] + target) > rad) return;
        }
    }
    if (j == n || budget == 0) {
        if (exact_norm && budget != 0) return;
        for (slong i = 0; i < n; i++) {
            double target = D.r[a * n + i] * D.r[b * n + i];
            if (mult[i] == 2) { if (std::fabs(p[i]) > target + DIS_TOL) return; }
            else if (std::fabs(std::fabs(p[i]) - target) > DIS_TOL) return;
        }
        std::vector<slong> full(nv);
        std::fill(full.begin() + j, full.end(), 0);
        pair_opt_t o;
        if (dis_verify(D, a, b, full, o)) out.push_back(o);
        if (budget == 0 || j == n) return;
    }
    for (slong val = 0; val * val <= budget; val++) {
        nv[j] = val;
        if (val) for (slong i = 0; i < n; i++) p[i] += val * D.Sd[i * n + j];
        dis_norm_dfs(D, a, b, exact_norm, j + 1, budget - val * val, nv, p, out);
        if (val) for (slong i = 0; i < n; i++) p[i] -= val * D.Sd[i * n + j];
    }
    nv[j] = 0;
}

/* All verified n_ab between distinct states (n^0 = 0), cached; n_ab = n_ba. */
static const std::vector<pair_opt_t>& dis_options(dis_t& D, slong a, slong b)
{
    if (a > b) std::swap(a, b);
    const slong n = D.n, idx = a * D.M + b;
    if (D.have[idx]) return D.opts[idx];
    D.have[idx] = 1;
    D.pairs_done++;
    dis_progress(D, 0);
    std::vector<pair_opt_t>& out = D.opts[idx];

    const std::vector<slong>& xa = D.e->sols[a];
    const std::vector<slong>& xb = D.e->sols[b];
    slong c = 0;
    for (slong j = 0; j < n; j++) c += xa[j] * xb[j];

    std::vector<slong> C;
    int any2 = 0;
    for (slong i = 0; i < n; i++)
        if (D.supp[a * n + i] && D.supp[b * n + i]) { C.push_back(i); any2 |= (D.e->mult[i] == 2); }

    if (C.empty()) {                                   /* then p = 0, so n = 0 */
        pair_opt_t o; o.n.assign(n, 0); o.sg.assign(n, 0);
        out.push_back(o);
        return out;
    }

    if (any2 || count_sq(c, n) < std::ldexp(1.0, (int) FLINT_MIN((slong) C.size(), 1000))) {
        std::vector<slong> nv(n, 0);
        std::vector<double> p(n, 0.0);
        dis_norm_dfs(D, a, b, !any2, 1, c, nv, p, out);    /* starts at j = 1: n^0 = 0 */
    } else {
        std::sort(C.begin(), C.end(), [&](slong u, slong w) {
            return D.r[a*n+u] * D.r[b*n+u] > D.r[a*n+w] * D.r[b*n+w];
        });
        std::vector<double> rr(C.size()), tail((C.size() + 1) * n, 0.0), v(n, 0.0);
        for (size_t t = 0; t < C.size(); t++) rr[t] = D.r[a * n + C[t]] * D.r[b * n + C[t]];
        for (slong d = (slong) C.size() - 1; d >= 0; d--)
            for (slong j = 0; j < n; j++)
                tail[d * n + j] = tail[(d + 1) * n + j] + std::fabs(D.Sd[C[d] * n + j]) * rr[d];
        dis_sign_dfs(D, a, b, C, rr, tail, 0, v, out);
    }
    return out;
}

/* The anchor itself: b_A^i = (|b_A^i|, 0). */
static vtx_t dis_anchor_vertex(dis_t& D, slong A)
{
    const slong n = D.n;
    vtx_t a;
    a.cand = A;
    a.nA = D.e->sols[A];
    a.s2.assign(n, 0);
    a.w1.assign(n, 0.0); a.w2.assign(n, 0.0);
    a.p = nf_vec(n, D.e->field);
    a.res = nf_vec(n, D.e->field);
    for (slong i = 0; i < n; i++) {
        nf_elem_set(a.p[i], D.q + A * n + i, D.e->field);
        nf_elem_zero(a.res[i], D.e->field);
        a.w1[i] = D.r[A * n + i];
    }
    return a;
}

/* States of candidate B with overlap n_AB with the anchor: one per sign choice of the
   second component on labels with m_i = 2 and positive residual. Labels outside
   supp A get zero components (no member of the family is supported there). */
static void dis_make_vertices(dis_t& D, slong A, slong B, const pair_opt_t& o, std::vector<vtx_t>& out)
{
    search_t* e = D.e;
    const slong n = D.n;
    nf_vec P(n, e->field), R(n, e->field);
    dis_Sn(P, D, o.n);
    nf_elem_t t; nf_elem_init(t, e->field);
    std::vector<slong> two;                          /* labels with a free second-component sign */
    arb_t g, v; arb_init(g); arb_init(v);
    qqbar_get_arb(g, e->gen, NUM_PREC);
    std::vector<double> w1(n, 0.0), rabs(n, 0.0);
    for (slong i = 0; i < n; i++) {
        nf_elem_mul(R[i], D.q + A * n + i, D.q + B * n + i, e->field);
        nf_elem_mul(t, P[i], P[i], e->field);
        nf_elem_sub(R[i], R[i], t, e->field);
        if (!D.supp[A * n + i]) continue;
        nf_elem_get_arb(v, P[i], g, e->field, NUM_PREC);
        w1[i] = arf_get_d(arb_midref(v), ARF_RND_NEAR) / D.r[A * n + i];
        if (e->mult[i] == 2 && dis_sign(R[i], D) > 0) {
            two.push_back(i);
            nf_elem_get_arb(v, R[i], g, e->field, NUM_PREC);
            arb_sqrt(v, v, NUM_PREC);
            rabs[i] = arf_get_d(arb_midref(v), ARF_RND_NEAR) / D.r[A * n + i];
        }
    }
    arb_clear(v); arb_clear(g);
    nf_elem_clear(t, e->field);

    for (ulong mask = 0; mask < (1ul << two.size()); mask++) {
        vtx_t x;
        x.cand = B;
        x.nA = o.n;
        x.s2.assign(n, 0);
        x.w1 = w1;
        x.w2.assign(n, 0.0);
        for (size_t k = 0; k < two.size(); k++) {
            slong i = two[k];
            x.s2[i] = ((mask >> k) & 1) ? -1 : 1;
            x.w2[i] = x.s2[i] * rabs[i];
        }
        x.p = nf_vec(n, e->field);
        x.res = nf_vec(n, e->field);
        for (slong i = 0; i < n; i++) {
            nf_elem_set(x.p[i], P[i], e->field);
            nf_elem_set(x.res[i], R[i], e->field);
        }
        out.push_back(std::move(x));
    }
}

/* Cardy's condition between two states in the gauge of anchor A:
   cheap numeric tests, then the exact identities above. */
static int dis_compatible(dis_t& D, slong A, const vtx_t& u, const vtx_t& v)
{
    search_t* e = D.e;
    const slong n = D.n;
    std::vector<double> p(n);
    for (slong i = 0; i < n; i++) p[i] = u.w1[i] * v.w1[i] + u.w2[i] * v.w2[i];
    std::vector<slong> nv(n);
    for (slong j = 0; j < n; j++) {
        double s = 0.0;
        for (slong i = 0; i < n; i++) s += D.Sd[i * n + j] * p[i];
        double k = std::nearbyint(s);
        if (k < 0 || std::fabs(s - k) > DIS_TOL) return 0;
        if (j == 0 && k != 0) return 0;               /* distinct elementary boundaries */
        nv[j] = (slong) k;
    }
    D.checked++;
    nf_vec P(n, e->field);
    dis_Sn(P, D, nv);
    nf_elem_t X, Y, Z; nf_elem_init(X, e->field); nf_elem_init(Y, e->field); nf_elem_init(Z, e->field);
    int ok = 1;
    for (slong i = 0; i < n && ok; i++) {
        if (!D.supp[A * n + i]) { ok = nf_elem_is_zero(P[i], e->field); continue; }
        nf_elem_mul(X, D.q + A * n + i, P[i], e->field);          /* q_A P - p_AB p_AC */
        nf_elem_mul(Y, u.p[i], v.p[i], e->field);
        nf_elem_sub(X, X, Y, e->field);
        nf_elem_mul(Y, u.res[i], v.res[i], e->field);             /* residual product */
        nf_elem_mul(Z, X, X, e->field);
        if (!nf_elem_equal(Z, Y, e->field)) { ok = 0; break; }
        if (!nf_elem_is_zero(Y, e->field)) ok = (dis_sign(X, D) == u.s2[i] * v.s2[i]);
    }
    nf_elem_clear(Z, e->field); nf_elem_clear(Y, e->field); nf_elem_clear(X, e->field);
    D.accepted += ok;
    return ok;
}

static nmat_t dis_nmatrix(dis_t& D, const std::vector<member_t>& F)
{
    const slong n = D.n;
    nmat_t N(F.size(), std::vector<std::vector<slong>>(F.size(), std::vector<slong>(n)));
    for (size_t k = 0; k < F.size(); k++)
        for (size_t l = 0; l < F.size(); l++)
            for (slong j = 0; j < n; j++) {
                double v = 0.0;
                for (slong i = 0; i < n; i++)
                    v += D.Sd[i * n + j] * (F[k].w1[i] * F[l].w1[i] + F[k].w2[i] * F[l].w2[i]);
                N[k][l][j] = (slong) std::nearbyint(v);   /* exact integers, already verified */
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
        if (used[g] || G[g].cand != F[k].cand || NG[g][g] != NF[k][k]) continue;
        int ok = 1;
        for (size_t l = 0; l < k && ok; l++) ok = (NF[k][l] == NG[g][pi[l]]);
        if (!ok) continue;
        used[g] = 1; pi[k] = (slong) g;
        if (dis_iso_rec(F, NF, G, NG, pi, used, k + 1)) return 1;
        used[g] = 0;
    }
    return 0;
}

typedef std::vector<uint64_t> bits_t;

static slong bits_count(const bits_t& b)
{
    slong c = 0;
    for (uint64_t w : b) c += __builtin_popcountll(w);
    return c;
}

/* All cliques of size T: vertices taken in increasing order, |R| + |cand| >= T. */
static void dis_cliques(dis_t& D, const std::vector<bits_t>& adj, std::vector<slong>& R,
                        bits_t cand, slong T, std::vector<std::vector<slong>>& out)
{
    if (D.truncated) return;
    if (++D.nodes > CLIQUE_NODE_LIMIT) { D.truncated = 1; return; }
    if ((D.nodes & 0xFFF) == 0) dis_progress(D, 0);
    if ((slong) R.size() == T) { out.push_back(R); return; }
    slong have = bits_count(cand);
    for (size_t w = 0; w < cand.size() && (slong) R.size() + have >= T; w++) {
        while (cand[w] && (slong) R.size() + have >= T) {
            slong v = (slong) (w * 64 + __builtin_ctzll(cand[w]));
            cand[w] &= cand[w] - 1;
            have--;
            bits_t c2(cand.size());
            for (size_t k = 0; k < cand.size(); k++) c2[k] = cand[k] & adj[v][k];
            R.push_back(v);
            dis_cliques(D, adj, R, c2, T, out);
            R.pop_back();
        }
    }
}

static member_t dis_member(const vtx_t& v, slong A)
{
    return member_t{v.cand, A, v.nA, v.s2, v.w1, v.w2};
}

/* Store a complete family (dropping gauge-equivalent repeats) and announce it. */
static void dis_add_family(dis_t& D, const std::vector<member_t>& F, int kind)
{
    nmat_t N = dis_nmatrix(D, F);
    for (size_t f = 0; f < D.found.size(); f++) {
        if (D.found[f].size() != F.size()) continue;
        std::vector<slong> pi(F.size());
        std::vector<char> used(F.size(), 0);
        if (dis_iso_rec(F, N, D.found[f], D.found_n[f], pi, used, 0)) return;
    }
    D.found.push_back(F);
    D.found_n.push_back(N);
    D.found_kind.push_back(kind);
    printf("  [%6.1fs] family %ld: %ld boundaries, %s\n", dis_now() - D.t0,
           (long) D.found.size(), (long) F.size(),
           kind == 2 ? "complete for I" : "complete and closed on a subset of I");
    fflush(stdout);
}

/* Complete families anchored at candidate A. */
static void dis_anchor(dis_t& D, slong A)
{
    const slong n = D.n, M = D.M;
    slong N = 0;                                      /* states in a basis over supp A */
    for (slong i = 0; i < n; i++) if (D.supp[A * n + i]) N += D.e->mult[i];

    std::vector<vtx_t> V;
    V.push_back(dis_anchor_vertex(D, A));
    for (slong B = 0; B < M; B++) {
        int inside = 1, equal = 1;
        for (slong i = 0; i < n; i++) {
            if (D.supp[B * n + i] && !D.supp[A * n + i]) inside = 0;
            if (D.supp[B * n + i] != D.supp[A * n + i]) equal = 0;
        }
        if (!inside || (equal && B < A)) continue;    /* such families have an earlier anchor */
        for (const pair_opt_t& o : dis_options(D, A, B)) dis_make_vertices(D, A, B, o, V);
    }
    if ((slong) V.size() < N) return;                 /* too few states for a basis */

    /* graph on V[1..]; every vertex is already compatible with the anchor */
    const slong nv = (slong) V.size() - 1, nw = (nv + 63) / 64, T = N - 1;
    std::vector<bits_t> adj(nv, bits_t(nw, 0));
    std::vector<slong> deg(nv, 0);
    for (slong k = 0; k < nv; k++)
        for (slong l = k + 1; l < nv; l++)
            if (dis_compatible(D, A, V[k + 1], V[l + 1])) {
                adj[k][l / 64] |= 1ull << (l % 64);
                adj[l][k / 64] |= 1ull << (k % 64);
                deg[k]++; deg[l]++; D.edges++;
            }

    /* k-core: a member of a T-clique has T - 1 neighbours inside it */
    bits_t alive(nw, 0);
    for (slong k = 0; k < nv; k++) alive[k / 64] |= 1ull << (k % 64);
    for (int changed = 1; changed; ) {
        changed = 0;
        for (slong k = 0; k < nv; k++) {
            if (!((alive[k / 64] >> (k % 64)) & 1) || deg[k] >= T - 1) continue;
            alive[k / 64] &= ~(1ull << (k % 64));
            for (slong l = 0; l < nv; l++)
                if ((adj[k][l / 64] >> (l % 64)) & 1) deg[l]--;
            changed = 1;
        }
    }
    if (bits_count(alive) < T) return;

    std::vector<std::vector<slong>> cliques;
    std::vector<slong> R;
    dis_cliques(D, adj, R, alive, T, cliques);

    for (const std::vector<slong>& c : cliques) {
        std::vector<member_t> F;
        F.push_back(dis_member(V[0], A));
        for (slong k : c) F.push_back(dis_member(V[k + 1], A));
        int kind = (N == D.e->dimI) ? 2 : 1;
        if (kind == 1) {
            /* closed: no other state, of any support, is compatible with all of F */
            int closed = 1;
            for (slong X = 0; X < M && closed; X++)
                for (const pair_opt_t& o : dis_options(D, A, X)) {
                    std::vector<vtx_t> xs;
                    dis_make_vertices(D, A, X, o, xs);
                    for (const vtx_t& x : xs) {
                        int member = 0;
                        for (slong k : c)
                            if (V[k + 1].cand == X && V[k + 1].nA == x.nA && V[k + 1].s2 == x.s2) member = 1;
                        if (member) continue;
                        int fits = 1;
                        for (size_t k = 0; k < c.size() && fits; k++) fits = dis_compatible(D, A, x, V[c[k] + 1]);
                        if (fits) { closed = 0; break; }
                    }
                    if (!closed) break;
                }
            if (!closed) { D.open_subset++; continue; }
        }
        dis_add_family(D, F, kind);
    }
}

static void disentangle(search_t* e)
{
    const slong n = e->n, M = (slong) e->sols.size();
    flint_printf("\n--- disentangling %wd candidates with the mutual Cardy condition ---\n", M);
    if (M == 0) return;
    for (slong i = 0; i < n; i++)
        if (e->mult[i] > 2) {
            flint_printf("label %wd is listed %d times; disentangling supports at most two "
                         "Ishibashi states per label, skipped\n", i + INDEX_BASE, e->mult[i]);
            return;
        }

    dis_t D;
    D.e = e; D.n = n; D.M = M;
    D.checked = D.accepted = D.nodes = D.edges = D.open_subset = 0;
    D.truncated = 0;
    D.pairs_done = D.anchors_done = 0;
    D.t0 = D.t_last = dis_now();
    D.opts.assign(M * M, std::vector<pair_opt_t>());
    D.have.assign(M * M, 0);

    /* exact |b_a^i|^2 for every i (it is exactly 0 for i not in I) */
    D.q = (nf_elem_struct*) flint_malloc(M * n * sizeof(nf_elem_struct));
    for (slong k = 0; k < M * n; k++) nf_elem_init(D.q + k, e->field);
    D.r.assign(M * n, 0.0);
    D.supp.assign(M * n, 0);
    arb_t g, v; arb_init(g); arb_init(v);
    qqbar_get_arb(g, e->gen, NUM_PREC);
    {
        nf_vec P(n, e->field);
        for (slong a = 0; a < M; a++) {
            dis_Sn(P, D, e->sols[a]);
            for (slong i = 0; i < n; i++) {
                nf_elem_set(D.q + a * n + i, P[i], e->field);
                if (nf_elem_is_zero(P[i], e->field)) continue;
                D.supp[a * n + i] = 1;
                nf_elem_get_arb(v, P[i], g, e->field, NUM_PREC);
                arb_sqrt(v, v, NUM_PREC);
                D.r[a * n + i] = arf_get_d(arb_midref(v), ARF_RND_NEAR);
            }
        }
    }
    D.Sd.assign(n * n, 0.0);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) {
            nf_elem_get_arb(v, entry(e->S, i, j, e->ctx), g, e->field, NUM_PREC);
            D.Sd[i * n + j] = arf_get_d(arb_midref(v), ARF_RND_NEAR);
        }
    D.Qs.assign((n + 1) * n, 0.0);
    for (slong j = n - 1; j >= 0; j--)
        for (slong i = 0; i < n; i++)
            D.Qs[j * n + i] = std::sqrt(D.Qs[(j + 1) * n + i] * D.Qs[(j + 1) * n + i]
                                        + D.Sd[i * n + j] * D.Sd[i * n + j]);
    arb_clear(v); arb_clear(g);

    printf("  [%6.1fs] |b|^2 computed exactly for %ld candidates, %ld Ishibashi states "
           "over %ld labels, searching families\n",
           dis_now() - D.t0, (long) M, (long) e->dimI, (long) e->nI);
    fflush(stdout);
    for (slong A = 0; A < M && !D.truncated; A++) {
        dis_anchor(D, A);
        D.anchors_done++;
        dis_progress(D, 0);
    }
    dis_progress(D, 1);

    /* report: complete for I first, then closed subset families, larger first */
    std::vector<int>&  kind = D.found_kind;
    std::vector<slong> ncov(D.found.size(), 0), in_complete(M, 0);
    for (size_t f = 0; f < D.found.size(); f++) {
        for (slong i = 0; i < n; i++)
            if (D.supp[D.found[f][0].anchor * n + i]) ncov[f] += e->mult[i];
        if (kind[f] == 2) for (const member_t& m : D.found[f]) in_complete[m.cand]++;
    }
    std::vector<size_t> ord(D.found.size());
    for (size_t f = 0; f < ord.size(); f++) ord[f] = f;
    std::stable_sort(ord.begin(), ord.end(), [&](size_t u, size_t w) {
        if (kind[u] != kind[w]) return kind[u] > kind[w];
        return D.found[u].size() > D.found[w].size();
    });
    slong nk[3] = {0, 0, 0};
    for (size_t f = 0; f < D.found.size(); f++) nk[kind[f]]++;

    flint_printf("pair checks: %wd exact, %wd accepted; compatibility edges: %wd; clique nodes: %wd%s\n",
                 D.checked, D.accepted, D.edges, D.nodes,
                 D.truncated ? "  (CLIQUE NODE LIMIT HIT: list below may be incomplete)" : "");
    flint_printf("complete families (up to gauge): %wd for I, %wd closed on a subset of I"
                 " (%wd more on a subset but extendable, not listed)\n",
                 nk[2], nk[1], D.open_subset);
    flint_printf("(searched: families containing a member whose support covers all the others')\n");

    for (size_t c = 0; c < ord.size() && c < (size_t) MAX_FAMILIES_SHOWN; c++) {
        const std::vector<member_t>& F = D.found[ord[c]];
        if (kind[ord[c]] == 2)
            flint_printf("\nfamily %wd: %wd boundaries, complete for I\n", (slong) c + 1, (slong) F.size());
        else
            flint_printf("\nfamily %wd: %wd boundaries, complete on %wd of the %wd Ishibashi states\n",
                         (slong) c + 1, (slong) F.size(), ncov[ord[c]], e->dimI);
        for (const member_t& m : F) {
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
       with one Ishibashi state, (b^{i,1}, b^{i,2}) for labels with two, in the basis
       where the family's first boundary (its anchor) lies along the first axis. */
    int has2 = 0;
    for (slong i = 0; i < n; i++) has2 |= (e->mult[i] == 2);
    arb_t gb, c1, c2, t1, bnum, bden;
    arb_init(gb); arb_init(c1); arb_init(c2); arb_init(t1); arb_init(bnum); arb_init(bden);
    qqbar_get_arb(gb, e->gen, NUM_PREC);
    const slong inum = TRANS_NUM - INDEX_BASE, iden = TRANS_DEN - INDEX_BASE;
    const int have_T = inum >= 0 && inum < n && iden >= 0 && iden < n
                       && e->inI[inum] && e->inI[iden];
    auto absb = [&](arb_t out, slong cand, slong i) {
        nf_elem_get_arb(out, D.q + cand * n + i, gb, e->field, NUM_PREC);
        arb_sqrt(out, out, NUM_PREC);
    };

    flint_printf("\n=== boundaries ===\n");
    if (has2)
        flint_printf("labels with two Ishibashi states are printed as (b^{i,1}, b^{i,2}), in the basis\n"
                     "where the first boundary of the family lies along the first axis\n");
    if (nk[2] + nk[1] == 0) flint_printf("no complete family found\n");
    if (nk[2] == 0 && e->dimI == e->nI)
        flint_printf("note: every label was given one Ishibashi state. For a non-diagonal modular\n"
                     "invariant, list each label as many times as it has Ishibashi states in the\n"
                     "Ishibashi index input (e.g. sigma and psi twice for Potts under Virasoro)\n");
    slong shown = 0;
    nf_vec P(n, e->field), R(n, e->field);
    nf_elem_t t; nf_elem_init(t, e->field);
    for (size_t c = 0; c < ord.size(); c++) {
        if (shown++ >= MAX_FAMILIES_SHOWN) {
            flint_printf("... further complete families not shown\n");
            break;
        }
        const std::vector<member_t>& F = D.found[ord[c]];
        if (kind[ord[c]] == 2)
            flint_printf("\nfamily %wd (complete for I): %wd boundaries\n", (slong) c + 1, (slong) F.size());
        else
            flint_printf("\nfamily %wd (complete on %wd of the %wd Ishibashi states): %wd boundaries\n",
                         (slong) c + 1, ncov[ord[c]], e->dimI, (slong) F.size());
        for (size_t k = 0; k < F.size(); k++) {
            const member_t& m = F[k];
            const slong A = m.anchor;
            dis_Sn(P, D, m.nA);                       /* p_A,m; for the anchor, q_A */
            flint_printf("boundary %wd.%wd (candidate %wd)\nx = ", (slong) c + 1, (slong) k + 1, m.cand + 1);
            for (slong j = 0; j < n; j++) flint_printf("%wd ", e->sols[m.cand][j]);
            flint_printf("\n");
            for (slong i = 0; i < n; i++) {
                if (!e->inI[i]) continue;
                if (!D.supp[m.cand * n + i]) { arb_zero(c1); arb_zero(c2); }
                else {
                    absb(t1, A, i);                                      /* |b_A^i| */
                    nf_elem_get_arb(c1, P[i], gb, e->field, NUM_PREC);
                    arb_div(c1, c1, t1, NUM_PREC);
                    nf_elem_mul(R[i], D.q + A * n + i, D.q + m.cand * n + i, e->field);
                    nf_elem_mul(t, P[i], P[i], e->field);
                    nf_elem_sub(R[i], R[i], t, e->field);
                    if (nf_elem_is_zero(R[i], e->field)) arb_zero(c2);
                    else {
                        nf_elem_get_arb(c2, R[i], gb, e->field, NUM_PREC);
                        arb_sqrt(c2, c2, NUM_PREC);
                        arb_div(c2, c2, t1, NUM_PREC);
                        if (m.s2[i] < 0) arb_neg(c2, c2);
                    }
                }
                flint_printf("  b^%wd = ", i + INDEX_BASE);
                if (e->mult[i] == 2) {
                    flint_printf("(");
                    arb_printn(c1, 15, 0);
                    flint_printf(", ");
                    arb_printn(c2, 15, 0);
                    flint_printf(")");
                } else {
                    arb_printn(c1, 15, 0);
                }
                flint_printf("\n");
            }
            if (have_T) {
                flint_printf("  T = 1 - |b^%d|/|b^%d| = ", TRANS_NUM, TRANS_DEN);
                if (!D.supp[m.cand * n + iden]) flint_printf("undefined (b^%d = 0)\n", TRANS_DEN);
                else {
                    absb(bnum, m.cand, inum);
                    absb(bden, m.cand, iden);
                    arb_div(bnum, bnum, bden, NUM_PREC);
                    arb_neg(bnum, bnum);
                    arb_add_ui(bnum, bnum, 1, NUM_PREC);
                    arb_printn(bnum, 15, 0);
                    flint_printf("\n");
                }
            }
        }
    }
    nf_elem_clear(t, e->field);
    arb_clear(bden); arb_clear(bnum); arb_clear(t1); arb_clear(c2); arb_clear(c1); arb_clear(gb);
    fflush(stdout);

    for (slong k = 0; k < M * n; k++) nf_elem_clear(D.q + k, e->field);
    flint_free(D.q);
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
