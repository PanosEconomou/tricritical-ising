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
#include <cmath>
#include <cstdio>
#include <vector>

#include "boundary_solver.h"

#define BOUND_PREC 256
#define NUM_PREC   128
#define LP_TOL     1e-6   /* LP ranges are widened by this before rounding */
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
    std::vector<slong> lb, ub;
    std::vector<slong> x;       /* current assignment */

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

    slong              nodes, leaves, nsolutions;
} search_t;

static void set_col_box(glp_prob* lp, int col, slong lo, slong hi)
{
    glp_set_col_bnds(lp, col, lo == hi ? GLP_FX : GLP_DB, (double) lo, (double) hi);
}

/*
 * Range of column `col` over the LP relaxation with the current fixings.
 * Returns 0 if the LP is infeasible. On numerical trouble it falls back to
 * the column's box, so a true solution is never lost.
 */
static int lp_range(search_t* e, int col, double* lo, double* hi)
{
    int feasible = 1;
    glp_set_obj_coef(e->lp, col, 1.0);

    for (int dir = 0; dir < 2 && feasible; dir++) {
        glp_set_obj_dir(e->lp, dir == 0 ? GLP_MIN : GLP_MAX);

        int ret = glp_simplex(e->lp, &e->parm);
        if (ret != 0) { glp_std_basis(e->lp); ret = glp_simplex(e->lp, &e->parm); }

        double val = (dir == 0) ? (double) e->lb[col - 1] : (double) e->ub[col - 1];
        if (ret == 0) {
            int st = glp_get_status(e->lp);
            if (st == GLP_OPT)         val = glp_get_obj_val(e->lp);
            else if (st == GLP_NOFEAS) feasible = 0;
        }
        if (dir == 0) *lo = val; else *hi = val;
    }

    glp_set_obj_coef(e->lp, col, 0.0);
    return feasible;
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
    double lo, hi;
    if (!lp_range(e, col, &lo, &hi)) return;

    slong a = FLINT_MAX(e->lb[j], (slong) std::ceil(lo - LP_TOL));
    slong b = FLINT_MIN(e->ub[j], (slong) std::floor(hi + LP_TOL));

    for (slong val = a; val <= b; val++) {
        glp_set_col_bnds(e->lp, col, GLP_FX, (double) val, (double) val);
        e->x[j] = val;
        dfs(e, d + 1);
    }
    set_col_box(e->lp, col, e->lb[j], e->ub[j]);
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
 * columns, taking each one's range from an LP over all constraints (floating
 * point, widened by LP_TOL), and recovers the pivot columns exactly at the
 * leaves. Strict positivity is then checked exactly.
 * Returns the number of solutions, or -1 on error.
 */
int solve_boundary_bootstrap(const slong* indices, slong nindices,
                             gr_mat_struct* matrix, nf_struct* field,
                             qqbar_struct* field_generator, gr_ctx_struct* context)
{
    const slong n      = gr_mat_nrows(matrix, context);
    const slong degree = fmpq_poly_degree(field->pol);

    if (gr_mat_ncols(matrix, context) != n || nindices <= 0) return -1;

    /* I as a mask */
    std::vector<char> inI(n, 0);
    for (slong k = 0; k < nindices; k++) {
        slong idx = indices[k] - INDEX_BASE;
        if (idx < 0 || idx >= n) {
            flint_fprintf(stderr, "ERROR: Ishibashi index %wd out of range\n", indices[k]);
            return -1;
        }
        inI[idx] = 1;
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
    e.n  = n;
    e.Z  = Z;
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

        /* 5. LP:  Z x = 0,  sum_j S_ij x^j >= -LP_TOL (i in I),  lb <= x <= ub */
        slong nI = n - nout;
        glp_term_out(GLP_OFF);
        e.lp = glp_create_prob();
        glp_add_cols(e.lp, (int) n);
        for (slong j = 0; j < n; j++) set_col_box(e.lp, (int) j + 1, e.lb[j], e.ub[j]);
        glp_add_rows(e.lp, (int) (rank + nI));

        std::vector<int>    ia(1, 0), ja(1, 0);
        std::vector<double> ar(1, 0.0);

        for (slong r = 0; r < rank; r++) {
            glp_set_row_bnds(e.lp, (int) r + 1, GLP_FX, 0.0, 0.0);
            for (slong j = 0; j < n; j++) {
                if (fmpz_is_zero(fmpz_mat_entry(Z, r, j))) continue;
                ia.push_back((int) r + 1); ja.push_back((int) j + 1);
                ar.push_back(fmpz_get_d(fmpz_mat_entry(Z, r, j)));
            }
        }
        for (slong i = 0, k = rank; i < n; i++) {
            if (!inI[i]) continue;
            glp_set_row_bnds(e.lp, (int) k + 1, GLP_LO, -LP_TOL, 0.0);
            for (slong j = 0; j < n; j++) {
                ia.push_back((int) k + 1); ja.push_back((int) j + 1);
                ar.push_back(arf_get_d(arb_midref(Sa + i * n + j), ARF_RND_NEAR));
            }
            k++;
        }
        glp_load_matrix(e.lp, (int) ia.size() - 1, ia.data(), ja.data(), ar.data());
        glp_scale_prob(e.lp, GLP_SF_AUTO);

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
        status = (int) e.nsolutions;
    } while (0);

    if (e.lp) glp_delete_prob(e.lp);
    arb_clear(e.v);
    fmpz_clear(e.acc);
    arb_clear(g);
    _arb_vec_clear(Sa, n * n);
    fmpz_mat_clear(Z);
    return status;
}
