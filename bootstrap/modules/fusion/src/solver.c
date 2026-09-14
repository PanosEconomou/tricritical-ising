#include <gmp.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>
#include <flint/fmpq_vec.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>
#include <flint/arb.h>
#include <flint/arb_fmpz_poly.h>
#include <4ti2/4ti2.h>
#include <stdio.h>
#include <time.h>

#include "solver.h"

#define BOUND_PREC 256

static int build_rhs(fmpq_mat_t rhs, const gr_mat_t matrix, slong i, slong j,
                       const nf_t field, gr_ctx_t context)
{
    slong rows   = gr_mat_nrows(matrix, context); 
    slong degree = fmpq_poly_degree(field->pol);
    if (fmpq_mat_nrows(rhs)!= rows * degree || fmpq_mat_ncols(rhs) != 1) return -1;

    nf_elem_t   product;        nf_elem_init(product, field);
    fmpq_poly_t coeffiecients;  fmpq_poly_init(coeffiecients);

    for (slong k = 0; k < rows; k++) {
        nf_elem_mul(product, 
                    gr_mat_entry_srcptr(matrix, k, i, context),
                    gr_mat_entry_srcptr(matrix, k, j, context),
                    field
                    );

        nf_elem_get_fmpq_poly(coeffiecients, product, field);

        for (slong l = 0; l < degree; l++) {
            fmpq_poly_get_coeff_fmpq(
                fmpq_mat_entry(rhs, degree * k + l, 0), 
                coeffiecients,
                l
            );
        }
    }

    fmpq_poly_clear(coeffiecients);
    nf_elem_clear(product, field);

    return 0;
}

static int build_lhs(fmpq_mat_t lhs, const gr_mat_t matrix, 
                     const nf_t field, gr_ctx_t context) 
{
    slong rows   = gr_mat_nrows(matrix, context);
    slong cols   = gr_mat_ncols(matrix, context);
    slong degree = fmpq_poly_degree(field->pol);
    fmpq_poly_t coefficients; fmpq_poly_init(coefficients);

    if (fmpq_mat_nrows(lhs)!= rows * degree || fmpq_mat_ncols(lhs) != cols) return -1;
    
    for (slong i = 0; i < rows; i++) {
        for (slong j = 0; j < cols; j++) {
            nf_elem_get_fmpq_poly(coefficients, 
                                  gr_mat_entry_srcptr(matrix, i, j, context), 
                                  field);

            for (slong k = 0; k < degree; k++) {
                fmpq_poly_get_coeff_fmpq(fmpq_mat_entry(lhs, degree * i + k, j), 
                                         coefficients,
                                         k);
            }
        }
    }

    fmpq_poly_clear(coefficients);

    return 0;
}

static slong nf_elem_ceil(const nf_elem_t a, const arb_t generator,
                          const nf_t field, slong prec)
{
    fmpq_poly_t p; fmpq_poly_init(p);
    arb_t v;       arb_init(v);

    nf_elem_get_fmpq_poly(p, a, field);
    _arb_fmpz_poly_evaluate_arb(v, fmpq_poly_numref(p), fmpq_poly_length(p), generator, prec);
    arb_div_fmpz(v, v, fmpq_poly_denref(p), prec);

    slong n = arf_get_si(arb_midref(v), ARF_RND_CEIL) + 1;

    arb_clear(v);
    fmpq_poly_clear(p);
    return n;
}

static int build_bounds(slong* bounds, gr_mat_t matrix, slong i, slong j, slong row, 
                        nf_t field, qqbar_t field_generator, gr_ctx_t context)
{
    slong unknowns = gr_mat_ncols(matrix, context);

    nf_elem_t product;   nf_elem_init(product, field);
    nf_elem_t ratio;     nf_elem_init(ratio,   field);
    arb_t     generator; arb_init(generator);

    qqbar_get_arb(generator, field_generator, BOUND_PREC);

    nf_elem_mul(product,
                gr_mat_entry_srcptr(matrix, row, i, context),
                gr_mat_entry_srcptr(matrix, row, j, context),
                field);

    for (slong k = 0; k < unknowns; k++) {
        nf_elem_div(ratio, 
                    product, 
                    gr_mat_entry_srcptr(matrix, row, k, context),
                    field);
        bounds[k] = nf_elem_ceil(ratio, generator, field, BOUND_PREC);
    }

    arb_clear(generator);
    nf_elem_clear(product, field);
    nf_elem_clear(ratio,   field);

    return 0;
}

// hehe get it? It converts the problem to one with integer coefficients tehee
static int zfy(fmpz_mat_t zrhs, fmpz_mat_t zlhs, fmpq_mat_t rhs, fmpq_mat_t lhs) 
{
    slong rows = fmpq_mat_nrows(lhs);
    slong cols = fmpq_mat_ncols(lhs);

    if (fmpq_mat_nrows(rhs)  != rows || fmpq_mat_ncols(rhs)  != 1)    return -1;
    if (fmpz_mat_nrows(zlhs) != rows || fmpz_mat_ncols(zlhs) != cols) return -1;
    if (fmpz_mat_nrows(zrhs) != rows || fmpz_mat_ncols(zrhs) != 1)    return -1;

    fmpq_mat_get_fmpz_mat_rowwise_2(zlhs, zrhs, NULL, lhs, rhs);

    return 0;
}

static int setup_zsolve(_4ti2_state** state, const fmpz_mat_t zlhs, const fmpz_mat_t zrhs, 
                        const slong* bounds)
{
    int status    = -1;

    int equations = (int) fmpz_mat_nrows(zlhs);
    int unknowns  = (int) fmpz_mat_ncols(zlhs);

    _4ti2_state*  s = _4ti2_zsolve_create_state(_4ti2_PREC_INT_ARB);
    _4ti2_matrix* m = NULL;
    mpz_t tmp; mpz_init(tmp);

    if (fmpz_mat_nrows(zrhs) != equations || fmpz_mat_ncols(zrhs) != 1) goto cleanup;
    if (s == NULL) goto cleanup;

    // Solve Ax = b with certain bounds for x
    // Setup matrix A
    if (_4ti2_state_create_matrix(s, equations, unknowns, "mat", &m) != _4ti2_OK)
        goto cleanup;

    for (int i = 0; i < equations; i++) {
        for (int j = 0; j < unknowns; j++) {
            fmpz_get_mpz(tmp, fmpz_mat_entry(zlhs, i, j));
            if (_4ti2_matrix_set_entry_mpz_ptr(m, i, j, tmp) != _4ti2_OK) 
                goto cleanup;
        }
    }

    // Set up b
    if (_4ti2_state_create_matrix(s, 1, equations, "rhs", &m) != _4ti2_OK)
        goto cleanup;

    for (int i = 0; i < equations; i++) {
        fmpz_get_mpz(tmp, fmpz_mat_entry(zrhs, i, 0));
        if (_4ti2_matrix_set_entry_mpz_ptr(m, 0, i, tmp) != _4ti2_OK) 
            goto cleanup;
    }

    // Set up bounds 
    // sign
    if (_4ti2_state_create_matrix(s, 1, unknowns, "sign", &m) != _4ti2_OK)
        goto cleanup;

    for (int i = 0; i < unknowns; i++) {
        if (_4ti2_matrix_set_entry_int32_t(m, 0, i, _4ti2_DB) != _4ti2_OK) 
            goto cleanup;
    }

    // lowe bound to 0
    if (_4ti2_state_create_matrix(s, 1, unknowns, "lb", &m) != _4ti2_OK)
        goto cleanup;

    for (int i = 0; i < unknowns; i++) {
        if (_4ti2_matrix_set_entry_int32_t(m, 0, i, 0) != _4ti2_OK) 
            goto cleanup;
    }

    // Upper bound by the vector
    if (_4ti2_state_create_matrix(s, 1, unknowns, "ub", &m) != _4ti2_OK)
        goto cleanup;

    for (int i = 0; i < unknowns; i++) {
        if (_4ti2_matrix_set_entry_int32_t(m, 0, i, (int32_t) bounds[i]) != _4ti2_OK) 
            goto cleanup;
    }

    *state = s;
    s      = NULL;
    status = 0;

cleanup:
    mpz_clear(tmp);
    if (s) _4ti2_state_delete(s);
    return status;
}

int solve(gr_mat_struct* matrix, nf_t field, qqbar_t field_generator, gr_ctx_struct *context)
{

    slong rows    =  gr_mat_nrows(matrix, context); 
    slong cols    =  gr_mat_ncols(matrix, context); 
    slong degree  =  fmpq_poly_degree(field->pol);
    fmpq_mat_t lhs;  fmpq_mat_init(lhs, rows * degree, cols);
    fmpq_mat_t rhs;  fmpq_mat_init(rhs, rows * degree, 1);
    fmpz_mat_t zlhs; fmpz_mat_init(zlhs, rows * degree, cols);
    fmpz_mat_t zrhs; fmpz_mat_init(zrhs, rows * degree, 1);
    slong* bounds =  flint_malloc(cols * sizeof(slong));

    build_lhs(lhs, matrix, field, context);
    build_rhs(rhs, matrix, 12, 24, field, context);
    build_bounds(bounds, matrix, 1, 5, 0, field, field_generator, context);
    zfy(zrhs, zlhs, rhs, lhs);

    // fmpq_mat_print(lhs); printf("\n");
    // fmpz_mat_print(zrhs); printf("\n");
    
    _4ti2_state* state = NULL;
    if (setup_zsolve(&state, zlhs, zrhs, bounds)) {
        fprintf(stderr, "ERROR: could not set up zsolve\n");
        goto cleanup;
    }

    clock_t start = clock();
    if (_4ti2_state_compute(state) != _4ti2_OK) {
        fprintf(stderr, "ERROR: zsolve failed\n");
        goto cleanup;
    }

    flint_printf("zsolve took %.2f sec\n",
                 (double)(clock() - start) / CLOCKS_PER_SEC);

    _4ti2_matrix* out;

    if (_4ti2_state_get_matrix(state, "zhom", &out) == _4ti2_OK)
        flint_printf("zhom:   %d rows (expect 0)\n", _4ti2_matrix_get_num_rows(out));

    if (_4ti2_state_get_matrix(state, "zfree", &out) == _4ti2_OK)
        flint_printf("zfree:  %d rows (expect 0)\n", _4ti2_matrix_get_num_rows(out));

    if (_4ti2_state_get_matrix(state, "zinhom", &out) == _4ti2_OK) {
        int nsol = _4ti2_matrix_get_num_rows(out);
        flint_printf("zinhom: %d solutions\n", nsol);
        if (nsol > 0 && nsol <= 50)
            _4ti2_matrix_write_to_stdout(out);
        else if (nsol > 50)
            flint_printf("(too many to print)\n");
    }

cleanup:
    if (state) _4ti2_state_delete(state);
    flint_free(bounds);
    fmpz_mat_clear(zrhs);
    fmpz_mat_clear(zlhs);
    fmpq_mat_clear(rhs);
    fmpq_mat_clear(lhs);
    return 0;
}
