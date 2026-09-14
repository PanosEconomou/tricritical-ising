#include "solver.h"

#include <flint/gr.h>
#include <flint/gr_mat.h>
#include <flint/fmpq_vec.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>
#include <stdio.h>


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

int solve(gr_mat_struct* matrix, nf_t field, gr_ctx_struct *context)
{

    slong rows   =   gr_mat_nrows(matrix, context); 
    slong cols   =   gr_mat_ncols(matrix, context); 
    slong degree =   fmpq_poly_degree(field->pol);
    fmpq_mat_t lhs;  fmpq_mat_init(lhs, rows * degree, cols);
    fmpq_mat_t rhs;  fmpq_mat_init(rhs, rows * degree, 1);
    fmpz_mat_t zlhs; fmpz_mat_init(zlhs, rows * degree, cols);
    fmpz_mat_t zrhs; fmpz_mat_init(zrhs, rows * degree, 1);

    build_lhs(lhs, matrix, field, context);
    build_rhs(rhs, matrix, 1, 5, field, context);
    zfy(zrhs, zlhs, rhs, lhs);

    // fmpq_mat_print(lhs); printf("\n");
    // fmpz_mat_print(zrhs); printf("\n");
    
    fmpz_mat_clear(zrhs);
    fmpz_mat_clear(zlhs);
    fmpq_mat_clear(rhs);
    fmpq_mat_clear(lhs);
    return 0;
}
