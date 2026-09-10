#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/qqbar.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>
#include <stdio.h>

#include "interface.h"

#define FIELD_MAX_BITS  1000
#define FIELD_PRECISION 256
#define FIELD_MAX_N     64

static int read_fmpz(FILE* file, fmpz_t out)
{
    char buffer[8192];
    if (fscanf(file, "%8191s", buffer) != 1) return 0;
    return fmpz_set_str(out, buffer, 10) == 0;
}

static int read_slong(FILE* file, slong* out)
{
    return fscanf(file, "%ld", out) == 1;
}

typedef struct {
    slong      ndim;
    slong*     dims;
    slong      total;
    slong      nvals;
    slong      precision;
    qqbar_ptr  values;
    slong*     idx;

} input_data_t;

void free_input_data(input_data_t* in) 
{
    if (in == NULL) return;

    if (in->dims)    { flint_free(in->dims);                    in->dims   = NULL; }
    if (in->idx)     { flint_free(in->idx);                     in->idx    = NULL; }
    if (in->values)  { _qqbar_vec_clear(in->values, in->nvals); in->values = NULL; }

    in->ndim = in->total = in->nvals = in->precision = 0;
}

static int read_input_metadata(FILE* file, input_data_t* in) 
{
    if (!read_slong(file, &in->ndim)) { 
        fprintf(stderr, "ERROR: Input file is not formatted correctly\n"); 
        return -1;
    }

    if (in->ndim <= 0) goto error;

    in->dims  = flint_malloc((in->ndim) * sizeof(slong));
    in->total = 1;
    for (slong k = 0; k < in->ndim; k++) { 
        if (!read_slong(file, &in->dims[k])) goto error;

        in->total *= in->dims[k]; 
    }

    if (!read_slong(file, &in->nvals))      goto error;
    if (!read_slong(file, &in->precision))  goto error;

    return 0;

error:
    fprintf(stderr, "ERROR: Input file is not formatted correctly\n"); 
    return -1;
}

static int read_input_polynomials(FILE* file, input_data_t* in) 
{
    int   status = -1;
    slong wp     = (slong) (in->precision * 3.33) + 64; // in base 2 
    
    in->values = _qqbar_vec_init(in->nvals);

    // Scratch
    fmpz_poly_t polynomial;  fmpz_poly_init(polynomial);
    fmpz_t      scale;       fmpz_init(scale);
    fmpz_t      coefficient; fmpz_init(coefficient);
    fmpz_t      target;      fmpz_init(target);
    arb_t       ball;        arb_init(ball);
    acb_t       z;           acb_init(z);
    qqbar_ptr   roots  = NULL;
    slong       nroots = 0;

    fmpz_set_ui(scale, 10);
    fmpz_pow_ui(scale, scale, in->precision);

    for (slong i = 0; i < in->nvals; i++) {

        slong degree;
        if (!read_slong(file, &degree) || degree < 1) goto cleanup;

        // Read the polynomial
        fmpz_poly_zero(polynomial);
        for (slong j = 0; j <= degree; j++) {
            if (!read_fmpz(file, coefficient)) goto cleanup;
            fmpz_poly_set_coeff_fmpz(polynomial, j, coefficient);
        }
        
        // store the numerical root for the target
        if (!read_fmpz(file, target)) goto cleanup;
        arb_set_fmpz(ball, target);
        arb_add_error_2exp_si(ball, -1);        // set precision to +- 1/2
        arb_div_fmpz(ball, ball, scale, wp);

        // Extract the polynomial's roots
        nroots = degree;
        roots  = _qqbar_vec_init(nroots);
        qqbar_roots_fmpz_poly(roots, polynomial, QQBAR_ROOTS_IRREDUCIBLE);

        slong match = -1, hits = 0;
        for (slong j = 0; j < nroots; j++) {
            qqbar_get_acb(z, roots + j, wp);
            if (arb_contains_zero(acb_imagref(z)) && 
                arb_overlaps(acb_realref(z), ball)) {
                match = j;
                hits++;
            }
        }

        if (hits != 1) {
            flint_fprintf(stderr, "ERROR: value %wd matched %wd roots of ", i, hits);
            fmpz_poly_print_pretty(polynomial, "x");
            flint_fprintf(stderr, "\n");
            status = -1;
            goto cleanup;
        }

        qqbar_set(in->values + i, roots + match);

        _qqbar_vec_clear(roots, nroots);
        roots  = NULL;

        // qqbar_print(in->values + i); flint_printf("\n");
    }

    status = 0;
    

cleanup:
    if (roots) _qqbar_vec_clear(roots, nroots);
    acb_clear(z);
    arb_clear(ball);
    fmpz_clear(target);
    fmpz_clear(coefficient);
    fmpz_clear(scale);
    fmpz_poly_clear(polynomial);
    return status;
} 

static int read_input_entries(FILE* file, input_data_t* in) 
{

    if (!in->total) goto error; 
    in->idx = flint_malloc(in->total * sizeof(slong));

    for (slong entry = 0; entry < in->total; entry++) {
        if (!read_slong(file, &in->idx[entry])) goto error;
        in->idx[entry]--; /* mathematica typically exports to 1 based arrays */
        if (in->idx[entry] < 0 || in->idx[entry] >= in->nvals) goto error;
    }

    return 0;

error:
    fprintf(stderr, "ERROR: Input file is not formatted correctly\n"); 
    return -1;
}


/* 
 * There is something called the primitive element theorem 
 * [https://en.wikipedia.org/wiki/Primitive_element_theorem], and the idea is that these
 * finitely generated modules can actually be generated by one element a, 
 * so this function looks for it.
 */
static int build_number_field(nf_t field, gr_ctx_t context, fmpq_poly_struct** reps,
                              const input_data_t* in) 
{
    int status = -1;

    qqbar_t     a;          qqbar_init(a);
    qqbar_t     candidate;  qqbar_init(candidate);
    fmpq_t      rational;   fmpq_init(rational);
    fmpq_poly_t scratch;    fmpq_poly_init(scratch);

    // start with a = 1
    qqbar_one(a);

    for (slong i = 0; i < in->nvals; i++) {

        // Skip value if it is already in the ring made by a
        if (qqbar_express_in_field(scratch, a, in->values + i, 
                                   FIELD_MAX_BITS, 0, FIELD_PRECISION)) continue;

        slong n;
        for (n=1; n <= FIELD_MAX_N; n++) {
            qqbar_mul_si(candidate, in->values + i, n);
            qqbar_add(candidate, candidate, a);

            // Check if this field works.
            if (qqbar_express_in_field(scratch, candidate, a, 
                                       FIELD_MAX_BITS, 0, FIELD_PRECISION) && 
                qqbar_express_in_field(scratch, candidate, in->values + i,
                                       FIELD_MAX_BITS, 0, FIELD_PRECISION)) {
               break; 
            }
        }

        if (n > FIELD_MAX_N) {
            fprintf(stderr, "ERROR: can't generate algebraic field for entry.\n");
            goto cleanup;
        }

        qqbar_set(a, candidate);
        // qqbar_print(a); printf("\n");
    }

    // Build the actual algebraic ring
    if (qqbar_degree(a) == 1) {
        fprintf(stderr, "ERROR: all values are rational; no number field to build.\n");
        goto cleanup;
    } else {
        fmpq_poly_set_fmpz_poly(scratch, QQBAR_POLY(a));
        gr_ctx_init_nf(context, scratch);
        nf_init(field, scratch);
    }

    // Convert the entries to the new algebraic ring
    *reps = flint_malloc(in->nvals * sizeof(fmpq_poly_struct));
    for (slong i = 0; i < in->nvals; i++) fmpq_poly_init(*reps + i);

    for (slong i = 0; i < in->nvals; i++) {
        if (qqbar_is_rational(in->values + i)) {
            qqbar_get_fmpq(rational, in->values + i);
            fmpq_poly_set_fmpq(*reps + i, rational);
            continue;
        } 

         if (!qqbar_express_in_field(*reps + i, a, in->values + i,
                                    FIELD_MAX_BITS, 0, FIELD_PRECISION)) {
            flint_fprintf(stderr,
                "ERROR: value %wd not expressible in the final field "
                "(raise FIELD_MAX_BITS or FIELD_PREC)\n", i);
            flint_free(*reps);
            *reps = NULL;
            nf_clear(field);
            gr_ctx_clear(context);
            goto cleanup;
        }

        // fmpq_poly_print_pretty(*reps+i, "x"); printf("\n");
    }

    status = 0;

cleanup:
    fmpq_clear(rational);
    fmpq_poly_clear(scratch);
    qqbar_clear(candidate);
    qqbar_clear(a);
    return status;
}

static int build_matrix(gr_mat_t matrix, gr_ctx_t context, const nf_t field, 
                        const fmpq_poly_struct* reps, const input_data_t* in)
{
    gr_mat_init(matrix, in->dims[0], in->dims[1], context);

    for (slong i = 0; i < in->total; i++) {
        slong row = i / in->dims[1];
        slong col = i % in->dims[1];
        gr_ptr entry = gr_mat_entry_ptr(matrix, row, col, context);

        nf_elem_set_fmpq_poly(entry, reps + in->idx[i], field);
    }

    return 0;
}

int parse_algebraic_matrix(gr_mat_t matrix, gr_ctx_t context, char* filename) 
{
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "ERROR: %s can't be opened", filename);
        return -1;
    }

    input_data_t in = {0};

    fmpq_poly_struct* reps = NULL;
    nf_t              field;
    int               have_field = 0;

    if (read_input_metadata(file, &in)) {
        fprintf(stderr, "ERROR: can't parse metadata from %s", filename);
        goto error;
    }

    if (in.ndim != 2) {
        flint_fprintf(stderr, 
                      "ERROR: %s does not contain a matrix, instead it contains a %wd-tensor",
                      filename, in.ndim
                      );
        goto error;
    }

    if (read_input_polynomials(file, &in)) {
        fprintf(stderr, "ERROR: can't parse polynomials from %s", filename);
        goto error;
    }

    if (read_input_entries(file, &in)) {
        fprintf(stderr, "ERROR: can't parse entries from %s", filename);
        goto error;
    }

    if (build_number_field(field, context, &reps, &in)) {
        fprintf(stderr, "ERROR: can't build the number field for %s", filename);
        goto error;
    }
    have_field = 1;

    if (build_matrix(matrix, context, field, reps, &in)) {
        fprintf(stderr, "ERROR: can't build the matrix for %s", filename);
        gr_mat_clear(matrix, context);
        gr_ctx_clear(context);
        goto error;
    }

    for (slong i = 0; i < in.nvals; i++) fmpq_poly_clear(reps + i);
    flint_free(reps);
    nf_clear(field);

    fclose(file);
    free_input_data(&in);
    return 0;

error:
    if (reps) {
        for (slong i = 0; i < in.nvals; i++) fmpq_poly_clear(reps + i);
        flint_free(reps);
    }
    if (have_field) nf_clear(field);
    fclose(file);
    free_input_data(&in);

    return -1;
}
