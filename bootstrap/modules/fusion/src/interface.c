#include "interface.h"

#include <flint/acb.h>
#include <flint/flint.h>
#include <flint/fmpq_types.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/gr_types.h>
#include <flint/qqbar.h>
#include <flint/gr.h>
#include <stdio.h>

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
    slong       ndim;
    slong*      dims;
    slong       total;
    slong       nvals;
    slong       precision;
    qqbar_ptr   values;
    slong*      idx;

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
        degree = 0;

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
static int build_number_field(gr_ctx_t context, fmpq_poly_t** reps, 
                              const input_data_t* in) 
{
    int status = -1;

    qqbar_t     a;          qqbar_init(a);
    qqbar_t     candidate;  qqbar_init(candidate);
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
        }
        
    }

cleanup:
    fmpq_poly_clear(scratch);
    qqbar_clear(candidate);
    qqbar_clear(a);
    return status;
}

int parse_algebraic_matrix(char *filename) 
{

    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "ERROR: %s can't be opened", filename);
        return -1;
    }

    input_data_t in = {0};

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

    fclose(file);
    free_input_data(&in);
    return 0;

error:
    fclose(file);
    free_input_data(&in);

    return -1;
}
