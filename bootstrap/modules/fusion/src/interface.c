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

#define FIELD_START_PREC     64                   /* first LLL precision tried (bits)   */
#define FIELD_QUICK_PREC     4096                 /* cap for the cheap membership probe */
#define FIELD_MAX_PREC       (((slong) 1) << 17)  /* hard cap: 131072 bits              */
#define FIELD_MAX_N          64
#define FIELD_DEGREE_PROBES  3

int load_ishibashi_indices(slong** indices, char* filename)
{
    FILE* file = fopen(filename, "r");
    slong n = -1;
    slong i = 0;

    if (!file) return -1;
    if (fscanf(file, "%ld", &n) == 1 && n > 0) {
        *indices = flint_malloc(n * sizeof(slong)); 
        while(i < n && fscanf(file, "%ld", *indices + i) == 1) i++;
    }

    fclose(file);

    if (n <= 0 || i != n) {
        flint_free(*indices);
        *indices = NULL;
        return -1;
    }

    return (int)n;
}

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
 * qqbar_express_in_field() does ONE LLL attempt at a fixed precision/height and
 * never raises them itself. A relation of degree d with B-bit coefficients needs
 * roughly (d+1)*B bits of working precision, so we double until it is found.
 * Returning 1 is a proof: FLINT verifies f(alpha) == x exactly.
 * Returning 0 only means "not found up to max_prec".
 */
static int express_adaptive(fmpq_poly_t res, const qqbar_t alpha,
                            const qqbar_t x, slong max_prec)
{
    slong d = qqbar_degree(alpha);

    if (qqbar_is_rational(x)) {
        fmpq_t q;
        fmpq_init(q);
        qqbar_get_fmpq(q, x);
        fmpq_poly_set_fmpq(res, q);
        fmpq_clear(q);
        return 1;
    }

    if (d % qqbar_degree(x) != 0) return 0;   /* x cannot lie in Q(alpha) */

    for (slong prec = FIELD_START_PREC; prec <= max_prec; prec *= 2) {
        /* LLL on d+1 numbers at `prec` bits can only certify ~prec/(d+1)-bit
           coefficients; a larger bound just yields spurious candidates that
           then fail the exact check. */
        slong max_bits = FLINT_MAX(8, prec / (d + 1));
        if (qqbar_express_in_field(res, alpha, x, max_bits, 0, prec))
            return 1;
    }
    return 0;
}

/* 
 * There is something called the primitive element theorem 
 * [https://en.wikipedia.org/wiki/Primitive_element_theorem], and the idea is that these
 * finitely generated modules can actually be generated by one element a, 
 * so this function looks for it.
 */
static int build_number_field(nf_t field, qqbar_t generator, gr_ctx_t context,
                              fmpq_poly_struct** reps, const input_data_t* in)
{
    int status = -1;

    qqbar_t     a;          qqbar_init(a);
    qqbar_t     candidate;  qqbar_init(candidate);
    fmpq_poly_t scratch;    fmpq_poly_init(scratch);
    qqbar_ptr   probes = _qqbar_vec_init(FIELD_DEGREE_PROBES);

    qqbar_zero(a);   /* Q(0) = Q */

    for (slong i = 0; i < in->nvals; i++) {
        const qqbar_struct* v = in->values + i;
        slong da = qqbar_degree(a);
        slong dv = qqbar_degree(v);

        if (dv == 1) continue;   /* rational */

        /* 1. Fast path: most entries are already in Q(a). */
        if (express_adaptive(scratch, a, v, FIELD_QUICK_PREC)) continue;

        /* 2. Exact degree test. If v is in Q(a), every a + n v is too, so
              deg(a + n v) <= deg(a). If v is not in Q(a), all but finitely
              many n give a primitive element of Q(a, v), of degree > deg(a). */
        slong best = da;
        for (slong n = 1; n <= FIELD_DEGREE_PROBES; n++) {
            qqbar_mul_si(probes + n - 1, v, n);
            qqbar_add(probes + n - 1, probes + n - 1, a);
            best = FLINT_MAX(best, qqbar_degree(probes + n - 1));
        }

        if (best == da) {
            /* No evidence v is outside Q(a): it just needs more precision. */
            if (express_adaptive(scratch, a, v, FIELD_MAX_PREC)) continue;
            flint_fprintf(stderr,
                "ERROR: could not express value %wd in the current field of "
                "degree %wd (raise FIELD_MAX_PREC)\n", i, da);
            goto cleanup;
        }

        /* 3. v is provably outside Q(a). Find n with Q(a + n v) = Q(a, v).
              It suffices to show a is in Q(a + n v): then
              v = ((a + n v) - a) / n is too. */
        int found = 0;
        for (slong n = 1; n <= FIELD_MAX_N && !found; n++) {
            if (n <= FIELD_DEGREE_PROBES) {
                qqbar_set(candidate, probes + n - 1);
            } else {
                qqbar_mul_si(candidate, v, n);
                qqbar_add(candidate, candidate, a);
            }

            slong dc = qqbar_degree(candidate);
            if (dc < best) continue;                 /* in a proper subfield */
            best = dc;
            if (dc % da != 0 || dc % dv != 0) continue;

            if (dc == da * dv)                       /* [Q(a,v):Q] <= da*dv */
                found = 1;
            else
                found = express_adaptive(scratch, candidate, a, FIELD_MAX_PREC);
        }

        if (!found) {
            flint_fprintf(stderr,
                "ERROR: can't generate algebraic field for entry %wd, root of ", i);
            fmpz_poly_fprint_pretty(stderr, QQBAR_POLY(v), "x");
            flint_fprintf(stderr, "\n");
            goto cleanup;
        }

        qqbar_set(a, candidate);
    }

    if (qqbar_degree(a) == 1) {
        fprintf(stderr, "ERROR: all values are rational; no number field to build.\n");
        goto cleanup;
    }

    /* Build the number field */
    fmpq_poly_set_fmpz_poly(scratch, QQBAR_POLY(a));
    gr_ctx_init_nf(context, scratch);
    nf_init(field, scratch);

    /* Convert every value to a polynomial in a */
    *reps = flint_malloc(in->nvals * sizeof(fmpq_poly_struct));
    for (slong i = 0; i < in->nvals; i++) fmpq_poly_init(*reps + i);

    for (slong i = 0; i < in->nvals; i++) {
        if (!express_adaptive(*reps + i, a, in->values + i, FIELD_MAX_PREC)) {
            flint_fprintf(stderr,
                "ERROR: value %wd not expressible in the final field "
                "(raise FIELD_MAX_PREC)\n", i);
            for (slong j = 0; j < in->nvals; j++) fmpq_poly_clear(*reps + j);
            flint_free(*reps);
            *reps = NULL;
            nf_clear(field);
            gr_ctx_clear(context);
            goto cleanup;
        }
    }

    qqbar_set(generator, a);
    status = 0;

cleanup:
    _qqbar_vec_clear(probes, FIELD_DEGREE_PROBES);
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

int parse_algebraic_matrix(gr_mat_t matrix, nf_t field, qqbar_t field_generator, gr_ctx_t context, char* filename) 
{
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "ERROR: %s can't be opened", filename);
        return -1;
    }

    input_data_t in = {0};

    fmpq_poly_struct* reps = NULL;
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

    if (build_number_field(field, field_generator, context, &reps, &in)) {
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
