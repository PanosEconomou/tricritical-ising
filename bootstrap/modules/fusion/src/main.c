#include <stdio.h>
#include <flint/flint.h>
#include <flint/gr_types.h>
#include <flint/nf.h>

#include "boundary_solver.h"
#include "interface.h"
#include "solver.h"

#define S_MATRIX_INPUT_FILE "in/S_B_DE.txt"
#define INDEX_INPUT_FILE "in/I_B_DE.txt"
#define OUTPUT_DIR "out/"

int main() {
    printf("Welcome to fuse! Let's bootstrap some fusion rules\n");

    slong* indices; 
    gr_mat_t matrix;
    gr_ctx_t context;
    nf_t     field;
    qqbar_t  field_generator; qqbar_init(field_generator);

    slong n = load_ishibashi_indices(&indices, INDEX_INPUT_FILE);
    parse_algebraic_matrix(matrix, field, field_generator, context, S_MATRIX_INPUT_FILE);
    // solve(matrix, field, field_generator, context);
    solve_boundary_bootstrap(indices, n, matrix, field, field_generator, context);


    return 0;
}
