#include <stdio.h>
#include <flint/flint.h>
#include <flint/gr_types.h>
#include <flint/nf.h>

#include "interface.h"
#include "solver.h"

#define INPUT_FILE "in/A.txt"
#define OUTPUT_DIR "out/"

int main() {
    printf("Welcome to fuse! Let's bootstrap some fusion rules\n");

    gr_mat_t matrix;
    gr_ctx_t context;
    nf_t     field;

    parse_algebraic_matrix(matrix, field, context, INPUT_FILE);
    solve(matrix, field, context);

    return 0;
}
