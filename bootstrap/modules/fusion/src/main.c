#include <stdio.h>
#include <flint/gr.h>

#include "interface.h"

#define INPUT_FILE "in/A.txt"
#define OUTPUT_DIR "out/"

int main() {
    printf("Welcome to fuse! Let's bootstrap some fusion rules\n");

    gr_mat_t matrix;
    gr_ctx_t context;

    parse_algebraic_matrix(matrix, context, INPUT_FILE);

    return 0;
}
