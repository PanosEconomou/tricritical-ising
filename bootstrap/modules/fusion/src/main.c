#include <flint/flint.h>
#include <stdio.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

#include "interface.h"

#define INPUT_FILE "in/A.txt"
#define OUTPUT_DIR "out/"

int main() {
    printf("Welcome to fuse! Let's bootstrap some fusion rules\n");

    parse_algebraic_matrix(INPUT_FILE);

    return 0;
}
