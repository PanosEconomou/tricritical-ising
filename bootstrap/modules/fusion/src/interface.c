#include "interface.h"

#include <pari/pari.h>
#include <flint/gr.h>

int pari_available = 0;

static void close_pari() {
    if (pari_available) {
        pari_close();
    }
}

static void ensure_pari() {
    if (!pari_available) {
        pari_init(256 * 1024 * 1024, 500000);
        atexit(close_pari);
        pari_available = 1;
    }
}

int parse_algebraic_matrix(gr_mat_struct *matrix, gr_ctx_struct *context, 
                           char *filename) 
{
    ensure_pari();

    GEN file = gp_readvec_file(filename);
    output(file);

    return 0;
}
