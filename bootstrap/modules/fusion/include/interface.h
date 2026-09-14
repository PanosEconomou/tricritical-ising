#ifndef FUSE_INTERFACE
#define FUSE_INTERFACE

#include <flint/gr.h>
#include <flint/nf.h>
#include <flint/qqbar.h>

int parse_algebraic_matrix(gr_mat_t matrix, nf_t field, qqbar_t field_generator, 
                           gr_ctx_t context, char* filename);

#endif
