#ifndef FUSE_SOLVER
#define FUSE_SOLVER

#include <flint/nf.h>
#include <flint/gr_mat.h>
#include <flint/qqbar.h>

int solve(gr_mat_t matrix, nf_t field, qqbar_t field_generator, gr_ctx_t context);

#endif
