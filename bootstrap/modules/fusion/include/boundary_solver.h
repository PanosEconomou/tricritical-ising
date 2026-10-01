#ifndef BOUNDARY_SOLVER_H
#define BOUNDARY_SOLVER_H

#include <flint/gr_mat.h>
#include <flint/nf.h>
#include <flint/qqbar.h>

#ifdef __cplusplus
extern "C" {
#endif

int solve_boundary_bootstrap(const slong* indices, slong nindices,
                             gr_mat_struct* matrix, nf_struct* field,
                             qqbar_struct* field_generator, gr_ctx_struct* context);

#ifdef __cplusplus
}
#endif

#endif
