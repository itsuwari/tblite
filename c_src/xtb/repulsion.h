#ifndef TBLITE_XTB_REPULSION_H
#define TBLITE_XTB_REPULSION_H

#include "../repulsion/effective.h"
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int nat;              /* number of atoms */
    const int *id;        /* atomic type indices (0-based) */
    const double *xyz;    /* coordinates, size 3*nat */
    int nid;              /* number of distinct atomic types */
    double cutoff;        /* distance cutoff */
    double *alpha;        /* size nid*nid */
    double *zeff;         /* size nid*nid */
    double *kexp;         /* size nid*nid */
    double *rexp;         /* size nid*nid */
} tb_repulsion;

void tb_repulsion_new(tb_repulsion *r, int nat, const int *id,
                      const double *xyz, int nid);
void tb_repulsion_free(tb_repulsion *r);
void tb_repulsion_energy(const tb_repulsion *r, double *eat, double *energy);

#ifdef __cplusplus
}
#endif

#endif /* TBLITE_XTB_REPULSION_H */
