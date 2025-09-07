#ifndef TBLITE_NCOORD_GFN_H
#define TBLITE_NCOORD_GFN_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t nid;        /* number of element types */
    double *rcov;      /* covalent radii per element */
    double cutoff;     /* real-space cutoff */
    double directed_factor; /* not used but kept for parity */
} gfn_ncoord_type;

/* initialise container with optional rcov array (length nid) */
void new_gfn_ncoord(gfn_ncoord_type *self, size_t nid, const double *rcov,
                    double cutoff);

/* double exponential counting function */
double ncoord_count(const gfn_ncoord_type *self, size_t izp, size_t jzp,
                    double r);

/* derivative of counting function */
double ncoord_dcount(const gfn_ncoord_type *self, size_t izp, size_t jzp,
                      double r);

#ifdef __cplusplus
}
#endif

#endif /* TBLITE_NCOORD_GFN_H */
