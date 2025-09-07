#include "gfn.h"
#include <math.h>
#include <stdlib.h>

/* constants from Fortran implementation */
static const double ka = 10.0; /* steepness of first function */
static const double kb = 20.0; /* steepness of second function */
static const double r_shift = 2.0; /* offset for second function */

void new_gfn_ncoord(gfn_ncoord_type *self, size_t nid, const double *rcov,
                    double cutoff) {
    self->nid = nid;
    self->cutoff = cutoff > 0.0 ? cutoff : 25.0; /* default cutoff */
    self->directed_factor = 1.0; /* not used */
    self->rcov = (double *)malloc(nid * sizeof(double));
    if (rcov) {
        for (size_t i = 0; i < nid; ++i) self->rcov[i] = rcov[i];
    } else {
        for (size_t i = 0; i < nid; ++i) self->rcov[i] = 0.0; /* caller may fill later */
    }
}

static inline double gfn_count(double k, double r, double r0) {
    return 1.0 / (1.0 + exp(-k * (r0 / r - 1.0)));
}

static inline double gfn_dcount(double k, double r, double r0) {
    double expterm = exp(-k * (r0 / r - 1.0));
    return (-k * r0 * expterm) / (r * r * pow(expterm + 1.0, 2));
}

double ncoord_count(const gfn_ncoord_type *self, size_t izp, size_t jzp,
                    double r) {
    double rc = self->rcov[izp] + self->rcov[jzp];
    return gfn_count(ka, r, rc) * gfn_count(kb, r, rc + r_shift);
}

double ncoord_dcount(const gfn_ncoord_type *self, size_t izp, size_t jzp,
                      double r) {
    double rc = self->rcov[izp] + self->rcov[jzp];
    double c1 = gfn_count(ka, r, rc);
    double c2 = gfn_count(kb, r, rc + r_shift);
    double d1 = gfn_dcount(ka, r, rc);
    double d2 = gfn_dcount(kb, r, rc + r_shift);
    return d1 * c2 + c1 * d2;
}
