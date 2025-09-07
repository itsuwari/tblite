#ifndef TBLITE_SCF_MIXER_TYPE_H
#define TBLITE_SCF_MIXER_TYPE_H

typedef struct mixer_type mixer_type;

typedef void (*mixer_set_fn)(mixer_type *self, const double *q, int n);
typedef void (*mixer_diff_fn)(mixer_type *self, const double *q, int n);
typedef void (*mixer_get_fn)(mixer_type *self, double *q, int n);
typedef void (*mixer_next_fn)(mixer_type *self);
typedef double (*mixer_error_fn)(const mixer_type *self);

typedef void (*mixer_free_fn)(mixer_type *self);

struct mixer_type {
    mixer_set_fn set;
    mixer_diff_fn diff;
    mixer_get_fn get;
    mixer_next_fn next;
    mixer_error_fn get_error;
    mixer_free_fn del;
};

#endif
