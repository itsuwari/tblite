#ifndef TBLITE_SCF_MIXER_H
#define TBLITE_SCF_MIXER_H

#include "mixer_type.h"
#include "mixer_broyden.h"

typedef struct {
    broyden_input *broyden;
} mixer_input;

void new_mixer(mixer_type **self, int memory, int ndim, double damp);

#endif
