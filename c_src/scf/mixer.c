#include "mixer.h"
#include <stdlib.h>

void new_mixer(mixer_type **self, int memory, int ndim, double damp){
    broyden_mixer *mix = malloc(sizeof(broyden_mixer));
    broyden_input in = {memory, damp};
    new_broyden(mix, ndim, in);
    *self = (mixer_type*)mix;
}
