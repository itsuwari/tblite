#ifndef TBLITE_COULOMB_CHARGE_H
#define TBLITE_COULOMB_CHARGE_H

/* Simple isotropic Coulomb charge model */

void effective_charge_energy(int nat, const double *q, const double *gamma,
                             double *energy);
void effective_charge_potential(int nat, const double *q, const double *gamma,
                                double *vat);

#endif /* TBLITE_COULOMB_CHARGE_H */
