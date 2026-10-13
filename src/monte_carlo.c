#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "sim.h"
#include "ran.h"
extern int nnsize;


int monte_carlo_sweep(Par *par, double *pos) {
  int i, j, d, naccept = 0;
  double r2_new, r2_old, ediff;
  double newpos[D], dist[D];

  for (i = 0; i < par->n; i++) {
    for (d = 0; d < D; d++) {
      double delta = par->b * dran_sign(); // random step in [-b, b)
      newpos[d] = inrange(pos[D * i + d] + delta, par->L[d]);
    }

    // Energy difference, only pairs involving particle i contribute
    ediff = 0.0;
    for (j = 0; j < par->n; j++) {
      if (j == i)
	continue;
      r2_new = dist2(par->L, newpos, &pos[D * j], dist);
      r2_old = dist2(par->L, &pos[D * i], &pos[D * j], dist);
      ediff += pair_interaction(r2_new) - pair_interaction(r2_old);
    }

    // Accept with probability min(1, exp(-ediff / T))
    if (ediff < 0.0 || dran() < exp(-ediff / par->T)) {
      for (d = 0; d < D; d++)
	      pos[D * i + d] = newpos[d];
        naccept++;
    }
  }

#ifdef FAST
  check_neighbor_list(par, pos);
#endif

  return naccept;
}
