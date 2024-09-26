#ifndef __TASK_H__
#define __TASK_H__

#include "config.h"
#include "math.h"

static const Scalar task_a = 1.3;
static const Scalar task_b = 2.2;
static const Scalar task_alpha = 0;
static const Scalar task_beta = 5. / 6;
static const Scalar task_ans = 3.0765665777241504;
static const Scalar task_anp = 1.5923668372743820;

Scalar task_f(Scalar x) {
  return 4. * cos(0.5 * x) * exp(-5. * x / 4.) +
         2. * sin(4.5 * x) * exp(x / 8.) + 2.;
}

#endif // __TASK_H__
