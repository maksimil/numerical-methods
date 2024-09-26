#include "config.h"
#include "task.h"

int main(int argc, char *argv[]) {
  LOG_INFO("\e[1;32m>> Part 1.1\e[0m\n");

  const Scalar f_a = task_f(task_a);
  const Scalar f_b = task_f(task_b);
  const Scalar f_m = task_f((task_a + task_b) / 2.);
  const Scalar h = task_b - task_a;

  const Scalar left_quad = f_a * h;
  const Scalar right_quad = f_b * h;
  const Scalar middle_quad = f_m * h;
  const Scalar trapezoid = (f_a + f_b) / 2. * h;
  const Scalar simpson = (f_a + 4. * f_m + f_b) / 6. * h;

  LOG_INFO("%10s  %20s  %11s  %11s\n", "Method", "Answer", "Error",
           "Rel Error");

  LOG_INFO("%10s  %20.16f  %11.4e  %11.4e\n", "scipy quad", task_ans, 0., 0.);
  LOG_INFO("%10s  %20.16f  %11.4e  %11.4e\n", "LQuad", left_quad,
           task_ans - left_quad, (task_ans - left_quad) / task_ans);
  LOG_INFO("%10s  %20.16f  %11.4e  %11.4e\n", "RQuad", right_quad,
           task_ans - right_quad, (task_ans - right_quad) / task_ans);
  LOG_INFO("%10s  %20.16f  %11.4e  %11.4e\n", "MQuad", middle_quad,
           task_ans - middle_quad, (task_ans - middle_quad) / task_ans);
  LOG_INFO("%10s  %20.16f  %11.4e  %11.4e\n", "Trapezoid", trapezoid,
           task_ans - trapezoid, (task_ans - trapezoid) / task_ans);
  LOG_INFO("%10s  %20.16f  %11.4e  %11.4e\n", "Simpson", simpson,
           task_ans - simpson, (task_ans - simpson) / task_ans);
}
