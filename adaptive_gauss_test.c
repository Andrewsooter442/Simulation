#include "./helper/adaptive_gauss_kronrod.c"
#include <math.h>
#include <stdio.h>

extern long double adaptive_gauss_kronrod(long double (*func)(long double),
                                          long double a, long double b,
                                          long double tol);

// Define test functions
long double trigonometric(long double x) { return sin(x) * cos(x); }

long double exponential(long double x) { return exp(-x * x); }

long double logarithmic(long double x) { return log(1 + x); }

long double singularity(long double x) { return 1.0 / sqrt(x); }

long double oscillatory(long double x) { return sinl(50 * x) * expl(-x); }

long double log_power(long double x) { return powl(x, 3) * logl(x); }

long double shifted_gaussian(long double x) { return expl(-powl(x - 2, 2)); }

long double high_degree_polynomial(long double x) {
  return powl(x, 10) - 5 * powl(x, 7) + 2 * powl(x, 4) - powl(x, 2) + 3;
}

long double complicated_function(long double x) {
  long double numerator = powl(x, 24) + 3 * powl(x, 10) - 5 * powl(x, 8) +
                          2 * powl(x, 6) - powl(x, 4) + 7 * powl(x, 2);
  long double denominator = powl(x, 2) * cosl(powl(x, 3));

  return numerator / denominator;
}

// Wrapper function to test adaptive_gauss_kronrod with a given function
void test_integration(long double (*func)(long double), long double a,
                      long double b, long double expected,
                      const char *test_name) {
  long double tol = 1e-9;
  long double result = adaptive_gauss_kronrod(func, a, b, tol);
  printf("Test: %s\n", test_name);
  printf("Computed Integral: %.10Lf\n", result);
  printf("Expected Integral: %.10Lf\n", expected);
  printf("Error: %.10Lf\n\n\n", fabsl(result - expected));
}

int main() {
  long double trig_expected = 0.25;
  test_integration(trigonometric, 0.0, M_PI / 2, trig_expected,
                   "Trigonometric Function");

  long double exp_expected = 1.493648265624854;
  test_integration(exponential, -1.0, 1.0, exp_expected,
                   "Exponential Function");

  long double log_expected = (2 * log(2) - 1) - (0);
  test_integration(logarithmic, 0.0, 1.0, log_expected, "Logarithmic Function");

  long double singular_expected = 2.0;
  test_integration(singularity, 0.000001, 1.0, singular_expected,
                   "Singular Function (1/sqrt(x))");

  test_integration(high_degree_polynomial, 0.0, 1.0,
                   3 + 1.0 / 11 - 5.0 / 8 + 2.0 / 5 - 1.0 / 3,
                   "High Degree Polynomial");
  test_integration(oscillatory, 0.0, 10.0,
                   (expl(-10) * (sinl(500) - 50 * cosl(500)) + 50) / 2501,
                   "Oscillatory Function");
  test_integration(log_power, 1.0, 3.0,
                   ((pow(3, 4) / 4) * log(3) - pow(3, 4) / 16) -
                       ((pow(1, 4) / 4) * log(1) - pow(1, 4) / 16),
                   "Log-Power Function");
  test_integration(shifted_gaussian, -10.0, 10.0, sqrtl(M_PI),
                   "Shifted Gaussian");

  test_integration(complicated_function, 0.1, 2.0, 0,
                   "Complicated Rational Function");

  return 0;
}
