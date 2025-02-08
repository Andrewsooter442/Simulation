#include "./helper/magnetic_field.c"
#include <stdio.h>
// #include <mathttps: //www.youtube.com/watch?v=szJn6ayPv-4h.h>

double finite_solenoid_single_winding(double s, double a, double N, double I,
                                      double L);
double ideal_long_solenoid_field(double N, double I, double L);

int main() {
  double L = 10;        // Length of solenoid (meters)
  double a = 0.05;      // Inner radius (meters)
  double b = a + 0.01;  // Outer radius (meters)
  double N = 100000000; // Total number of turns
  double I = 5.0;       // Current (amperes)
  double n = N / L;
  int num_test_points = 10;

  // Points along the axis to calculate the magnetic field
  double test_points[num_test_points + 1];
  double interval = L / num_test_points;
  for (int i = 0; i < num_test_points; i++) {
    test_points[i] = interval * i;
  }
  test_points[num_test_points] = L;

  // Calculate the ideal long solenoid magnetic field
  double ideal_B = ideal_long_solenoid_field(N, I, L);

  // Compare finite solenoid field to the derived formula and ideal field
  printf("\nComparison of Derived Finite Solenoid Formula:\n");
  printf("s (m)\tDerived B (T)\t\t\tIdeal B "
         "(T)\tError (%%)\n");

  for (int i = 0; i <= num_test_points; i++) {
    double s = test_points[i];
    double derived_magnetic_field_val =
        non_ideal_solenoid_magnetic_field_along_axis(n, I, a, b, L, 0, s);
    double finite_solenoid_single_winding_val =
        finite_solenoid_single_winding(s, a, N, I, L);
    double error_with_ideal =
        fabs((derived_magnetic_field_val - ideal_B) / ideal_B) * 100.0;

    printf("%.2f\t%.16f\t%.16f\t%.2f\n", s, derived_magnetic_field_val, ideal_B,
           error_with_ideal);
  }

  // Testing the non_ideal_solenoid_magnetic_field_along_axis with multiple
  // turns of winding with formula based for single layer of winding
  printf("\nComparison of Derived Finite Solenoid Formula:\n");
  printf("s (m)\tDerived B (T)\t\t\tFinite single winding "
         "(T)\tError (%%)\n");

  for (int i = 0; i <= num_test_points; i++) {
    double s = test_points[i];
    double derived_magnetic_field_val =
        non_ideal_solenoid_magnetic_field_along_axis(n, I, a, b, L, 0, s);
    double finite_solenoid_single_winding_val =
        finite_solenoid_single_winding(s, a, N, I, L);
    double error_with_finite =
        fabs((derived_magnetic_field_val - finite_solenoid_single_winding_val) /
             finite_solenoid_single_winding_val) *
        100.0;
    printf("%.2f\t%.20f\t\t%.20f\t\t%.2f\n", s, derived_magnetic_field_val,
           finite_solenoid_single_winding_val, error_with_finite);
  }

  // Compareing the magnetic field at any point with the formula based for
  // single layer of winding

  printf("\nComparison of Derived Finite Solenoid Formula:\n");
  printf("s (m)\tDerived B for any point (T)\t\t\tFinite single winding "
         "(T)\tError (%%)\n");

  for (int i = 0; i <= num_test_points; i++) {
    double s = test_points[i];
    double bz = 0;
    double br = 0;
    double b_theta = 0;
    compute_magnetic_field_cylindrical(a, L, I, N, 0, 0, s, &br, &b_theta, &bz);
    printf("br %.20f\nb_theta %.20f\nbz %.20f\n", br, b_theta, bz);
    double finite_solenoid_single_winding_val =
        finite_solenoid_single_winding(s, a, N, I, L);
    double error_with_finite = fabs((bz - finite_solenoid_single_winding_val) /
                                    finite_solenoid_single_winding_val) *
                               100.0;
    printf("%.2f\t%.20f\t\t%.20f\t\t%.2f\n", s, bz,
           finite_solenoid_single_winding_val, error_with_finite);
  }
  return 0;
}

// Methods to find the magnetic files in ideal cases.

// Research paper for this -
// https://www.sciencedirect.com/science/article/pii/S0304885317334662?via%3Dihub
// wikipedia article
// https://en.wikipedia.org/wiki/Solenoid#Finite_continuous_solenoid
double finite_solenoid_single_winding(double s, double a, double N, double I,
                                      double L) {
  double const_term = MU_0 * N * I * pow(a, 2) / (2 * L);
  double c =
      -((MU_0 * N * I) / (2 * L)) *
      (-s / sqrt(s * s + a * a) - (L - s) / sqrt((L - s) * (L - s) + a * a));
  return c;
}

// Simple ideal solenoid
double ideal_long_solenoid_field(double N, double I, double L) {
  double MU_0 = 4 * M_PI * 1e-7;
  double n_turns_per_length = (double)N / L;
  return MU_0 * n_turns_per_length * I;
}
