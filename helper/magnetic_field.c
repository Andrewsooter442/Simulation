#include <math.h>
#include <stdio.h>
double MU_0 = (4 * M_PI * 1e-7); // Permeability of free space

// Reference https://ieeexplore.ieee.org/document/8008997
double derived_magnetic_field_s(double s, double a, double b, double N,
                                double L, double I) {
  // Calculate current density J (current per unit cross-sectional area)
  double numerator = (double)N * I;
  double denumerator = (M_PI * (b * b - a * a));
  double J = numerator / denumerator;
  // printf("J %f\t %f \n", numerator, denumerator);

  // Define the bounds of the solenoid (centered at z = 0)
  double z1 = -L / 2; // Start of the solenoid
  double z2 = L / 2;  // End of the solenoid

  // Compute terms for the derived formula
  double term1 = log((b + sqrt(pow(b, 2) + pow((s - L), 2))) /
                     (a + sqrt(pow(a, 2) + pow((s - L), 2))));

  double term2 = log((b + sqrt(pow(b, 2) + pow(s, 2))) /
                     (a + sqrt(pow(a, 2) + pow(s, 2))));
  // printf("test %f\t\n", term1);

  // Magnetic field calculation using the derived formula
  double B = MU_0 * (J / 2) * (term2 + term1);
  return B;
}

// n  Number of turns per unit length (turns per meter)
// I  Current in amperes (A)
// R1 Inner radius of solenoid (meters)
// R2 Outer radius of solenoid (meters)
// b  Lower bound of solenoid along z-axis (meters) = 0 if center is not at
// origin but at L/2 a  Upper bound of solenoid along z-axis (meters) z  Point
// along the z-axis where Bz is calculated (meters)

// Reference https://ieeexplore.ieee.org/document/8008997
double non_ideal_solenoid_magnetic_field_along_axis(double n, double I,
                                                    double R1, double R2,
                                                    double a, double b,
                                                    double z) {
  if (R2 <= R1) {
    printf("Error: Outer radius R2 must be greater than inner radius R1.\n");
    return -1;
  }

  // Compute logarithmic terms
  double term1 = (z - a) * log((sqrt(R2 * R2 + (z - a) * (z - a)) + R2) /
                               (sqrt(R1 * R1 + (z - a) * (z - a)) + R1));
  double term2 = (z + b) * log((sqrt(R2 * R2 + (z + b) * (z + b)) + R2) /
                               (sqrt(R1 * R1 + (z + b) * (z + b)) + R1));

  // Compute the magnetic field
  double Bz = (0.5 * MU_0 * n * I / (R2 - R1)) * (term2 - term1);

  return Bz;
}

// https://ieeexplore.ieee.org/abstract/document/10057086
//
// Define constants
int GAUSS_POINTS = 10; // Number of quadrature points
void compute_magnetic_field_cylindrical(double R, double L, double I, int N,
                                        double r_obs, double phi_obs,
                                        double z_obs, double *Br, double *Bphi,
                                        double *Bz) {
  *Br = 0.0;
  *Bphi = 0.0;
  *Bz = 0.0;

  double gauss_legendre_x[] = {
      -0.9739065285, -0.8650633666, -0.6794095682, -0.4333953941, -0.1488743389,
      0.1488743389,  0.4333953941,  0.6794095682,  0.8650633666,  0.9739065285};
  double gauss_legendre_w[] = {
      0.0666713443, 0.1494513491, 0.2190863625, 0.2692667193, 0.2955242247,
      0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443};

  double dz = L / N;    // Step size along solenoid axis
  double I_dz = I * dz; // Current element contribution

  // Loop over solenoid turns (along z-axis)
  for (int k = 0; k < N; k++) {
    double zk = -L / 2.0 + k * dz; // z-position of the kth loop

    // Gauss-Legendre quadrature over θ (0 to 2π)
    for (int i = 0; i < GAUSS_POINTS; i++) {
      double theta = M_PI * (gauss_legendre_x[i] + 1.0); // Map to [0, 2π]
      double w = gauss_legendre_w[i];

      // Current element position in cylindrical coordinates
      double r_k = R;
      double phi_k = theta;
      double z_k = zk;

      // Observation point in cylindrical coordinates
      double r_p = r_obs;
      double phi_p = phi_obs;
      double z_p = z_obs;

      // Compute distance components
      double r_diff = r_p - r_k * cos(phi_k - phi_p);
      double z_diff = z_p - z_k;
      double r_sq = r_diff * r_diff + z_diff * z_diff;
      double r_mag = sqrt(r_sq);

      if (r_mag < 1e-6)
        continue; // Avoid singularities

      // Biot-Savart law: dB = (μ₀/4π) * (Idl × r) / r³
      double dl_phi = 1.0; // dl along φ direction
      double dBr = -sin(phi_k - phi_p) * dl_phi / (r_sq * r_mag);
      double dBphi = cos(phi_k - phi_p) * dl_phi / (r_sq * r_mag);
      double dBz = (r_k / r_mag) * dl_phi / (r_sq * r_mag);

      // Sum weighted contributions
      *Br += MU_0 / (4 * M_PI) * I_dz * w * dBr;
      *Bphi += MU_0 / (4 * M_PI) * I_dz * w * dBphi;
      *Bz += MU_0 / (4 * M_PI) * I_dz * w * dBz;
    }
  }
}

// Function to compute Biot-Savart law integral using Gaussian quadrature
// void compute_magnetic_field(double R, double L, double I, int N, double x,
//                             double y, double z, double *Bx, double *By,
//                             double *Bz) {
//   *Bx = 0.0;
//   *By = 0.0;
//   *Bz = 0.0;
//
//   double dtheta = 2.0 * M_PI / GAUSS_POINTS; // Step size in theta
//   double dz = L / N;                         // Step size along solenoid axis
//   double I_dz = I * dz;                      // Current element contribution
//
//   // Loop over solenoid turns
//   for (int k = 0; k < N; k++) {
//     double zk = -L / 2.0 + k * dz; // z-position of the kth loop
//
//     // Loop over Gauss quadrature points
//     for (int i = 0; i < GAUSS_POINTS; i++) {
//       double theta = M_PI * (gauss_legendre_x[i] + 1.0); // Map to [0, 2π]
//       double w = gauss_legendre_w[i];
//
//       // Current element position
//       double xk = R * cos(theta);
//       double yk = R * sin(theta);
//
//       // Vector r from current element to observation point
//       double rx = x - xk;
//       double ry = y - yk;
//       double rz = z - zk;
//       double r_squared = rx * rx + ry * ry + rz * rz;
//       double r = sqrt(r_squared);
//
//       if (r < 1e-6)
//         continue; // Avoid singularities
//
//       // Biot-Savart law: dB = (μ₀/4π) * (Idl × r) / r³
//       double dlx = -sin(theta) * R * dtheta;
//       double dly = cos(theta) * R * dtheta;
//       double dlz = 0.0;
//
//       // Cross product dl × r
//       double dBx = (dly * rz - dlz * ry) / (r_squared * r);
//       double dBy = (dlz * rx - dlx * rz) / (r_squared * r);
//       double dBz = (dlx * ry - dly * rx) / (r_squared * r);
//
//       // Sum weighted contributions
//       *Bx += MU_0 / (4 * M_PI) * I_dz * w * dBx;
//       *By += MU_0 / (4 * M_PI) * I_dz * w * dBy;
//       *Bz += MU_0 / (4 * M_PI) * I_dz * w * dBz;
//     }
//   }
// }
