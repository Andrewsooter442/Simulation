#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Function prototype
double func(double x) {
  return (pow(x, 24) + 3 * pow(x, 10) - 5 * pow(x, 8) + 2 * pow(x, 6) -
          pow(x, 4) + 7 * pow(x, 2)) /
         (pow(x, 2))*cos(x * x * x);
}
double gaussian_quadrature(int n, double (*f)(double), double a, double b);
void get_nodes_weights(int n, double **nodes, double **weights);

int main() {
  // Example: Integrating the function f(x) = x^2 over the interval [-1, 1]
  double result = 0;
  for (int i = 0; i < 100; i++) {
    result = gaussian_quadrature(5, func, 0, 2.0);
  }
  printf("Integral result: %.10f\n", result);
  return 0;
}

// Example function: f(x) = x^2

// Gaussian Quadrature for n-point rule
double gaussian_quadrature(int n, double (*f)(double), double a, double b) {
  // Allocate memory for nodes and weights
  double *nodes = (double *)malloc(n * sizeof(double));
  double *weights = (double *)malloc(n * sizeof(double));

  if (nodes == NULL || weights == NULL) {
    printf("Memory allocation failed\n");
    exit(1);
  }

  // Get the Gauss-Legendre nodes and weights
  get_nodes_weights(n, &nodes, &weights);

  double sum = 0.0;

  // Loop through nodes and compute weighted sum
  for (int i = 0; i < n; i++) {
    // Change of variable: map node from [-1, 1] to [a, b]
    double x = 0.5 * (a + b) + 0.5 * (b - a) * nodes[i];
    double w = 0.5 * (b - a) * weights[i];

    sum += w * f(x); // Add the weighted value of the function
  }

  // Free allocated memory
  free(nodes);
  free(weights);

  return sum;
}

void get_nodes_weights(int n, double **nodes, double **weights) {
  // Predefined nodes and weights for different n values
  static double nodes_arr[12][12] = {
      {0.0},
      {-0.5773502691896257, 0.5773502691896257},
      {-0.7745966692414834, 0.0, 0.7745966692414834},
      {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
       0.8611363115940526},
      {-0.9061798459386639, -0.5384693101056831, 0.0, 0.5384693101056831,
       0.9061798459386639},
      // Add more nodes here for 6, 7, ..., 12 point quadratures
      // Ensure the arrays are padded with zeros for unfilled positions
  };

  static double weights_arr[12][12] = {
      {2.0},
      {1.0, 1.0},
      {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0},
      {0.3478548451374539, 0.6521451548625461, 0.6521451548625461,
       0.3478548451374539},
      {0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
       0.4786286704993665, 0.2369268850561891},
      // Add corresponding weights for other cases
  };

  if (n < 1 || n > 12) {
    printf("Number of nodes: %d, out of range. Supported range: 1 to 12.\n", n);
    exit(1);
  }

  // Set the nodes and weights for the requested number of nodes
  for (int i = 0; i < n; i++) {
    (*nodes)[i] = nodes_arr[n - 1][i];
    (*weights)[i] = weights_arr[n - 1][i];
  }
}

// Function to get nodes and weights for Gaussian quadrature
// void get_nodes_weights(int n, double **nodes, double **weights) {
//   switch (n) {
//   case 1:
//     (*nodes)[0] = 0.0;
//     (*weights)[0] = 2.0;
//     break;
//
//   case 2:
//     (*nodes)[0] = -1.0 / sqrt(3.0);
//     (*nodes)[1] = 1.0 / sqrt(3.0);
//     (*weights)[0] = 1.0;
//     (*weights)[1] = 1.0;
//     break;
//
//   case 3:
//     (*nodes)[0] = -sqrt(3.0 / 5.0);
//     (*nodes)[1] = 0.0;
//     (*nodes)[2] = sqrt(3.0 / 5.0);
//     (*weights)[0] = 5.0 / 9.0;
//     (*weights)[1] = 8.0 / 9.0;
//     (*weights)[2] = 5.0 / 9.0;
//     break;
//
//   case 4:
//     (*nodes)[0] = -0.8611363115940526;
//     (*nodes)[1] = -0.3399810435848563;
//     (*nodes)[2] = 0.3399810435848563;
//     (*nodes)[3] = 0.8611363115940526;
//     (*weights)[0] = 0.3478548451374539;
//     (*weights)[1] = 0.6521451548625461;
//     (*weights)[2] = 0.6521451548625461;
//     (*weights)[3] = 0.3478548451374539;
//     break;
//
//   case 5:
//     (*nodes)[0] = -0.9061798459386639;
//     (*nodes)[1] = -0.5384693101056831;
//     (*nodes)[2] = 0.0;
//     (*nodes)[3] = 0.5384693101056831;
//     (*nodes)[4] = 0.9061798459386639;
//     (*weights)[0] = 0.2369268850561891;
//     (*weights)[1] = 0.4786286704993665;
//     (*weights)[2] = 0.5688888888888889;
//     (*weights)[3] = 0.4786286704993665;
//     (*weights)[4] = 0.2369268850561891;
//     break;
//
//   default:
//     printf("number of nodes: %d\n > 5", n);
//     exit(1);
//   }
// }
