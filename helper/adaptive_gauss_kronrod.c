#include <math.h>
#include <stdio.h>

long double TOLERANCE = 1e-6;

// Function to integrate (modify as needed)
long double f(long double x) { return 2 * x; }
// long double f(long double x) { return (pow(x, 24) * exp(-x) * cos(x)); }
// long double f(long double x) { return (exp(-x) * cos(x)); }

// Gauss weights
long double G10W[] = {2.955242247147528701738929946513383e-01L,
                      2.692667193099963550912269215694694e-01L,
                      2.190863625159820439955349342281632e-01L,
                      1.494513491505805931457763396576973e-01L,
                      6.667134430868813759356880989333179e-02L};

// Kronrod nodes
long double GK21N[] = {0.000000000000000000000000000000000e+00L,
                       1.488743389816312108848260011297200e-01L,
                       2.943928627014601981311266031038656e-01L,
                       4.333953941292471907992659431657842e-01L,
                       5.627571346686046833390000992726941e-01L,
                       6.794095682990244062343273651148736e-01L,
                       7.808177265864168970637175783450424e-01L,
                       8.650633666889845107320966884234930e-01L,
                       9.301574913557082260012071800595083e-01L,
                       9.739065285171717200779640120844521e-01L,
                       9.956571630258080807355272806890028e-01L};

// Kronrod weights
long double GK21W[] = {1.494455540029169056649364683898212e-01L,
                       1.477391049013384913748415159720680e-01L,
                       1.427759385770600807970942731387171e-01L,
                       1.347092173114733259280540017717068e-01L,
                       1.234919762620658510779581098310742e-01L,
                       1.093871588022976418992105903258050e-01L,
                       9.312545458369760553506546508336634e-02L,
                       7.503967481091995276704314091619001e-02L,
                       5.475589657435199603138130024458018e-02L,
                       3.255816230796472747881897245938976e-02L,
                       1.169463886737187427806439606219205e-02L};

long double adaptive_gauss_kronrod(long double (*f)(long double), long double a,
                                   long double b, long double tol) {
  long double c = (a + b) / 2.0;
  long double h = (b - a) / 2.0;
  // long double function_values[node];

  // Evaluate function at Gauss-Kronrod 21-point nodes
  // for (int i = 0; i < node; i++){
  //   function_values[i];
  // }
  long double f0 = f(c);
  long double f1 = f(c - h * GK21N[1]);
  long double f2 = f(c + h * GK21N[1]);
  long double f3 = f(c - h * GK21N[2]);
  long double f4 = f(c + h * GK21N[2]);
  long double f5 = f(c - h * GK21N[3]);
  long double f6 = f(c + h * GK21N[3]);
  long double f7 = f(c - h * GK21N[4]);
  long double f8 = f(c + h * GK21N[4]);
  long double f9 = f(c - h * GK21N[5]);
  long double f10 = f(c + h * GK21N[5]);
  long double f11 = f(c - h * GK21N[6]);
  long double f12 = f(c + h * GK21N[6]);
  long double f13 = f(c - h * GK21N[7]);
  long double f14 = f(c + h * GK21N[7]);
  long double f15 = f(c - h * GK21N[8]);
  long double f16 = f(c + h * GK21N[8]);
  long double f17 = f(c - h * GK21N[9]);
  long double f18 = f(c + h * GK21N[9]);
  long double f19 = f(c - h * GK21N[10]);
  long double f20 = f(c + h * GK21N[10]);

  // Compute 21-point Kronrod estimate
  long double I_kronrod =
      h * (GK21W[0] * f0 + GK21W[1] * (f1 + f2) + GK21W[2] * (f3 + f4) +
           GK21W[3] * (f5 + f6) + GK21W[4] * (f7 + f8) + GK21W[5] * (f9 + f10) +
           GK21W[6] * (f11 + f12) + GK21W[7] * (f13 + f14) +
           GK21W[8] * (f15 + f16) + GK21W[9] * (f17 + f18) +
           GK21W[10] * (f19 + f20));

  // Compute 10-point Gauss estimate
  long double I_gauss =
      h * (G10W[0] * (f1 + f2) + G10W[1] * (f6 + f5) + G10W[2] * (f9 + f10) +
           G10W[3] * (f13 + f14) + G10W[4] * (f17 + f18));

  if (TOLERANCE / tol == 1) {
    printf("Gauss: %Lf,\nKronrod: %Lf\n\n", I_gauss, I_kronrod);
  }

  // Estimate error
  long double error = fabsl(I_kronrod - I_gauss);

  // Check if within tolerance
  if (error < tol) {
    return I_kronrod;
  } else {
    // Subdivide interval
    long double left = adaptive_gauss_kronrod(f, a, c, tol / 2);
    long double right = adaptive_gauss_kronrod(f, c, b, tol / 2);
    return left + right;
  }
}
//
// int main() {
//   long double a = 0.0, b = 1.0; // Integration limits
//   long double tol = 1e-6;       // Tolerance for error
//   long double result = adaptive_gauss_kronrod(a, b, tol);
//
//   printf("Integral result: %.25Lf\n", result);
//   return 0;
// }
