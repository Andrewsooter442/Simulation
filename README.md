# How the project is structured 
The methods are stored in the helper folder (For now only the adaptive_gauss_kronrod.c is usefull)

In the base directory there are tests which can be run to test the method on functions more functions can be added or removed in the test file.
For now only the adaptive_gauss_test.c file and adaptive_gauss_kronrod.c is important


### How to run
#### On mac and linux.

```gcc adaptive_gauss_test.c  -o prog.out && ./prog.out```
