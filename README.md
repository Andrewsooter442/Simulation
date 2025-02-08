# Project Structure  
- The `helper/` folder contains supporting methods. Currently, only `adaptive_gauss_kronrod.c` is relevant.  
- The base directory contains test files for evaluating methods. You can add or remove functions in the test file as needed.  
- For now, only `adaptive_gauss_test.c` and `adaptive_gauss_kronrod.c` are important.  

## Running the Project  
### On macOS and Linux  
Compile and run the test program with:  
```sh
gcc adaptive_gauss_test.c -o prog.out && ./prog.out
