# **BK Solver**

## Theoretical summary
This solver is to solve the following Balitsky-Kovchegov equation for the S-matrix element numerically:
$$
\partial_Y S(r,Y) = \int d^{2}r K(r,r_1,r_2) \left[S(r_1,Y)S(r_2,Y) - S(r,Y)\right]
$$
where $K$ is the kernel, and $S(r,Y)$ is the S-matrix element at dipole size $r$ and rapidity $Y$. Dipole amplitude: $N(r,Y) = 1 - S(r,Y)$.

### **Integral kernel kernels**

1. Leading-order kernel. 

2. Running-coupling kernel with parent-dipole prescription.

3. Running-coupling kernel with minimum-dipole prescription.

4. Running-coupling kernel with Balitsky prescription.

5. Running-coupling kernel with Kovchegov-Weigert prescription.

### **Initial conditions**

1. Step function
$$
S(r,Y=0) = 1, ~ r < r_0,\\
S(r,Y=0) = 0, ~ r \le r_0  
$$

2. Golec-Biernat-Wusthoff function
$$
S(r,Y=0) = \exp\left[-\left(\frac{r^2Q_{s0}^2}{4}\right)^{\gamma}\right]  
$$

3. McLerran-Venugopalan function
$$
S(r,Y=0) = \exp\left[-\left(\frac{r^2Q_{s0}^2}{4}\right)^{\gamma}\ln \left(e + \frac{1}{r\Lambda_{QCD}}\right)\right]  
$$

## Numerical schemes

- RK4


## Requirements
- C++11 or later
- CMake 3.10+
- Make

## Building

```bash
# Create build directory
mkdir build && cd build

# Configure
cmake ..

# Build
make

# Run
./bk_solver
```

## Project Structure
```
fkpp_solver/
├── CMakeLists.txt
├── README.md
├── src/
│   ├── main.cpp
│   ├── solver.cpp
│   ├── init_condition.cpp
│   ├── interp.cpp
│   └── utils.cpp
├── include/
│   ├── solver.h  
│   ├── init_condition.h
│   ├── interp.h
│   └── utils.h
└── bin/
    ├── grid_t.txt
    ├── grid_x.txt
    └── solution_x.xxx.txt
        ....
```