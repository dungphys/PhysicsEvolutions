# **FKPP Solver**

## Theoretical summary
This solver is to solve the following 1D FKPP equation numerically:
$$
\partial_t u = D\partial^2_x u + f(u),
$$
where $D$ is the diffusion constant, $f(u)$ is the reaction function.

### **Reaction kernels**

1. Canonical FKPP kernel $~(r>0)$
$$
f(u) = ru(1-u), 
$$

2. Zeldovich kernel $(r>0, \theta \ge 0)$
$$
f(u) = 0, ~~ (x\le \theta), \\
f(u) = ru(1-u), ~~ (x>\theta)
$$

3. Power-law kernel $(r>0,k>1)$
$$
f(u) = ru^k(1-u), 
$$

4. Weak-alle kernel $(r>0,\delta>0)$
$$
f(u) = r u^{1 + \delta}(1 - u),
$$

5. Quadratic kernel $(r>0)$
$$
f(u) = r u(1 - u^n),
$$

### **Initial conditions**

1. Heaviside step function
$$
u(t=0,x) = 1-\Theta(x)  
$$

2. Exponetial-type function
$$
u(t=0,x) = 1-\exp(-e^{-kx})  
$$

3. Logistic function
$$
u(t=0,x) = \frac{1}{1 + \exp(-kx)}
$$

## Numerical schemes

### **Lattice $(t,x)$**

$$
NX {\rm ~points}:~~ X_{min} \le x \le X_{max}, ~{\rm spacing:}~dX,\\
NT {\rm ~points}:~~ 0\le t \le T_{max}, ~{\rm spacing:}~dT
$$

### **Explicit Euler scheme (conditionally stable)**

$$
u(t+dT, x) \\ = u(t, x) + \underbrace{D\frac{dT}{dX^2} \left[u(t, x-dX) - 2u(t, x) + u(t, x+dX)\right]}_{diffusion} + \underbrace{dTf(u(x))}_{reaction}
$$

Stability condition: $D\frac{dT}{dX^2} \le 0.5$

### **Implicit Euler scheme (unconditionally stable)**

$$
u(t+dT, x) - D\frac{dT}{dX^2}\left[u(t+dT, x-dX) - 2u(t+dT, x) + u(t+dT, x+dX)\right] \\[1.2em] 
= u(t,x) + dT f\left(u(t+dT, x)\right)
$$

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
./fkpp_solver
```

## Project Structure
```
fkpp_solver/
├── CMakeLists.txt
├── requirements.txt
├── README.md
├── plot/
|   ├── figures/
│   ├── front_pos.py     /*front position*/ 
│   ├── plot_front.py    /*plot solutions at different times*/
│   ├── py_requirements.txt
|   └── utils.cpp
├── src/
│   ├── main.cpp
│   ├── solver.cpp
│   ├── init.cpp
│   └── utils.cpp
├── include/
│   ├── solver.h  
│   ├── init.h
│   └── utils.h
└── bin/
    ├── grid_t.txt
    ├── grid_x.txt
    └── solution_x.xxx.txt
        ....
```