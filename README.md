# NAMF - Numerical Analysis MATLAB Functions

A comprehensive collection of MATLAB functions for numerical analysis, organized by mathematical method categories.

## 📋 Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Function Categories](#function-categories)
- [Installation and Usage](#installation-and-usage)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## 🎯 Overview

NAMF is a collection of MATLAB implementations of fundamental numerical analysis algorithms, developed for educational and research purposes. The functions cover three main areas:

- **Quadrature Formulas**: Numerical methods for computing definite integrals
- **Root Finding**: Iterative algorithms for finding function zeros
- **Linear Systems**: Direct and iterative methods for system solving

## 📁 Project Structure

```
NAMF/
├── Quadrature_Formulas/           # Numerical integration methods
├── Root_Finding/                  # Root finding algorithms
├── Linear_Systems/                # Linear system solving
│   ├── Factorizations/           # Factorization methods
│   └── Iterative_Methods/        # Iterative methods
├── LICENSE
└── README.md
```

## 🔧 Function Categories

### 📊 Quadrature Formulas

Methods to estimate definite integral values without computing the antiderivative:

| Function | Description |
|----------|-------------|
| `trap_comp_integr_dop.m` | Double integration with composite trapezoids |
| `NewCot_chiu_semp.m` | Simple closed Newton-Cotes |
| `NewCot_chiu_comp.m` | Composite closed Newton-Cotes |
| `NewCot_aper_semp.m` | Simple open Newton-Cotes |
| `NewCot_aper_comp.m` | Composite open Newton-Cotes |
| `NCcc_Num_iter.m` | Newton-Cotes iteration number calculation |
| `NewCot_chiu_semp_pol_carat.m` | Simple closed Newton-Cotes for characteristic polynomial |

### 🎯 Root Finding Functions

Iterative algorithms to find α ∈ [a,b] such that f(α) = 0:

| Function | Description |
|----------|-------------|
| `Bisect.m` | Bisection method |
| `rad_newton_vett.m` | Vectorial Newton method |
| `rad_pol_newton.m` | Newton method for polynomials |
| `Horner.m` | Horner's algorithm |
| `Frobenius.m` | Frobenius method |
| `rad_minmax.m` | Min/max root search |
| `max_iter_trapezi.m` | Maximum iterations for trapezoidal rule |
| `newton_mat_trid_sim.m` | Newton method for symmetric tridiagonal matrix eigenvalues |
| `val_pol_carat_trid.m` | Characteristic polynomial evaluation for tridiagonal matrices |

### 🔢 Linear Systems

#### Factorizations
| Function | Description |
|----------|-------------|
| `Fatt_LU.m` | LU factorization with pivoting |
| `Fatt_LU_NOPIV.m` | LU factorization without pivoting |
| `Fatt_QR.m` | QR factorization |

#### Iterative Methods
| Function | Description |
|----------|-------------|
| `Jacobi.m` | Jacobi method |
| `GaussSeidel.m` | Gauss-Seidel method |

#### System Solving
| Function | Description |
|----------|-------------|
| `RSL_SI.m` | Lower triangular system solving (forward substitution) |
| `RSL_SA.m` | Upper triangular system solving (backward substitution) |
| `RisolSisMatTrid.m` | Tridiagonal matrix system solving |

## 🚀 Installation and Usage

1. **Clone the repository:**
   ```bash
   git clone https://github.com/username/NAMF.git
   ```

2. **Add to MATLAB path:**
   ```matlab
   addpath(genpath('path/to/NAMF'))
   ```

3. **Use the functions:**
   ```matlab
   % Example: Bisection method
   f = 'x^2 - 2';
   [x, fx, n] = Bisect(f, 0, 2, 1e-6);
   
   % Example: Iterative methods
   A = [4 -1; -1 4]; b = [1; 2]; x0 = [0; 0];
   [x, iter] = Jacobi(A, b, x0, 1e-6, 100);
   ```

4. **Run tests:**
   ```matlab
   NAMF_tests  % Run comprehensive test suite
   ```

## 📝 Examples

### Numerical Integration
```matlab
% Double integral with composite trapezoids
f = @(x,y) x.*y + sin(x.*y);
result = trap_comp_integr_dop(0, 1, 0, 1, 10, 10, f);
```

### Root Finding
```matlab
% Bisection to find √2
f = 'x^2 - 2';
[root, fval, iter] = Bisect(f, 1, 2, 1e-10);
fprintf('√2 ≈ %.10f (iterations: %d)\n', root, iter);
```

### Linear Systems
```matlab
% LU factorization
A = [4 3; 6 3];
[L, U] = Fatt_LU(A);
```

## 🤝 Contributing

Contributions are welcome! To contribute:

1. Fork the project
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 License

This project is distributed under the MIT License. See the `LICENSE` file for more details.

---

**Note**: This project is developed for educational and research purposes in numerical analysis.