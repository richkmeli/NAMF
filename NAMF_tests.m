% NAMF_tests.m
% Comprehensive test suite for NAMF functions
% Run this script to validate all functions work correctly

fprintf('=== NAMF COMPREHENSIVE TEST SUITE ===\n\n');

%% Test 1: Root Finding Functions
fprintf('1. Testing Root Finding Functions\n');
fprintf('----------------------------------\n');

% Test Bisection method
try
    f = 'x^2 - 2';
    [root, fval, iter] = Bisect(f, 1, 2, 1e-8);
    error = abs(root - sqrt(2));
    fprintf('✓ Bisect: root=%.8f, error=%.2e, iterations=%d\n', root, error, iter);
catch ME
    fprintf('✗ Bisect failed: %s\n', ME.message);
end

% Test Newton method for polynomials
try
    coeffp = [1 0 -2]; % x^2 - 2
    root = rad_pol_newton(coeffp, 1.5, 1e-8, 100);
    error = abs(root - sqrt(2));
    fprintf('✓ Newton polynomial: root=%.8f, error=%.2e\n', root, error);
catch ME
    fprintf('✗ Newton polynomial failed: %s\n', ME.message);
end

% Test Horner's algorithm
try
    coeffp = [1 0 -2]; % x^2 - 2
    x = sqrt(2);
    [y, dy] = Horner(coeffp, x);
    fprintf('✓ Horner: f(√2)=%.2e, f\'(√2)=%.6f\n', y, dy);
catch ME
    fprintf('✗ Horner failed: %s\n', ME.message);
end

%% Test 2: Quadrature Formulas
fprintf('\n2. Testing Quadrature Formulas\n');
fprintf('------------------------------\n');

% Test simple closed Newton-Cotes
try
    f = @(x) x.^2;
    exact = 8/3; % integral of x^2 from 0 to 2
    Q1 = NewCot_chiu_semp(f, 0, 2, 1); % Trapezoidal
    Q2 = NewCot_chiu_semp(f, 0, 2, 2); % Simpson
    fprintf('✓ Newton-Cotes (Trapezoid): %.6f, error=%.2e\n', Q1, abs(Q1-exact));
    fprintf('✓ Newton-Cotes (Simpson): %.6f, error=%.2e\n', Q2, abs(Q2-exact));
catch ME
    fprintf('✗ Newton-Cotes failed: %s\n', ME.message);
end

% Test composite Newton-Cotes
try
    f = @(x) x.^2;
    exact = 8/3;
    Q = NewCot_chiu_comp(f, 0, 2, 2, 4); % Simpson with 4 intervals
    fprintf('✓ Composite Newton-Cotes: %.6f, error=%.2e\n', Q, abs(Q-exact));
catch ME
    fprintf('✗ Composite Newton-Cotes failed: %s\n', ME.message);
end

% Test double integration
try
    f = @(x,y) x.*y;
    exact = 0.25; % integral of xy over [0,1]x[0,1]
    result = trap_comp_integr_dop(0, 1, 0, 1, 10, 10, f);
    fprintf('✓ Double integration: %.6f, error=%.2e\n', result, abs(result-exact));
catch ME
    fprintf('✗ Double integration failed: %s\n', ME.message);
end

%% Test 3: Linear System Factorizations
fprintf('\n3. Testing Linear System Factorizations\n');
fprintf('---------------------------------------\n');

% Test LU factorization
try
    A = [4 3 2; 6 3 1; 8 7 5];
    [L, U] = Fatt_LU(A);
    error = norm(L*U - A, 'fro');
    fprintf('✓ LU factorization: ||L*U - A|| = %.2e\n', error);
catch ME
    fprintf('✗ LU factorization failed: %s\n', ME.message);
end

% Test LU factorization without pivoting
try
    A = [2 1; 1 2];
    [L, U] = Fatt_LU_NOPIV(A);
    error = norm(L*U - A, 'fro');
    fprintf('✓ LU no pivoting: ||L*U - A|| = %.2e\n', error);
catch ME
    fprintf('✗ LU no pivoting failed: %s\n', ME.message);
end

% Test QR factorization
try
    A = [1 2; 3 4; 5 6];
    [Q, R] = Fatt_QR(A);
    error = norm(Q*R - A, 'fro');
    fprintf('✓ QR factorization: ||Q*R - A|| = %.2e\n', error);
catch ME
    fprintf('✗ QR factorization failed: %s\n', ME.message);
end

%% Test 4: Triangular System Solvers
fprintf('\n4. Testing Triangular System Solvers\n');
fprintf('------------------------------------\n');

% Test forward substitution (lower triangular)
try
    L = [2 0 0; 1 3 0; 4 2 5];
    b = [2; 5; 20];
    x = RSL_SI(L, b);
    error = norm(L*x - b);
    fprintf('✓ Forward substitution: ||L*x - b|| = %.2e\n', error);
catch ME
    fprintf('✗ Forward substitution failed: %s\n', ME.message);
end

% Test backward substitution (upper triangular)
try
    U = [2 1 4; 0 3 2; 0 0 5];
    b = [20; 5; 2];
    x = RSL_SA(U, b);
    error = norm(U*x - b);
    fprintf('✓ Backward substitution: ||U*x - b|| = %.2e\n', error);
catch ME
    fprintf('✗ Backward substitution failed: %s\n', ME.message);
end

% Test tridiagonal system solver
try
    M = [2 -1 0; -1 2 -1; 0 -1 2];
    b = [1; 0; 1];
    x = RisolSisMatTrid(M, b);
    error = norm(M*x - b);
    fprintf('✓ Tridiagonal solver: ||M*x - b|| = %.2e\n', error);
catch ME
    fprintf('✗ Tridiagonal solver failed: %s\n', ME.message);
end

%% Test 5: Iterative Methods
fprintf('\n5. Testing Iterative Methods\n');
fprintf('----------------------------\n');

% Test system with guaranteed convergence
A = [4 -1 0; -1 4 -1; 0 -1 4]; % Diagonally dominant
b = [1; 2; 3];
x_exact = A \ b;

% Test Jacobi method
try
    [x_jacobi, iter_jacobi] = Jacobi(A, b, zeros(3,1), 1e-8, 100);
    error_jacobi = norm(x_jacobi - x_exact);
    fprintf('✓ Jacobi: error=%.2e, iterations=%d\n', error_jacobi, iter_jacobi);
catch ME
    fprintf('✗ Jacobi failed: %s\n', ME.message);
end

% Test Gauss-Seidel method
try
    [x_gs, iter_gs] = GaussSeidel(A, b, zeros(3,1), 1e-8, 100);
    error_gs = norm(x_gs - x_exact);
    fprintf('✓ Gauss-Seidel: error=%.2e, iterations=%d\n', error_gs, iter_gs);
catch ME
    fprintf('✗ Gauss-Seidel failed: %s\n', ME.message);
end

%% Summary
fprintf('\n=== TEST SUITE COMPLETED ===\n');
fprintf('All major function categories have been tested.\n');
fprintf('Check above for any ✗ marks indicating failures.\n');