% NAMF_examples.m
% Usage examples for NAMF functions
% Run this script after configuring NAMF with NAMF_setup()

%% Example 1: Root finding with bisection
fprintf('=== EXAMPLE 1: Root finding ===\n');
fprintf('Finding sqrt(2) using the bisection method\n');

% Define the function f(x) = x^2 - 2
f = 'x^2 - 2';
[root, fval, iterations] = Bisect(f, 1, 2, 1e-10);

fprintf('Root found: %.10f\n', root);
fprintf('Function value: %.2e\n', fval);
fprintf('Number of iterations: %d\n', iterations);
fprintf('Error compared to sqrt(2): %.2e\n\n', abs(root - sqrt(2)));

%% Example 2: Numerical integration
fprintf('=== EXAMPLE 2: Numerical integration ===\n');
fprintf('Computing double integral of f(x,y) = x*y over [0,1]x[0,1]\n');

% Define the function to integrate
f_integrand = @(x,y) x.*y;

% Calculate the integral with different precision levels
m_values = [5, 10, 20];
for i = 1:length(m_values)
    m = m_values(i);
    result = trap_comp_integr_dop(0, 1, 0, 1, m, m, f_integrand);
    exact = 0.25; % Exact value of the integral
    error = abs(result - exact);
    
    fprintf('m = %2d: Result = %.6f, Error = %.2e\n', m, result, error);
end
fprintf('Exact value: %.6f\n\n', exact);

%% Example 3: Linear systems
fprintf('=== EXAMPLE 3: Linear systems ===\n');
fprintf('LU factorization of a 3x3 matrix\n');

% Define a test matrix
A = [4, 3, 2; 6, 3, 1; 8, 7, 5];
fprintf('Matrix A:\n');
disp(A);

% Calculate LU factorization
[L, U] = Fatt_LU(A);

fprintf('Matrix L (lower triangular):\n');
disp(L);
fprintf('Matrix U (upper triangular):\n');
disp(U);

% Verification: L*U should equal A (up to permutations)
product = L * U;
fprintf('Product L*U:\n');
disp(product);
fprintf('Error ||A - L*U||: %.2e\n\n', norm(A - product, 'fro'));

%% Example 4: Iterative methods
fprintf('=== EXAMPLE 4: Iterative methods ===\n');
fprintf('Solving system Ax = b with Jacobi and Gauss-Seidel\n');

% Test system
A = [4, -1, 0; -1, 4, -1; 0, -1, 4];
b = [1; 2; 3];
x0 = [0; 0; 0]; % Initial guess
tol = 1e-6;
max_iter = 100;

fprintf('System: Ax = b\n');
fprintf('A =\n'); disp(A);
fprintf('b = [%.1f; %.1f; %.1f]\n', b(1), b(2), b(3));

% Exact solution for comparison
x_exact = A \ b;
fprintf('Exact solution: [%.6f; %.6f; %.6f]\n', x_exact(1), x_exact(2), x_exact(3));

try
    % Jacobi method
    [x_jacobi, iter_jacobi] = Jacobi(A, b, x0, tol, max_iter);
    error_jacobi = norm(x_jacobi - x_exact);
    fprintf('Jacobi - Iterations: %d, Error: %.2e\n', iter_jacobi, error_jacobi);
    
    % Gauss-Seidel method
    [x_gs, iter_gs] = GaussSeidel(A, b, x0, tol, max_iter);
    error_gs = norm(x_gs - x_exact);
    fprintf('Gauss-Seidel - Iterations: %d, Error: %.2e\n', iter_gs, error_gs);
catch ME
    fprintf('Error in executing iterative methods: %s\n', ME.message);
end

fprintf('\n=== End of NAMF examples ===\n');