clear; clc;

pairs = {
    'Sn(oxide)-Sn(melt)', 404423.180, 0.204925, 39.305410;
    'O(oxide)-Sn(melt)',   88619.439, 0.168992, 16.277801;
};

rcut = 12.0;
tol_newton = 1e-10;
tol_compare = 1e-5;
maxIter = 100;

options = optimset( ...
    'TolX', 1e-12, ...
    'Display', 'off' ...
);

results = cell(size(pairs, 1), 8);

for p = 1:size(pairs, 1)

    pairName = pairs{p, 1};
    A   = pairs{p, 2};
    rho = pairs{p, 3};
    C   = pairs{p, 4};

    UB = @(r) A .* exp(-r ./ rho) - C ./ r.^6;

    dUB = @(r) -(A ./ rho) .* exp(-r ./ rho) + 6 .* C ./ r.^7;

    d2UB = @(r) (A ./ rho.^2) .* exp(-r ./ rho) - 42 .* C ./ r.^8;

    r_newton = 2.5;

    for iter = 1:maxIter

        r_old = r_newton;

        first_derivative  = dUB(r_old);
        second_derivative = d2UB(r_old);

        if abs(second_derivative) < eps
            error('Second derivative is too small for Newton iteration: %s', pairName);
        end

        r_newton = r_old - first_derivative / second_derivative;

        if abs(r_newton - r_old) < tol_newton
            break;
        end

    end

    if iter == maxIter
        warning('Newton-Raphson iteration did not fully converge for %s.', pairName);
    end

    r_min = fminbnd(UB, 1.5, 5.0, options);

    diff_r0 = abs(r_newton - r_min);

    if diff_r0 > tol_compare
        warning('The two methods give different r0 values for %s.', pairName);
    end

    r0 = r_newton;

    D = -UB(r0);

    if D <= 0
        error('The calculated well depth D is non-positive for %s.', pairName);
    end

    curvature = d2UB(r0);

    if curvature <= 0
        error('The Buckingham potential minimum has non-positive curvature for %s.', pairName);
    end

    alpha = sqrt(curvature / (2 * D));

    results{p, 1} = pairName;
    results{p, 2} = A;
    results{p, 3} = rho;
    results{p, 4} = C;
    results{p, 5} = D;
    results{p, 6} = alpha;
    results{p, 7} = r0;
    results{p, 8} = rcut;

    fprintf('\nPair: %s\n', pairName);
    fprintf('  r0 from Newton-Raphson = %.10f Angstrom\n', r_newton);
    fprintf('  r0 from fminbnd        = %.10f Angstrom\n', r_min);
    fprintf('  difference             = %.3e Angstrom\n', diff_r0);
    fprintf('  D                      = %.6f eV\n', D);
    fprintf('  alpha                  = %.6f Angstrom^-1\n', alpha);
    fprintf('  rcut                   = %.3f Angstrom\n', rcut);

end

T = cell2table(results, ...
    'VariableNames', { ...
    'Pair', ...
    'A_eV', ...
    'rho_Angstrom', ...
    'C_eV_Angstrom6', ...
    'D_eV', ...
    'alpha_Angstrom_inv', ...
    'r0_Angstrom', ...
    'rcut_Angstrom'});

disp(' ');
disp('Derived Morse potential parameters:');
disp(T);

writetable(T, 'Buckingham_to_Morse_parameters.csv');

disp(' ');
disp('Results have been written to Buckingham_to_Morse_parameters.csv');