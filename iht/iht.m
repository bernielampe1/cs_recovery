% Implementation for IHT sparse signal recovery
% https://arxiv.org/pdf/0805.0510.pdf

% INPUTS:
% y = compressed samples
% A = sampling matrix
% k = sparsity
% itrs = maximum number of iterations
% tol = algorithm termination tolerance

% OUTPUTS
% s = reconstructed signal
% r_mag = magnitude of residual at each iteration

function [xhat, r_mag] = iht(y, A, k, itrs, tol)

    % sizes of compressed and uncompressed signal
    [~, N] = size(A);

    % initialize the residual
    r = y;

    % initial solution
    xhat = zeros(N, 1);

    % init residual magnitudes
    r_mag = zeros(1, itrs);

    % loop for max iterations
    for itr = 1:itrs
        % step the solution in the direction vector from residual
        beta = xhat + A.' * r;

        % find the k-th largest coefficient
        threshold = min(maxk(abs(beta), k));

        % threshold
        xhat = beta;
        xhat(abs(beta) < threshold) = 0;

        r_mag(itr) = norm(r);
        fprintf(1, "itr = %d, residual = %f\n", itr, norm(r));

        % update residual
        r = y - A * xhat;

        % stopping criteria
        if norm(r) / norm(y) < tol; break; end
    end
end
