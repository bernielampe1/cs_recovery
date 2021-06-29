% Simple test to demonstrate IHT reconstruction
function sse = iht_dense_test()
    k = 50;          % sparsity (number of non-zero dct coeff)
    n = 2000;        % dimension of signal vector
    amps = [-100, 100];  % amplitude of dct coeff

    % number of samples needed
    c = 10;  % NOTE: if this isn't large enough, the iteration will diverge
    m = ceil(c * k * log(n/k));   % number of samples needed
    A = normc(randn(m, n));       % sampling matrix

    % generate signal
    signal_f = zeros(n, 1);
    signal_f(randi([1, n], k, 1)) = randi(amps, k, 1);
    signal_t = idct(signal_f);

    % sample the signal
    y = A * signal_t;

    % iht reconstruction
    xp = iht(y, A * dctmtx(n)', k, 500, 1e-10);

    % compute the error
    sse = sum((signal_f - xp).^2);

    % plot it
    figure;
    stem(signal_f);
    hold on;
    stem(xp);
    title('plot of signal and reconstruction');

    figure;
    stem(abs(signal_f - xp));
    title('absolute error on each sample');
end

function [M, n] = normc(M)
   n = sqrt(sum(M.^2,1));     % Compute norms of columns
   M = bsxfun(@rdivide,M,n);  % Normalize M
   n = reshape(n,[],1);       % Store column vector of norms
end
