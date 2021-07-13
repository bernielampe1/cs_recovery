% Simple test to demonstrate StOMP reconstruction
function sse = stomp_test()
    k = 50;          % sparsity (number of non-zero dct coeff)
    n = 2000;        % dimension of signal vector
    amps = [-100, 100];  % amplitude of dct coeff

    % number of samples needed
    c = 4;  % NOTE: if this isn't large enough, the iteration will diverge
    m = ceil(c * k * log(n/k))    % number of samples needed
    A = normc(randn(m, n));       % sampling matrix

    % generate signal
    signal_f = zeros(n, 1);
    signal_f(randi([1, n], k, 1)) = randi(amps, k, 1);
    signal_t = idct(signal_f);

    % sample the signal
    y = A * signal_t;

    % stomp reconstruction
    xp = StOMP(A * dctmtx(n)', y, 5);

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
