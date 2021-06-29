rng(1);
rng_gauss = normrnd(0,1,1,100000);
rng_idx = ceil(141*142*unifrnd(0,1,1,1500));

csvwrite('rng_gauss.csv',rng_gauss);
csvwrite('rng_idx.csv',rng_idx);