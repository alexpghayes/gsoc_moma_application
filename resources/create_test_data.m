%% TEST 1

rng(1);
load fisheriris.mat;
x = meas;
clear meas species;
[n, p] = size(x);

lamu = 1;
lamv = 2;

alphau = 3;
alphav = 4;

Omegu = rand(n, n);
Omegu = 0.5 * (Omegu + Omegu') + n * eye(n);

Omegv = rand(p, p);
Omegv = 0.5 * (Omegv + Omegv') + p * eye(p);

startu = 0;
startv = 0;
posu = 0;
posv = 0;
maxit = 100;

[U, V, d, ~] = sfpca_fixed(x, lamu, lamv, alphau, alphav, Omegu, ...
    Omegv, startu, startv, posu, posv, maxit);

save('../tests/testthat/test1.mat');
clear;

%% TEST 2

rng(2);
load fisheriris.mat;
x = meas;
clear meas species;
[n, p] = size(x);

lamu = 1;
lamv = 1;

alphau = 2;
alphav = 3;

Omegu = eye(n);
Omegv = eye(p);

startu = 0;
startv = 0;
posu = 0;
posv = 0;
maxit = 100;

[U, V, d, ~] = sfpca_fixed(x, lamu, lamv, alphau, alphav, Omegu, ...
    Omegv, startu, startv, posu, posv, maxit);

save('../tests/testthat/test2.mat');
clear;

%% TEST 3

rng(2);
load fisheriris.mat;
x = meas;
clear meas species;
[n, p] = size(x);

lamu = 1.7;
lamv = 0.2;

alphau = 2;
alphav = 0;

Omegu = rand(n, n);
Omegu = 0.5 * (Omegu + Omegu') + n * eye(n);

Omegv = rand(p, p);
Omegv = 0.5 * (Omegv + Omegv') + p * eye(p);

startu = 0;
startv = 0;
posu = 0;
posv = 0;
maxit = 100;

[U, V, d, ~] = sfpca_fixed(x, lamu, lamv, alphau, alphav, Omegu, ...
    Omegv, startu, startv, posu, posv, maxit);

save('../tests/testthat/test3.mat');
clear;