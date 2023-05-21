% Simulation setting
% CFA
% N
% J
% K
% within == 1
% maxiter
% correlation

person = 200;% 500, 1000, 2000;
item = 55;
domain = 5;
r = "high"; % "low";
within = 1; % 0;
maxiter = 10; %, 50;
Nrep = 100;
lr = 0.01; %, 0.1, 1;
S = 10; % 10, 50;
M = 10; % 10, 50;

simulation_IS_parallel(person, item, domain, r, within, maxiter, Nrep, S, M);


