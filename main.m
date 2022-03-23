clear; clc; rng(1005);
cd('/Users/nacho/Dropbox/Job Ladder/CPVR');


%----------------------------------------------------------------------
%------------------------         PARAMETERS        --------------------------
%----------------------------------------------------------------------


aa = [0 .1 .2 .3 .35];
%aa = .35;

for i = 1:length(aa)

freq = 1; % .25 = monthly, 1 = quarterly, 4 = yearly

p.z = 1;     % unemployment flow value
p.a = aa(i);   % rent sharing. Range = (0, .35)
p.sW = 4;     % Dispersion of worker heterogeniety (around ~U[0,1]) 

% Calibration (quarterly) from Moscarini & Postel Vinay (2020)
p.b = .98 ^ freq;                        % discount rate
p.d = 1 - (1 - .024) ^ freq;        % Separation rate (BLS = .04 - quarterly -)
p.l0 = 1 - (1 - .45) ^ freq;        % UE contact rate (Jarosch's = .06 - monthly -)
p.l1 = 1 - (1 - .02) ^ freq;        % EE contact rate (Jarosch's = .08 - monthly -)


% Parametrize firm's productivity with logN
p.mu = log(200);
p.s2 = .4;
p.n0 = 60;

% Parameters for simulation
p.nW = 5e5;                          % # workers
p.nT = 36 * (4 / freq);        % # periods (35 years)


%----------------------------------------------------------------------
%-------------------------         RUN MODEL        ---------------------------
%----------------------------------------------------------------------


[p, sol] = model_solve(p);
[hist, stats] = simulate(p, sol);

run('plots.m')


end