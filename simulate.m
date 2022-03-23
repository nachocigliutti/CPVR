function [hist, stats] = simulate(p, sol)

% Set seed
rng(1408)

% Unpack parameters
z = p.z;
a = p.a;
b = p.b; 
d = p.d;
l0 = p.l0; 
l1 = p.l1;
mu = p.mu;
s2 = p.s2;
n = p.n;
w = p.w;
y = p.y;
lb = p.lb;
nT = p.nT;
nW = p.nW;
sW = p.sW;

% Unpack value and policy functions
U = sol.U;
W = sol.W;
S = sol.S;
wage = sol.wage;
nU = sol.nU;

% Store history of jobs and wages
hist_e = zeros(nW, nT);     % employment status
hist_y = zeros(nW, nT);     % current firm (=0 if unemployed)
hist_b = zeros(nW, nT);     % benchmark firm (=0 if unemployed)
hist_w = zeros(nW, nT);     % wage history
hist_f = zeros(nW, nT);     % #firms

cdf = [0; cumsum(w)];
step = 0:(1/n):1;

%----------------------------------------------------------------------
%-------------------------         SIMULATION        -------------------------
%----------------------------------------------------------------------

% Sample potential wage offers for all workers
val = interp1(cdf, step, rand(nW, 1)); 
ind = ceil(val * n);

% Start most workers with a wage offer (95%)
e0 = (rand(nW, 1) <= .95) .* (ind >= nU);
hist_e(:, 1) = e0;
hist_y(:, 1) = e0 .* ind;
hist_b(:, 1) = e0 .* nU;
hist_w(:, 1) = e0 .* wage(ind, nU) + (1 - e0) * z;
hist_f(:, 1) = e0;

for t = 2:nT
    
    % Sample shocks
    eps_d = (rand(nW, 1) <= d);
    eps_l0 = (rand(nW, 1) <= l0);
    eps_l1 = (rand(nW, 1) <= l1);

    % Potential wage offers for all workers
    val = interp1(cdf, step,rand(nW, 1)); 
    ind = ceil(val * n);
   
    % Partition of the possible cases (destruction, offers and 'no offer')
    dest = hist_e(:, t-1) .* eps_d;
    offerU = (1 - hist_e(:, t-1)) .* eps_l0;
    offerE = hist_e(:, t-1) .* (1 - eps_d) .* eps_l1;
    offerN = hist_e(:, t-1) .* (1 - eps_d) .* (1 - eps_l1) + (1 - hist_e(:, t-1)) .* (1 - eps_l0);

    % Cases with potential changes in states
    move = offerE .* (ind > hist_y(:, t-1)) + offerU .* (ind >= nU);
    barg = offerE .* (ind <= hist_y(:, t-1)) .* (ind > hist_b(:, t-1));
    same = offerE .* (ind <= hist_b(:, t-1)) + offerU .* (ind < nU) + offerN;

    % Update histories
    hist_e(:, t) = hist_e(:, t-1) .* (1 - dest) + (1 - hist_e(:, t-1)) .* move;
    hist_y(:, t) = move .* ind + hist_e(:, t) .* (barg + same) .* hist_y(:, t-1);
    hist_b(:, t) = hist_e(:, t) .* move .* hist_y(:, t-1) + barg .* ind + ...
                            same .* hist_b(:, t-1) + (1 - hist_e(:, t-1)) .* move .* nU;
    hist_w(:, t) = hist_e(:, t) .* wage(sub2ind(size(wage), max(hist_y(:, t), 1), max(hist_b(:, t), 1))) + ...
                            (1 - hist_e(:, t)) .* z;
    hist_f(:, t) = hist_f(:, t-1) + move;
end

per = 10:10:100;
per_idx = [0 per];

% NO WORKER HETEROGENEITY (p<30 concentrated around UI)
% % Construct percentiles
% per_ini = prctile(hist_w(:, 1), per);
% per_end = prctile(hist_w(:, end), per);
% 
% % Find MIN to account for bunching at the bottom
% ini_min = find(per_ini > min(per_ini), 1, 'first');
% end_min = find(per_end > min(per_end), 1, 'first');
% 
% per_ini = [0 per_ini];
% per_end = [0 per_end];
% 
% Assign workers to income deciles
% val_ini = interp1(per_ini, per_idx, hist_w(:, 1));
% val_end = interp1(per_end, per_idx, hist_w(:, end));
% val_ini = ceil(val_ini / 10);
% val_end = ceil(val_end / 10);
% 
% % Compute statistics - NO worker heterogeneity
% med_w_ini = zeros(10, nT);
% med_w_end = zeros(10, nT);
% med_f_ini = zeros(10, nT);
% med_f_end = zeros(10, nT);
% 
% for p = 1:10
%     med_w_ini(p, :) = median(hist_w(val_ini == p, :));
%     med_w_end(p, :) = median(hist_w(val_end == p, :));
%     med_f_ini(p, :) = median(hist_f(val_ini == p, :));
%     med_f_end(p, :) = median(hist_f(val_end == p, :));
% end

%------------------------------------

% WORKER HETEROGENEITY (everything linear in x --> scale wages)

% Worker het. is fixed over time
eps_p = sW * rand(nW, 1);
hist_w = repmat(eps_p, 1, nT) .* hist_w;

per_ini = prctile(hist_w(:, 1), per);
per_end = prctile(hist_w(:, end), per);
per_mean = prctile(mean(hist_w, 2), per);

% Assign workers to income deciles
val_ini = interp1([0 per_ini], per_idx, hist_w(:, 1));
val_end = interp1([0 per_end], per_idx, hist_w(:, end));
val_mean = interp1([0 per_mean], per_idx, mean(hist_w, 2));

val_ini = ceil(val_ini / 10);
val_end = ceil(val_end / 10);
val_mean = ceil(val_mean / 10);

%------------------------------------

% Compute stats

% Medians
med_w_ini = zeros(10, nT); 
med_p_ini = zeros(10, nT); 
med_f_ini = zeros(10, nT); 

med_w_end = zeros(10, nT);
med_p_end = zeros(10, nT);
med_f_end = zeros(10, nT);

% Initial wage
mean_w_ini = zeros(10, nT); 
mean_p_ini = zeros(10, nT); 
unemp_ini = zeros(10, nT); 
mean_f_ini = zeros(10, nT);

% Final wage
mean_w_end = zeros(10, nT);
mean_p_end = zeros(10, nT);
unemp_end = zeros(10, nT); 
mean_f_end = zeros(10, nT);

% Mean wage
mean_w_mean = zeros(10, nT); 
mean_p_mean = zeros(10, nT);
unemp_mean = zeros(10, nT);
mean_f_mean = zeros(10, nT); 

for p = 1:10
    med_w_ini(p, :) = median(hist_w(val_ini == p, :));
    med_w_end(p, :) = median(hist_w(val_end == p, :));
    mean_w_ini(p, :) = mean(hist_w(val_ini == p, :));
    mean_w_end(p, :) = mean(hist_w(val_end == p, :));
    mean_w_mean(p, :) = mean(hist_w(val_mean == p, :));

    med_f_ini(p, :) = median(hist_f(val_ini == p, :));
    med_f_end(p, :) = median(hist_f(val_end == p, :));
    mean_f_ini(p, :) = mean(hist_f(val_ini == p, :));
    mean_f_end(p, :) = mean(hist_f(val_end == p, :));
    mean_f_mean(p, :) = mean(hist_f(val_mean == p, :));

    unemp_ini(p, :) = sum(1 - hist_e(val_ini == p, :)) ./ sum(val_ini == p);
    unemp_end(p, :) = sum(1 - hist_e(val_end == p, :)) ./ sum(val_end == p);
    unemp_mean(p, :) = sum(1 - hist_e(val_mean == p, :)) ./ sum(val_mean == p);

    for t = 1:nT
        med_p_ini(p, t) = median(hist_y(val_ini == p & hist_e(:, t) == 1, t));
        med_p_end(p, t) = median(hist_y(val_end == p & hist_e(:, t) == 1, t));
        mean_p_ini(p, t) = mean(hist_y(val_ini == p & hist_e(:, t) == 1, t));
        mean_p_end(p, t) = mean(hist_y(val_end == p & hist_e(:, t) == 1, t));
        mean_p_mean(p, t) = mean(hist_y(val_mean == p & hist_e(:, t) == 1, t));
    end
end

% Save everything in cells
hist = cell(4, 2);
hist{1, 1} = 'Current firm';
hist{2, 1} = 'Benchmark firm';
hist{3, 1} = 'Wages';
hist{4, 1} = 'Cum. # firms';

hist{1, 2} = hist_y;
hist{2, 2} = hist_b;
hist{3, 2} = hist_w;
hist{4, 2} = hist_f;

stats = cell(4, 5);
stats{1, 1} = 'Wages';
stats{2, 1} = 'Firm productivity';
stats{3, 1} = 'Unemployment rate';
stats{4, 1} = '# firms';

stats{1, 2} = med_w_ini; stats{1, 3} = med_w_end; stats{1, 4} = mean_w_ini; stats{1, 5} = mean_w_end; stats{1, 6} = mean_w_mean;
stats{2, 2} = med_p_ini; stats{2, 3} = med_p_end; stats{2, 4} = mean_p_ini; stats{2, 5} = mean_p_end; stats{2, 6} = mean_p_mean;
stats{3, 2} = unemp_ini; stats{3, 3} = unemp_end; stats{3, 4} = unemp_ini;   stats{3, 5} = unemp_end;  stats{3, 6} = unemp_mean;
stats{4, 2} = med_f_ini; stats{4, 3} = med_f_end; stats{4, 4} = mean_f_ini; stats{4, 5} = mean_f_end; stats{4, 6} = mean_f_mean;

%T = cell2table(stats, 'VariableNames',{'Stat', 'BOC - Median', 'EOC - Median', 'BOC - Mean', 'EOC - Mean', 'Mean Career W'});

% Colum names
colnames = cell(11, 1);
colnames{1} = 'age';
for i = 1:10
    colnames{i+1} = strcat('g', num2str(i));
end

% 1st column
ages = 20 + (0:(nT-1))';

% Files
names = {'BOC', 'EOC', 'avgW'};

for i = 1:3
    T1 = splitvars(table(ages, stats{1, 2+i}'));
    T2 = splitvars(table(ages, stats{2, 2+i}'));
    T3 = splitvars(table(ages, stats{3, 2+i}'));
    T4 = splitvars(table(ages, stats{4, 2+i}'));

    T1.Properties.VariableNames = colnames;
    T2.Properties.VariableNames = colnames;
    T3.Properties.VariableNames = colnames;
    T4.Properties.VariableNames = colnames;
    
    writetable(T1, strcat(names{i}, '_model_wage','.txt'), 'Delimiter',' ');
    writetable(T2, strcat(names{i}, '_model_prod','.txt'), 'Delimiter',' ');
    writetable(T3, strcat(names{i}, '_model_urate','.txt'), 'Delimiter',' ');
    writetable(T4, strcat(names{i}, '_model_Nfirm','.txt'), 'Delimiter',' ');
end