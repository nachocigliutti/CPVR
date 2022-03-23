function [p, sol] = model_solve(p)

%----------------------------------------------------------------------
%------------------------         PARAMETERS        --------------------------
%----------------------------------------------------------------------

z = p.z;        % unemployment flow value
a = p.a;        % rent sharing. Range = (0, .35)
b = p.b;        % discount rate

% Contact rates
d = p.d;        % Separation rate (BLS = .04 - quarterly -)
l0 = p.l0;      % UE contact rate (Jarosch's = .06 - monthly -)
l1 = p.l1;      % EE contact rate (Jarosch's = .08 - monthly -)

% Firm's productivity (log Normal)
mu = p.mu;
s2 = p.s2;
n0 = p.n0;

% Approximate log Normal with qnwnorm
[x, w] = qnwnorm(n0, mu, s2);

% Trimmed state
lb = 1e-4;
x = x(w > lb); 
w = w(w > lb); 
w = w / sum(w);
y = exp(x);
n = length(y);

G = logncdf(y, mu, s2);
Gc = 1 - G;

% Parameters to initialize VF iteration
tol = 1e-7;    % tolerance
dif = 1;         % dif between S's
iter = 0;       % iteration counter

%----------------------------------------------------------------------
%------------------         VALUE FUNCTION ITERATION        -------------------
%----------------------------------------------------------------------

S = (y - z)  / (1 - b * d);
div = 1 - b * (1 - d) *(1 - a * l1 * Gc);

while dif > tol
    S = max(S, 0);
    Su = repmat(S, 1, n);
    Se = tril(Su);
    TS = ((y - z) + b * a * (1 - d) * (w' * (l1 * Se - l0 * Su))') ./ div;
    TS = max(TS, 0);

    dif = norm(abs(TS - S));
    iter = iter + 1;
    S = TS;

    if mod(iter, 50) == 0
        fprintf('Iteration = %2.0f, distance = %4.7f', iter, dif)
        fprintf('\n')
        disp('---------------------------------------------')
        fprintf('\n')
    end
    
    if dif <= tol
        fprintf('Model solved in %2.0f iterations', iter)
        fprintf('\n')
        disp('---------------------------------------------')
        fprintf('\n')
    end
end

% Value of unemployment
U = (z + b * l0 * a * w(S >= 0)' * S(S >= 0)) / (1 - b);

% Value of being employed at y (row) with benchmark y_ (col)
W = tril(U + repmat(S', n, 1) + a * (repmat(S, 1, n) - repmat(S', n, 1)));

% Choice sets
M1 = zeros(n, n, n);
M2 = zeros(n, n, n);

for i = 1:n          % y:   current
    for j = 1:n      % y_: benchmark
        for k = 1:n  % y': offer
            M1(i, j, k) = (k > i);                  % movers
            M2(i, j, k) = (k <= i) * (k > j);  % stayers + wage increase
        end
    end
end

% Wages
cv = zeros(n);  % continuation value

p1 = zeros(n);  
p2 = zeros(n);  
p3 = zeros(n);  

for i = 1:n         % current 
    for j = 1:n     % benchmark
        m1 = squeeze(M1(i, j, :));
        m2 = squeeze(M2(i, j, :));
        p1(i, j) = l1 * w' * (m1 .* W(:, i));
        p2(i, j) = l1 * w' * (m2 .* W(i, :)');
        p3(i, j) = (1 - l1 * w' * (m1 + m2)) * W(i, j);
        cv(i, j) = l1 * w' * (m1 .*  W(:, i) + m2 .* W(i, :)') + (1 - l1 * w' * (m1 + m2)) * W(i, j);
    end
end

% Wage of worker at y (row) with benchmark y_ (col)
wage = tril(W - b * d * U - b * (1 - d) * cv);
wage(wage == 0) = nan;
wage(wage < 0) = 1;

% Min productivity to accept offer
nU = find(S > 0, 1, 'first');

% Update structures
p.w = w;
p.y = y;
p.n = n;
p.lb = lb;

sol.U = U;
sol.W = W;
sol.S = S;
sol.wage = log(wage);
sol.nU = nU;
sol.p1 = p1;
sol.p2 = p2;
sol.p3 = p3;
end