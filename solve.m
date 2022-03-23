% OLD MODEL SOLVER
%
% function [sol, p] = solve(p)
% 
% %----------------------------------------------------------------------
% %------------------------         PARAMETERS        --------------------------
% %----------------------------------------------------------------------
% 
% z = p.z;        % unemployment flow value
% a = p.a;        % rent sharing. Range = (0, .35)
% b = p.b;        % discount rate
% 
% % Contact rates
% d = p.d;        % Separation rate (BLS = .04 - quarterly -)
% l0 = p.l0;      % UE contact rate (Jarosch's = .06 - monthly -)
% l1 = p.l1;      % EE contact rate (Jarosch's = .08 - monthly -)
% 
% % Firm's productivity (log Normal)
% mu = p.mu;
% s2 = p.s2;
% n0 = p.n0;
% 
% % Approximate log Normal with qnwnorm
% [x0, w0] = qnwnorm(n0, mu, s2);
% y0 = exp(x0);
% 
% G = logncdf(y0, mu, s2);
% Gc = 1 - G;
% 
% % Parameters to initialize VF iteration
% tol = 1e-7;    % tolerance
% dif = 1;         % dif between S's
% iter = 0;       % iteration counter
% 
% %----------------------------------------------------------------------
% %------------------         VALUE FUNCTION ITERATION        -------------------
% %----------------------------------------------------------------------
% 
% S = (y0 - z)  / (1 - b * d);
% div = 1 - b * (1 - d) *(1 - a * l1 * Gc);
% 
% while dif > tol
%     S = max(S, 0);
%     Su = repmat(S, 1, n0);
%     Se = tril(Su);
%     TS = ((y0 - z) + b * a * (1 - d) * (w0' * (l1 * Se - l0 * Su))') ./ div;
%     TS = max(TS, 0);
% 
%     dif = norm(abs(TS - S));
%     iter = iter + 1;
%     S = TS;
% end
% 
% % Value of unemployment
% U = (z + b * l0 * a * w0(S >= 0)' * S(S >= 0)) / (1 - b);
% 
% % Value of being employed at y (row) with benchmark y_ (col)
% W = tril(U + repmat(S', n0, 1) + a * (repmat(S, 1, n0) - repmat(S', n0, 1)));
% 
% % Choice sets
% M1 = zeros(n0, n0, n0);
% M2 = zeros(n0, n0, n0);
% 
% for i = 1:n0          % y:   current
%     for j = 1:n0      % y_: benchmark
%         for k = 1:n0  % y': offer
%             M1(i, j, k) = (k > i);                  % movers
%             M2(i, j, k) = (k <= i) * (k > j);  % stayers + wage increase
%         end
%     end
% end
% 
% % Wages
% cv = zeros(n0);  % continuation value
% 
% for i = 1:n0         % current 
%     for j = 1:n0     % benchmark
%         m1 = squeeze(M1(i, j, :));
%         m2 = squeeze(M2(i, j, :));
%         cv(i, j) = l1 * w0' * (m1 .*  W(:, i) + m2 .* W(i, :)') + (1 - l1 * w0' * (m1 + m2)) * W(i, j);
%     end
% end
% 
% % Wage of worker at y (row) with benchmark y_ (col)
% wage = tril(W - b * d * U - b * (1 - d) * cv);
% 
% p.w0 = w0;
% p.y0 = y0;
% 
% sol.U = U;
% sol.W = W;
% sol.S = S;
% sol.wage = wage;
% 
% end




%----------------------------------------------------------------------
%----------------------------------------------------------------------




% OLD PLOTS CODE

% set(groot, 'DefaultAxesLineWidth', 1.5);
% set(groot, 'DefaultLineLineWidth', 2);
% set(groot, 'DefaultTextInterpreter', 'latex');
% set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
% set(groot, 'DefaultLegendInterpreter', 'latex');
% set(groot, 'DefaultAxesFontSize', 20);
% 
% set(0,'DefaultLegendAutoUpdate','off')
% 
% cmap = [249,65,68; 
%             243,114,44;
%             248,150,30;
%             249,132,74;
%             249,199,79;
%             144,190,109;
%             67,170,139; 
%             77,144,142; 
%             87,117,144; 
%             39,125,161] / 255;
% 
% pick = [1 2 4 5 6 9 10];
% 
% fullfig
% 
% colororder(cmap)
% 
% t = (freq/4):(freq/4):25;
% min_w = sol.nU;
% max_w = max(sol.wage(p.w > .01, sol.nU));
% 
% subplot(2, 2, 1)
% plot(t, stats{1, 5}(pick, :)')
% title('Wages -- Percentiles at Mean')
% xlabel('Time (years)')
% ylabel('Wages')
% legend(strcat('$p=$',num2str(pick')),'location', 'NW')
% 
% subplot(2, 2, 2)
% plot(t, stats{2, 5}(pick, :)')
% title('Firm productivity -- Percentiles at Mean')
% xlabel('Time (years)')
% ylabel('Productivity')
% legend(strcat('$p=$',num2str(pick')),'location', 'NW')
% 
% subplot(2, 2, 3)
% plot(t, stats{3, 5}(pick, :)')
% title('Unemployment rate -- Percentiles with Mean')
% xlabel('Time (years)')
% ylabel('U rate')
% legend(strcat('$p=$',num2str(pick')),'location', 'NW')
% 
% subplot(2, 2, 4)
% plot(t, stats{4, 5}(pick, :)')
% title('\# firms -- Percentiles at Mean')
% xlabel('Time (years)')
% ylabel('\# firms')
% legend(strcat('$p=$',num2str(pick')),'location', 'NW')
% 
% sgt = sgtitle(strcat('Results for $\alpha=', num2str(p.a), '$'));
% sgt.FontSize = 26;