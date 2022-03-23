%----------------------------------------------------------------------
%----------------------------         PLOTS        ---------------------------
%----------------------------------------------------------------------

set(groot, 'DefaultAxesLineWidth', 1.5);
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesFontSize', 20);

set(0,'DefaultLegendAutoUpdate','off')

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

hexcol = {'#00429d', '#325da9', '#4e78b5', '#6694c1', '#80b1cc', '#9dced6', ...
                '#c0eade', '#ffffe0', '#ffdac4', '#ffb3a7', '#fb8a8c', '#eb6574', ...
                '#d5405e', '#b81b4a', '#93003a'};

cmap = hex2rgb(hexcol);
pick = 1:10;
t = 19 + ((freq/4):(freq/4):36);
pos = [.03 .35 .05 .25];

%----------------------------------------------------------------------
%----------------------------------------------------------------------


% 1. WAGE SCHEDULE
%----------------

fullfig

w_sup_min = sol.nU;
w_sup_max = find((p.w > .01) .* (1:p.n)' > 0, 1, 'last');
w_plot = w_sup_min:2:w_sup_max;
colororder(cmap)

plot(p.y(p.w > p.lb), sol.wage(p.w > p.lb, w_plot))

title(strcat('Wage schedule with different benchmark firms -- $\alpha=', num2str(p.a), '$'));
legend([strcat('$y_{-1}=$',num2str(w_plot')) ''],'location', 'best')
xlabel('Current firm productivity')
ylabel('Wage')

yyaxis right
area(p.y(p.w > p.lb), p.w(p.w > p.lb), 'FaceColor', 'k', 'LineStyle','none', 'FaceAlpha', .1);

saveas(gcf,strcat('wages_a', num2str(100*p.a), '.png'))


%----------------------------------------------------------------------
%----------------------------------------------------------------------


% 2. WAGE DECOMPOSITION (CONTRIBUTIONS)
%------------------------------------

X(:, 1) = sol.W(:, sol.nU);
X(:, 2) = sol.p3(:, sol.nU);
X(:, 3) = sol.p1(:, sol.nU);
X(:, 4) = sol.p2(:, sol.nU);

fullfig

colororder(cmap([1 5], :))

subplot(2, 1, 1)
plot(p.y(sol.nU:end), X(sol.nU:end, [1 2]))

legend('Entire Value', 'Current','location', 'NW', 'FontSize', 24)
title(strcat('Entry Wage Decomposition -- Value Function and NPV -- $\alpha=', num2str(p.a), '$'));
xlabel('Firm Productivity')
ylabel('Value')


subplot(2, 1, 2)
plot(p.y(sol.nU:end), X(sol.nU:end, [3 4]))

legend('Move', 'Bargain','location', 'NW', 'FontSize', 24)
title('Entry Wage Decomposition -- Move and Bargain')
xlabel('Firm Productivity')
ylabel('Value')

saveas(gcf,strcat('wage_contrib_a', num2str(100*p.a), '.png'))


% %----------------------------------------------------------------------
% %----------------------------------------------------------------------
% 
% 
% % 3. WAGES AND PRODUCTIVITY (MEDIAN)
% %---------------------------------
% 
% fullfig
% 
% colororder(cmap)
% 
% min_w = sol.nU;
% max_w = max(sol.wage(p.w > .01, sol.nU));
% 
% subplot(2, 2, 1)
% plot(t, stats{1, 2}(pick, :)')
% hold on
% patch([t fliplr(t)], [repmat(min_w, 1, length(t)) repmat(max_w, 1, length(t))], ...
%             'k', 'LineStyle','none', 'FaceAlpha', .1)
% hold off
% title('Wages -- BOC percentiles')
% xlabel('Time (years)')
% ylabel('Wage')
% legend(strcat('$p=$',num2str(pick')),'location', 'NW')
% 
% subplot(2, 2, 2)
% plot(t, stats{1, 3}(pick, :)')
% hold on
% patch([t fliplr(t)], [repmat(min_w, 1, length(t)) repmat(max_w, 1, length(t))], ...
%             'k', 'LineStyle','none', 'FaceAlpha', .1)
% hold off
% title('Wages -- EOC percentiles')
% xlabel('Time (years)')
% ylabel('Wage')
% legend(strcat('$p=$',num2str(pick')),'location', 'NW')
% 
% subplot(2, 2, 3)
% plot(t, stats{2, 2}(pick, :)')
% title('Firm productivity -- BOC percentiles')
% xlabel('Time (years)')
% ylabel('Productivity')
% legend(strcat('$p=$',num2str(pick')),'location', 'NW')
% 
% subplot(2, 2, 4)
% plot(t, stats{2, 3}(pick, :)')
% title('Firm productivity -- EOC percentiles')
% xlabel('Time (years)')
% ylabel('Productivity')
% legend(strcat('$p=$',num2str(pick')),'location', 'NW')
% 
% sgt = sgtitle(strcat('Median Wages and Productivity -- $\alpha=', num2str(p.a), '$'));
% sgt.FontSize = 26;
% 
% saveas(gcf,strcat('wages_prod_med_a', num2str(100*p.a), '.png'))
% 
% 
% %----------------------------------------------------------------------
% %----------------------------------------------------------------------
% 
% 
% % 4. UNEMPLOYMENT AND # FIRMS (MEDIAN)
% %----------------------------------
% 
% fullfig
% 
% colororder(cmap)
% 
% subplot(2, 2, 1)
% plot(t, stats{3, 2}(pick, :)')
% title('Unemployment rate -- BOC percentiles')
% xlabel('Time (years)')
% ylabel('U rate')
% legend(strcat('$p=$',num2str(pick')),'location', 'NW')
% 
% subplot(2, 2, 2)
% plot(t, stats{3, 3}(pick, :)')
% title('Unemployment rate -- EOC percentiles')
% xlabel('Time (years)')
% ylabel('U rate')
% legend(strcat('$p=$',num2str(pick')),'location', 'NW')
% 
% subplot(2, 2, 3)
% plot(t, stats{4, 2}(pick, :)')
% title('\# firms -- BOC percentiles')
% xlabel('Time (years)')
% ylabel('\# firms')
% legend(strcat('$p=$',num2str(pick')),'location', 'NW')
% 
% subplot(2, 2, 4)
% plot(t, stats{4, 3}(pick, :)')
% title('\# firms -- EOC percentiles')
% xlabel('Time (years)')
% ylabel('\# firms')
% legend(strcat('$p=$',num2str(pick')),'location', 'NW')
% 
% sgt = sgtitle(strcat('Unemployment and Median \# firms -- $\alpha=', num2str(p.a), '$'));
% sgt.FontSize = 26;
% 
% saveas(gcf,strcat('unemp_firms_med_a', num2str(100*p.a), '.png'))
% 
% 
% %----------------------------------------------------------------------
% %----------------------------------------------------------------------


% 3. WAGES AND PRODUCTIVITY (MEAN)
%-------------------------------

fullfig

colororder(cmap)

min_w = sol.nU;
max_w = max(sol.wage(p.w > .01, sol.nU));

subplot(2, 3, 1)
plot(t, stats{1, 4}(pick, :)')
% hold on
% patch([t fliplr(t)], [repmat(min_w, 1, length(t)) repmat(max_w, 1, length(t))], ...
%             'k', 'LineStyle','none', 'FaceAlpha', .1)
% hold off
title('Wages -- BOC percentiles')
xlabel('Time (years)')
ylabel('Wage')
legend(strcat('$p=$',num2str(pick')),'location', 'NW')

subplot(2, 3, 2)
plot(t, stats{1, 5}(pick, :)')
% hold on
% patch([t fliplr(t)], [repmat(min_w, 1, length(t)) repmat(max_w, 1, length(t))], ...
%             'k', 'LineStyle','none', 'FaceAlpha', .1)
% hold off
title('Wages -- EOC percentiles')
xlabel('Time (years)')
ylabel('Wage')

subplot(2, 3, 3)
plot(t, stats{1, 6}(pick, :)')
% hold on
% patch([t fliplr(t)], [repmat(min_w, 1, length(t)) repmat(max_w, 1, length(t))], ...
%             'k', 'LineStyle','none', 'FaceAlpha', .1)
% hold off
title('Wages -- Mean LW percentiles')
xlabel('Time (years)')
ylabel('Wage')

%----------------------------------------------

subplot(2, 3, 4)
plot(t, stats{2, 4}(pick, :)')
title('Firm productivity -- BOC percentiles')
xlabel('Time (years)')
ylabel('Productivity')

subplot(2, 3, 5)
plot(t, stats{2, 5}(pick, :)')
title('Firm productivity -- EOC percentiles')
xlabel('Time (years)')
ylabel('Productivity')

subplot(2, 3, 6)
plot(t, stats{2, 6}(pick, :)')
title('Firm productivity -- Mean LW percentiles')
xlabel('Time (years)')
ylabel('Productivity')

legend(strcat('$p=$', num2str(pick')), 'FontSize', 24, 'Position', pos);
sgt = sgtitle(strcat('Mean Wages and Productivity -- $\alpha=', num2str(p.a), '$'));

sgt.FontSize = 26;

saveas(gcf,strcat('wages_prod_mean_a', num2str(100*p.a), '.png'))


%----------------------------------------------------------------------
%----------------------------------------------------------------------


% 4. UNEMPLOYMENT AND # FIRMS (MEAN)
%---------------------------------

fullfig

colororder(cmap)

subplot(2, 3, 1)
plot(t, stats{3, 4}(pick, :)')
title('Unemployment rate -- BOC percentiles')
xlabel('Time (years)')
ylabel('U rate')

subplot(2, 3, 2)
plot(t, stats{3, 5}(pick, :)')
title('Unemployment rate -- EOC percentiles')
xlabel('Time (years)')
ylabel('U rate')

subplot(2, 3, 3)
plot(t, stats{3, 6}(pick, :)')
title('Unemployment rate -- Mean LW percentiles')
xlabel('Time (years)')
ylabel('U rate')

%----------------------------------------------

subplot(2, 3, 4)
plot(t, stats{4, 4}(pick, :)')
title('Mean \# firms -- BOC percentiles')
xlabel('Time (years)')
ylabel('\# firms')

subplot(2, 3, 5)
plot(t, stats{4, 5}(pick, :)')
title('Mean \# firms -- EOC percentiles')
xlabel('Time (years)')
ylabel('\# firms')

subplot(2, 3, 6)
plot(t, stats{4, 6}(pick, :)')
title('Mean \# firms -- Mean LW percentiles')
xlabel('Time (years)')
ylabel('\# firms')

legend(strcat('$p=$', num2str(pick')), 'FontSize', 24, 'Position', pos);
sgt = sgtitle(strcat('Unemployment and Mean \# firms -- $\alpha=', num2str(p.a), '$'));
sgt.FontSize = 26;

saveas(gcf,strcat('unemp_firms_mean_a', num2str(100*p.a), '.png'))

%close all