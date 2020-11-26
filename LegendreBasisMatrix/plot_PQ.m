close all;
load('figure_format.mat');
%%
if ~exist('N_basis','var') || ~exist('x','var') || ~exist('P','var') || ~exist('Q','var')
    T_Leg = tic;
    fprintf('Loading Legendre Functions. ');
    load('LegendreFunctions.mat','N_basis','x','P','Q');
    fprintf('Time: %s.\n\n',datestr(seconds(toc(T_Leg)),'MM:SS.FFF'));
end

N_plot = 5;

delta_x = 0.001;
xP_num = (-1: delta_x: 1)';
xQ_num = (-1 + delta_x: delta_x: 1-delta_x)';
P_num = zeros(length(xP_num), N_plot*2);
Q_num = zeros(length(xQ_num), N_plot*2);

for n = 0 : N_plot*2 -1
    P_num(:,indP(n)) = subs( P{indP(n)}, x, xP_num);
    Q_num(:,indP(n)) = subs( Q{indP(n)}, x, xQ_num);
end

%% PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
h_f = figure;
h_sp = subplot(20,2,[1,2]);
title(h_sp, 'Legendre functions of the first kind');
set(h_sp.Title, format_title);
set(h_sp,format_blank_axis);

tP_locX = [0,-0.55, 0, 0.42, 0.63, 0.72, 0.82, 0.85, 0.91, 0.95];
tP_locY = -0.6;
%%
h_sp = subplot(20,2,[3,39]);
hold on;

plot(h_sp, xP_num, P_num(:,indP(0:2:2*(N_plot-1))),'k');

text(h_sp, tP_locX(indP(0)), 0.95, '$P_0$');

for ii = 2 : 2 : 2* (N_plot-1)
    text(h_sp, tP_locX(indP(ii)), tP_locY, sprintf('$P_%d$',ii));
    [min_P, ind_P] = min(P_num(:,indP(ii)));
    length_arrow = sqrt((tP_locX(indP(ii))-(-xP_num(ind_P)))^2+(min_P - tP_locY)^2);
    plot_arrow(tP_locX(indP(ii)),tP_locY,-xP_num(ind_P),min_P,'headwidth', 0.04 / (length_arrow),'headheight', 0.04 / (length_arrow));
end


format_axis.Xlim = [-1.0,1.0];
format_axis.Ylim = [-1.05,1.05];
ylabel(h_sp, '$$P_{2n}(x)$$');
xlabel(h_sp, '$$x$$');

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel], format_axis_label);

%%
h_sp = subplot(20,2,[4,40]);
hold on;

plot(h_sp, xP_num, P_num(:,indP(1:2:2*N_plot-1)),'k');


for ii = 1 : 2 : 2 * N_plot -1 
    text(h_sp, tP_locX(indP(ii)), tP_locY, sprintf('$P_%d$',ii));
    if ii == 1 
        continue
    end
    [max_P, ind_P] = max(P_num(1:(length(xP_num)+1)/2,indP(ii)));
    length_arrow = sqrt((tP_locX(indP(ii))-(-xP_num(ind_P)))^2+(-max_P - tP_locY)^2);
    plot_arrow(tP_locX(indP(ii)),tP_locY,-xP_num(ind_P),-max_P,'headwidth', 0.04 / (length_arrow),'headheight', 0.04 / (length_arrow));
end


format_axis.Xlim = [-1.0,1.0];
format_axis.Ylim = [-1.05,1.05];
ylabel(h_sp, '$$P_{2n+1}(x)$$');
xlabel(h_sp, '$$x$$');

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel], format_axis_label);

%%
set(h_f,format_figure);
set(findobj(h_f, 'Type','line'), format_line);
set(findobj(h_f, 'Type','text'), format_text, 'VerticalAlignment','Top', 'HorizontalAlignment','Center');
%%
figure_name = 'Legendre_P';
saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(    imind,cm,fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%% QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
h_f = figure;
h_sp = subplot(20,2,[1,2]);
title(h_sp, 'Legendre functions of the second kind');
set(h_sp.Title, format_title);
set(h_sp,format_blank_axis);

tQ_locX = [-0.8, 0, 0.65, 0.65, 0.75, 0.75, 0.85, 0.85, 0.95, 0.95];
tQ_locY = -1.2;
%%
h_sp = subplot(20,2,[3,39]);
hold on;

plot(h_sp, xQ_num, Q_num(:,indP(0:2:2*(N_plot-1))),'k');

for ii = 0 : 2 : 2* (N_plot-1)
    text(h_sp, tQ_locX(indP(ii)), tQ_locY, sprintf('$Q_%d$',ii));
    if ii == 0 
        continue
    end
    [max_Q, ind_Q] = max(Q_num(1:(length(xQ_num)+1)/2,indP(ii)));
    length_arrow = sqrt((tQ_locX(indP(ii))-(-xQ_num(ind_Q)))^2+(-max_Q - tQ_locY)^2);
    plot_arrow(tQ_locX(indP(ii)),tQ_locY,-xQ_num(ind_Q),-max_Q,'headwidth', 0.04 / (length_arrow),'headheight', 0.08 / (length_arrow));
end

format_axis.Xlim = [-1,1];
format_axis.Ylim = [-2,2];
ylabel(h_sp, '$$Q_{2n}(x)$$');
xlabel(h_sp, '$$x$$');

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel], format_axis_label);

%%
h_sp = subplot(20,2,[4,40]);
hold on;

plot(h_sp, xQ_num, Q_num(:,indP(1:2:2*N_plot-1)),'k');

for ii = 1 : 2 : 2 * N_plot -1 
    text(h_sp, tQ_locX(indP(ii)), tQ_locY, sprintf('$Q_%d$',ii));
    [min_Q, ind_Q] = min(Q_num(1:(length(xQ_num)+1)/2,indP(ii)));
    length_arrow = sqrt((tQ_locX(indP(ii))-(-xQ_num(ind_Q)))^2+(min_Q - tQ_locY)^2);
    plot_arrow(tQ_locX(indP(ii)),tQ_locY,-xQ_num(ind_Q),min_Q,'headwidth', 0.04 / (length_arrow),'headheight', 0.08 / (length_arrow));
end

format_axis.Xlim = [-1,1];
format_axis.Ylim = [-2,2];
ylabel(h_sp, '$$Q_{2n+1}(x)$$');
xlabel(h_sp, '$$x$$');

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel], format_axis_label);

%%
set(h_f,format_figure);
set(findobj(h_f, 'Type','line'), format_line);
set(findobj(h_f, 'Type','text'), format_text, 'VerticalAlignment','Top', 'HorizontalAlignment','Center');

%%
figure_name = 'Legendre_Q';
saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(    imind,cm,fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

