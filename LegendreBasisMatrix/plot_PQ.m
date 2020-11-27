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

ind0 = (xP_num ==0);
offsetX = [0, -0.1, 0,  0, 0, -0.025, 0, -0.025, 0, +0.025];
offsetY = [-0.05, 0, -0.05,  -0.05, 0.05, -0.05, 0.05, -0.05, -0.05, -0.05];

%
h_sp = subplot(20,2,[3,39]);
hold on;

plot(h_sp, xP_num, P_num(:,indP(0:2:2*(N_plot-1))),'k');

    
for ii = 0 : 2 : 2* (N_plot-1)
    text(h_sp, 0, P_num( ind0, indP(ii)) + offsetY(indP(ii)), sprintf('$P_%d$',ii));
end

format_axis.Xlim = [-1.0,1.0];
format_axis.Ylim = [-1.05,1.05];
ylabel(h_sp, '$$P_{2n}(x)$$');
xlabel(h_sp, '$$x$$');

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel], format_axis_label);

%
h_sp = subplot(20,2,[4,40]);
hold on;

plot(h_sp, xP_num, P_num(:,indP(1:2:2*N_plot-1)),'k');

for ii = 1 : 2 : 2 * N_plot -1 
    if ii == 1
        [exP,indx] = min(P_num(xP_num<=0,indP(ii)));
    else
        [exP,indx] = max(P_num(xP_num<=0,indP(ii)));
    end
    text(h_sp, -xP_num(indx) + offsetX(indP(ii)), -exP + offsetY(indP(ii)), sprintf('$P_%d$',ii));
end


format_axis.Xlim = [-1.0,1.0];
format_axis.Ylim = [-1.05,1.05];
ylabel(h_sp, '$$P_{2n+1}(x)$$');
xlabel(h_sp, '$$x$$');

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel], format_axis_label);

%
set(h_f,format_figure);
set(findobj(h_f, 'Type','line'), format_line);
set(findobj(h_f, 'Type','text'), format_text, 'VerticalAlignment','Middle', 'HorizontalAlignment','Center');
%%
figure_name = 'Legendre_P';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(:,101:1400,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%% QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
h_f = figure;
h_sp = subplot(20,2,[1,2]);
title(h_sp, 'Legendre functions of the second kind');
set(h_sp.Title, format_title);
set(h_sp,format_blank_axis);

tQ_locX = [-0.8, 0, 0.65, 0.65, 0.75, 0.75, 0.85, 0.85, 0.95, 0.95];
tQ_locY = -1.2;
ind0 = (xQ_num ==0);
indxQ = xQ_num<=0 & xQ_num >-0.9;
offxQ = find(indxQ, 1);
            
offsetX = [0.05, 0,  0, 0, 0, 0, +0.05, 0, +0.05,0];
offsetY = [-0.05, -0.1, +0.1, +0.1,  +0.1, -0.1, +0.1, -0.1, +0.1, +0.1];
%
h_sp = subplot(20,2,[3,39]);
hold on;

plot(h_sp, xQ_num, Q_num(:,indP(0:2:2*(N_plot-1))),'k');

for ii = 0 : 2 : 2* (N_plot-1)
    switch ii
        case 0
            [exQ,indx] = min(Q_num(indxQ,indP(ii)));
            text(h_sp, xQ_num(indx+offxQ) + offsetX(indP(ii)), exQ + offsetY(indP(ii)), sprintf('$Q_%d$',ii));
        case {2}
            [exQ,indx] = max(Q_num(indxQ,indP(ii)));
            text(h_sp, xQ_num(indx+offxQ) + offsetX(indP(ii)), exQ + offsetY(indP(ii)), sprintf('$Q_%d$',ii));
        case {4,6,8}
            [exQ,indx] = min(Q_num(indxQ,indP(ii)));
            text(h_sp, -xQ_num(indx+offxQ) + offsetX(indP(ii)), -exQ + offsetY(indP(ii)), sprintf('$Q_%d$',ii));
        
    end
    
end

format_axis.Xlim = [-1,1];
format_axis.Ylim = [-2,2];
ylabel(h_sp, '$$Q_{2n}(x)$$');
xlabel(h_sp, '$$x$$');

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel], format_axis_label);

%
h_sp = subplot(20,2,[4,40]);
hold on;

plot(h_sp, xQ_num, Q_num(:,indP(1:2:2*N_plot-1)),'k');

for ii = 1 : 2 : 2 * N_plot -1 
    text(h_sp, xQ_num(ind0), Q_num( ind0, indP(ii)) + offsetY(indP(ii)), sprintf('$Q_%d$',ii));
end

format_axis.Xlim = [-1,1];
format_axis.Ylim = [-2,2];
ylabel(h_sp, '$$Q_{2n+1}(x)$$');
xlabel(h_sp, '$$x$$');

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel], format_axis_label);

%
set(h_f,format_figure);
set(findobj(h_f, 'Type','line'), format_line);
set(findobj(h_f, 'Type','text'), format_text, 'VerticalAlignment','Middle', 'HorizontalAlignment','Center');

%%
figure_name = 'Legendre_Q';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(:,101:1400,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

