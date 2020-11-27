close all;
load('figure_format.mat');
%%
if (   ~exist('N','var')  || ~exist('N_vec0','var') || ~exist('N_basis','var') ...
    || ~exist('CP_num','var') || ~exist('CQ_num','var') || ~exist('Mdiff_num','var') ...
    || ~exist('CP_sym','var') || ~exist('CQ_sym','var') || ~exist('Mdiff_sym','var') )
    T_A = tic;
    fprintf('Loading M0 matrix. ');
    load('MatrixData.mat','N','N_vec0','N_basis','CP_num','CQ_num','Mdiff_num','CP_sym','CQ_sym','Mdiff_sym','M0');
    fprintf('Time: %s.\n\n',datestr(seconds(toc(T_A)),'MM:SS.FFF'));   
end
%%
ind = indP(0:N_basis);
per_err_CP = abs(CP_num(ind) - CP_sym)./abs(CP_num(ind)) * 100;
per_err_CQ = abs(CQ_num(ind) - CQ_sym)./abs(CQ_num(ind)) * 100;
per_err_Mdiff = abs(Mdiff_num(ind) - Mdiff_sym)./abs(Mdiff_num(ind))* 100;

fprintf('Maximum and median percent error of symbolic evaluation vs. numeric calculation for n = 0 : %d:\n', N_basis); 
fprintf('\tCP: %1.3g%%, %1.3g%%.\n', max(per_err_CP), median(per_err_CP));
fprintf('\tCQ: %1.3g%%, %1.3g%%.\n', max(per_err_CQ), median(per_err_CQ));
fprintf('\tM'': %1.3g%%, %1.3g%%.\n', max(per_err_Mdiff), median(per_err_Mdiff));

M_diag = diag(M0);

%%
h_f = figure;
h_a = axes;

plot(h_a, N_vec0, CP_num, 'Marker','o','MarkerSize',2, 'MarkerFaceColor','k', 'MarkerEdgeColor','none','LineStyle','none');

format_axis.xlim = [min(N_vec0),max(N_vec0)];
max_y = ceil(max(abs(CP_num)/10))*10;
format_axis.ylim = [-1,1] * max_y;

set([h_a.XLabel,h_a.YLabel], format_axis_label);
xlabel(h_a, '$$n$$');
ylabel(h_a, '$$\mathrm{c}_{2n}^{\mathrm{MP}}$$');
set(h_a, format_axis);
set(h_a,'XTick',[0,N/8:N/8:N])

set(findobj(h_f, 'Type','line'), format_line);
set(h_f,format_figure,'Position',[50,50,1200,800]);

hold on;
text(h_a, N/2,  abs(CP_num(indP(N/2)))*0.93, 'Even $$n$$');
text(h_a, N/2, -abs(CP_num(indP(N/2)))*0.93, 'Odd $$n$$');
text(h_a, N/2, 0, '$$\mathrm{c}_{2n}^{\mathrm/8{MP}}\rightarrow (-1)^{n}\sqrt{\pi n}$$');
% plot(h_a, N_vec0(2:end), sqrt(pi*N_vec0(2:end)),'r--','LineWidth',2);
% plot(h_a, N_vec0(2:end),-sqrt(pi*N_vec0(2:end)),'r--','LineWidth',2);
set(findobj(h_f, 'Type','text'), format_text);

figure_name = 'CP_MP';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(:,51:1150,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%%
h_f = figure;
h_a = axes;

plot(h_a, N_vec0, Mdiff_num, 'Marker','o','MarkerSize',2, 'MarkerFaceColor','k', 'MarkerEdgeColor','none','LineStyle','none');

format_axis.xlim = [min(N_vec0),max(N_vec0)];
format_axis.ylim = [floor(min(Mdiff_num)), ceil(max(Mdiff_num))];
set([h_a.XLabel,h_a.YLabel], format_axis_label);
xlabel(h_a, '$$n$$');
ylabel(h_a, '$$M''_{2n}(0)$$');
set(h_a, format_axis);
set(h_a,'XTick',[0,N/8:N/8:N]);

set(findobj(h_f, 'Type','line'), format_line);
set(h_f,format_figure,'Position',[50,50,1200,800]);

hold on;
text(N/2, -N, '$$M''_{2n}(0)\rightarrow\ -2n$$','VerticalAlignment','Bottom');
set(findobj(h_f, 'Type','text'), format_text);


figure_name = 'MM_diff';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(:,51:1150,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%%
h_f = figure;
h_a = axes;

plot(h_a, N_vec0, M_diag,'Marker','o','MarkerSize',5,'MarkerFaceColor','k','Color','k', 'MarkerEdgeColor','none');

format_axis.xlim = [min(N_vec0),max(N_vec0)/32];
format_axis.ylim = [pi/8 * 0.99, max(M_diag)];
set([h_a.XLabel,h_a.YLabel], format_axis_label);
xlabel(h_a,'$$n$$');
ylabel(h_a,'$$m_{n,n}$$');
set(h_a, format_axis);
set(h_a,'XTick',[0,(1/8:1/8:1)*format_axis.xlim(2)]);

set(findobj(h_f, 'Type','line'), format_line);
set(h_f,format_figure,'Position',[50,50,1200,800]);

hold on;
plot(N_vec0,pi/8*ones(size(N_vec0)),'k--','LineWidth',0.5);

h_a.YTick = [pi/8, h_a.YTick];
h_a.YTickLabel{1} = '$$\pi/8$$';    
text(h_a,format_axis.xlim(2)/2, pi/8 * 1.01, ...
        '$$m_{n,n}\rightarrow \pi/8$$','VerticalAlignment','Bottom');
set(findobj(h_f, 'Type','text'), format_text);


figure_name = 'M_diag';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(:,51:1150,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');
