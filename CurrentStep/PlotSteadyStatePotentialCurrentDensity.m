close all;
addpath(fullfile('..','LegendreBasisMatrix'))
load('figure_format.mat');
%%

if  ~exist('G','var')  || ~exist('b0SS','var') || ~exist('VSS','var') || ~exist('G','var')
    load('I_SS.mat','b0SS','VSS','G');
end
load('BasisFunctions_U0J0_num.mat','N_basis','N_plot','r_num','J0_num','U0_num');

G_plot = 10.^(-2:0.5:2); 
[~,ind_G_plot] = intersect(G, G_plot);
G_plot_lines = kron(G_plot,ones(size(r_num))); 
r_plot_lines = repmat(r_num,size(G_plot));


PhiSS = zeros(length(r_num), length(G));
JSS = zeros(length(r_num)-1, length(G));

format_axis.XLim = 10.^[-2,2];
format_axis.XTick = 10.^(-2:2);
format_axis.XScale = 'log';
format_axis.XDir = 'reverse';
format_axis.YLim = [0,1];
format_axis.YTick = 0:0.2:1;
   
for ii = 1 : length(G)
    PhiSS(:,ii) = U0_num * b0SS(indP(0:N_basis), ii);
    JSS(:,ii) = J0_num * b0SS(indP(0:N_basis) ,ii) / b0SS(indP(0),ii) * pi/4;
end

save('I_SS.mat', 'G_plot', 'r_num', 'PhiSS', 'JSS', '-append');
%%
h_f = figure;
h_sp = subplot(20,2,[1,2]);
title(sprintf('Steady state of current step response'));
set(h_sp.Title, format_title);
set(h_sp,format_blank_axis);

%
h_sp = subplot(20,2,[3,39]);
surf(G, r_num, PhiSS,'LineStyle' ,'none','FaceAlpha',0.5);
hold on;box on;
plot3(G_plot_lines, r_plot_lines, PhiSS(:,ind_G_plot));


format_axis.ZLim = [0,1.5];
format_axis.ZTick = (0:0.25:1.5);

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel,h_sp.ZLabel], format_axis_label);
caxis(h_sp, [0,max(PhiSS(:))])
view(h_sp, [75, 25])
xlabel({'Conductance  $G$'});
ylabel('Radial position $r/r_{0}$');
zlabel({'Potential $\varphi^{\mathrm{SS}}_{0}(r)/V_{0}$'});
    
%
h_sp = subplot(20,2,[4,40]);
surf(G, r_num(1:end-1), JSS,'LineStyle' ,'none','FaceAlpha',0.5);     % normalization factor: 4*B_0/pi=4/pi
hold on;box on;
plot3(G_plot_lines(1:end-1,:), r_plot_lines(1:end-1,:), JSS(:,ind_G_plot));

format_axis.ZLim = [0,4];
format_axis.ZTick = (0:0.5:4);

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel,h_sp.ZLabel], format_axis_label);
caxis(h_sp, [0.5,4])
view(h_sp, [75,25])

xlabel({'Conductance  $G$'});ylabel('Radial position $r/r_{0}$');
zlabel({'Current density  $J^{\mathrm{SS}}_{0}(r)/\overline{J_{0}}$'});


set(h_f,format_figure,'Position',[0,0,1800,850]);
set(findobj(h_f, 'Type','line'), format_line, 'LineWidth', 1.5,'Color','k');
set(findobj(h_f, 'Type','text'), format_text);

%%
figure_name = 'U0J0_SS';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(:,51:1750,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%%
rmpath(fullfile('..','LegendreBasisMatrix'))
