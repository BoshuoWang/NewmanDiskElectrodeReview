%%
close all;
addpath(fullfile('..','LegendreBasisMatrix'));
addpath(fullfile('..','CurrentStep'));

load('figure_format.mat');
if ~exist('N','var')|| ~exist('M0','var')
    T_A = tic;
    fprintf('Loading matrix. ');
    load('MatrixData.mat', 'N', 'M0');
end
%%
load('BasisFunctions_U0J0_num.mat','N_basis','r_num','J0_num','U0_num');

dr = min(diff(r_num));      % 1e-4
r_num_dense = 0 : dr : r_num(end);

lp_filt = designfilt(   'lowpassfir', ...
                        'DesignMethod', 'equiripple', 'SampleRate', 1/dr,...
                        'PassbandFrequency', 20, 'PassbandRipple',1,...
                        'StopbandFrequency', 100, 'StopbandAttenuation',80);

% Direct-Form FIR, Equiripple
% Fs = 100000 cycle per unit length
% Fpass = 20; Apass = 1 ; 
% Fstop = 100;Astop = 80;

%%
load('I_SS.mat','b0SS','VSS','G');
load('V_TZ.mat','Lambda','b0TZ');

K = VSS.^(-1);

N_sum = N_basis;
G_plot = [0,logspace(-2,2,9)];
[~, ind_G_plot] = intersect(G, G_plot);
G_plot_lines = kron(G_plot,ones(size(r_num)));
r_plot_lines = repmat(r_num,size(G_plot));

C = zeros(N_basis+1, length(G));
C_alt = zeros(N_basis+1, length(G));

for jj = 0 : N_basis
    Denom = b0TZ{end}(:,indP(jj))' * M0 * b0TZ{end}(:,indP(jj));
    for ii = 1 : length(G)
        C(indP(jj), ii) = - Lambda{end}(indP(jj)) / ( (Lambda{end}(indP(jj)) + G(ii) ) * 2 * Denom );
        C_alt(indP(jj),ii) = ( (b0SS(:,ii)*K(ii))' * M0 * b0TZ{end}(:,indP(jj)) - 1/2 ) / ( Denom );
    end
end
C_err = C_alt ./ C - 1;

fprintf('Relative errors of C using different methods: Max: %2.2g. Min: %2.2g. Mean: %2.2g. StdDev: %2.2g.\n', ...
            nanmax(C_err(:)), nanmin(C_err(:)), nanmean(C_err(:)),  nanstd(C_err(:)));


% SS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PhiSS = zeros(length(r_num), length(G));
JSS = zeros(length(r_num)-1, length(G));

for ii = 1 : length(G)
    PhiSS(:,ii) = U0_num * b0SS(indP(0:N_basis), ii) * K(ii);
    JSS(:,ii)   = J0_num * b0SS(indP(0:N_basis), ii) * K(ii) * pi/4;
end

save('V_SS.mat', 'G', 'G_plot', 'K', 'r_num', 'PhiSS', 'JSS');

% TZi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PhiTZi = zeros(length(r_num), N_sum+1);
JTZi = zeros(length(r_num)-1, N_sum+1);

for jj = 0 : N_sum
    PhiTZi(:,indP(jj)) = U0_num(:,indP(0 : N_sum)) * b0TZ{end}(indP(0 : N_sum), indP(jj) );
    JTZi(:, indP(jj))  = J0_num(:,indP(0 : N_sum)) * b0TZ{end}(indP(0 : N_sum), indP(jj) ) * pi/4;
end

% TZ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PhiTZ = zeros(length(r_num), length(G));
JTZ = zeros(length(r_num)-1, length(G));
JTZ_filt = zeros(length(r_num)-1, length(G));

for ii = 1 : length(G)
    PhiTZ(:,ii) = PhiTZi * C(indP(0 : N_sum),ii);
    JTZ(:,ii)   = JTZi * C(indP(0 : N_sum),ii);
    temp = interp1(r_num(1:end-1), JTZ(:,ii), r_num_dense, 'linear', 'extrap');
    temp = filtfilt(lp_filt , temp);
    JTZ_filt(:,ii) = interp1(r_num_dense,temp, r_num(1:end-1));
end


save('V_TZ.mat', 'G_plot', 'C', 'r_num', 'PhiTZi', 'JTZi', 'PhiTZ', 'JTZ','JTZ_filt', '-append');

%%
h_f = figure;

format_axis.XLim = 10.^[-2,2];
format_axis.XTick = 10.^(-2:2);
format_axis.XScale = 'log';
format_axis.XDir = 'reverse';
format_axis.YLim = [0,1];
format_axis.YTick = 0:0.2:1;

h_sp = subplot(20,2,[1,2]);
title(sprintf('Steady state of voltage step response'));
set(h_sp.Title, format_title);
set(h_sp,format_blank_axis);

%
h_sp = subplot(20,2,[3,39]);
surf(G, r_num, PhiSS,'LineStyle' ,'none','FaceAlpha',0.5);
hold on;box on;
plot3(G_plot_lines, r_plot_lines, PhiSS(:,ind_G_plot));


format_axis.ZLim = [0,1.01];
format_axis.ZTick = (0:0.2:1);

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


set(h_f,format_figure,'Position',[0,0,1800,900]);
set(findobj(h_f, 'Type','line'), format_line, 'LineWidth', 1.5,'Color','k');
set(findobj(h_f, 'Type','text'), format_text);

figure_name = 'U0J0_SS';
saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(    imind,cm,fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%%
h_f = figure;
N_plot = 16;

%
h_sp = subplot(1,2,1);
plot(h_sp, r_num, PhiTZi(:, 1:N_plot), 'k');

hold on;
box on;
plot(h_sp, [0,1],[0,0],'k--');

text(h_sp, 0.98, 12, '$i \uparrow$',  'HorizontalAlignment','right');
text(h_sp, 0.98, 8, '$\mathbf{\uparrow}$',  'HorizontalAlignment','right');
plot(h_sp, [0.97,0.97],[-2,8],'k-');

xlabel(h_sp, 'Radial position $r/r_{0}$');
ylabel(h_sp, {'Potential on electrode surface'});

format_axis.XLim = [0,1];
format_axis.XTick = 0:0.2:1;
format_axis.XScale = 'lin';
format_axis.XDir = 'normal';
format_axis.YLim = 60*[-1,1];
if isfield(format_axis,'YTick')
    format_axis = rmfield(format_axis,'YTick');
end

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel], format_axis_label);

%
h_sp = subplot(1,2,2);
plot(h_sp, r_num(1:end-1), JTZi(:,1:N_plot)/(4/pi), 'k');     % normalization factor: 4*B_0/pi=4/pi
hold on;
box on;
plot(h_sp, [0,1],[0,0],'k--');

text(h_sp, 0.98, 325, '$i \uparrow$',  'HorizontalAlignment','right');
text(h_sp, 0.98, 250, '$\mathbf{\uparrow}$',  'HorizontalAlignment','right');
plot(h_sp, [0.97,0.97],[-25,250],'k-');

xlabel(h_sp, 'Radial position $r/r_{0}$');
ylabel(h_sp, {'Current density on electrode surface'});

format_axis.XLim = [0,1];
format_axis.YLim = 1000*[-1,1];

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel], format_axis_label);

set(h_f,format_figure);
set(findobj(h_f, 'Type','line'), format_line);
set(findobj(h_f, 'Type','text'), format_text);

figure_name = 'U0J0_eigen';
saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(    imind,cm,fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%%
h_f = figure;

format_axis.XLim = 10.^[-2,2];
format_axis.XTick = 10.^(-2:2);
format_axis.XScale = 'log';
format_axis.XDir = 'reverse';
format_axis.YLim = [0,1];
format_axis.YTick = 0:0.2:1;
   
h_sp = subplot(20,2,[1,2]);
title(sprintf('Transient of voltage step response'));
set(h_sp.Title, format_title);
set(h_sp,format_blank_axis);

%
h_sp = subplot(20,2,[3,39]);
surf(G, r_num, PhiTZ,'LineStyle' ,'none','FaceAlpha',0.5);
hold on;box on;
plot3(G_plot_lines, r_plot_lines, PhiTZ(:,ind_G_plot));


format_axis.ZLim = [-1,0];
format_axis.ZTick = (-1:0.25:0);

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel,h_sp.ZLabel], format_axis_label);
caxis(h_sp, [min(PhiTZ(:)),max(PhiTZ(:))])
view(h_sp, [75, 25])
xlabel({'Conductance  $G$'});
ylabel('Radial position $r/r_{0}$');
zlabel({'Potential $\varphi^{\mathrm{TZ}}_{0}(r,0^{+})/V_{0}$'});
    
%
h_sp = subplot(20,2,[4,40]);
surf(G, r_num(1:end-1), JTZ_filt,'LineStyle' ,'none','FaceAlpha',0.5);     % normalization factor: 4*B_0/pi=4/pi
hold on;box on;
plot3(G_plot_lines(1:end-1,:), r_plot_lines(1:end-1,:), JTZ_filt(:,ind_G_plot));

format_axis.ZLim = [-4,0.5];
format_axis.ZTick = (-4:0.5:0);

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel,h_sp.ZLabel], format_axis_label);
caxis(h_sp, [-2,0.5])
view(h_sp, [75,25])

xlabel({'Conductance  $G$'});ylabel('Radial position $r/r_{0}$');
zlabel({'Current density  $J^{\mathrm{TZ}}_{0}(r,0^{+})/\overline{J_{0}}$'});


set(h_f,format_figure,'Position',[0,0,1800,900]);
set(findobj(h_f, 'Type','line'), format_line, 'LineWidth', 1.5,'Color','k');
set(findobj(h_f, 'Type','text'), format_text);


figure_name = 'U0J0_TZ';
saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(    imind,cm,fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%%
h_f = figure;
h_sp = subplot(20,2,[1,2]);
title(sprintf('Transient of voltage step response'));
set(h_sp.Title, format_title);
set(h_sp,format_blank_axis);

%
h_sp = subplot(20,2,[3,39]);
surf(G, r_num, PhiTZ,'LineStyle' ,'none','FaceAlpha',0.5);
hold on;box on;
plot3(G_plot_lines, r_plot_lines, PhiTZ(:,ind_G_plot));


format_axis.ZLim = [-1,0];
format_axis.ZTick = (-1:0.25:0);

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel,h_sp.ZLabel], format_axis_label);
caxis(h_sp, [min(PhiTZ(:)),max(PhiTZ(:))])
view(h_sp, [75, 25])
xlabel({'Conductance  $G$'});
ylabel('Radial position $r/r_{0}$');
zlabel({'Potential $\varphi^{\mathrm{TZ}}_{0}(r,0^{+})/V_{0}$'});
    
%
h_sp = subplot(20,2,[4,40]);
h_s = surf(G, r_num(1:end-1), JTZ,'LineStyle' ,'none','FaceAlpha',0.5);     % normalization factor: 4*B_0/pi=4/pi
hold on;box on;
plot3(G_plot_lines(1:end-1,:), r_plot_lines(1:end-1,:), JTZ(:,ind_G_plot));

format_axis.ZLim = [-4,0.5];
format_axis.ZTick = (-4:0.5:0.5);

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel,h_sp.ZLabel], format_axis_label);
caxis(h_sp, [-2,0.5])
view(h_sp, [75,25])

xlabel({'Conductance  $G$'});ylabel('Radial position $r/r_{0}$');
zlabel({'Current density  $J^{\mathrm{TZ}}_{0}(r,0^{+})/\overline{J_{0}}$'});


set(h_f,format_figure,'Position',[0,0,1800,900]);
set(findobj(h_f, 'Type','line'), format_line, 'LineWidth', 1.5,'Color','k');
set(findobj(h_f, 'Type','text'), format_text);

figure_name = 'U0J0_TZ_unfilt';
saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(    imind,cm,fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%%
rmpath(fullfile('..','LegendreBasisMatrix'))
