%%
close all;
addpath(fullfile('..','LegendreBasisMatrix'))
load('figure_format.mat');

if (  ~exist('G','var')  || ~exist('Omega','var') || ~exist('bH','var') || ~exist('BH0','var') ...
   || ~exist('Reff','var') || ~exist('Ceff','var') || ~exist('Reff_inf','var') || ~exist('Ceff_inf','var')  )
    load('freq_disp_data.mat','G','Omega','bH', 'BH0','Reff','Ceff','Reff_inf','Ceff_inf');
    fprintf('Loading frequency dispesion results. \n');
end

%%
h_f = figure;

h_sp = subplot(1,2,1);

plot(h_sp, Omega, squeeze(Reff(:,:,end)));

axis tight;hold on;box on;

set([h_sp.XLabel,h_sp.YLabel], format_axis_label);
set(h_sp.Title, format_title);

format_axis.Xlim = [1e-3,1e3];
format_axis.Ylim = [0.9,111];
format_axis.XTick = 10.^(-3:3);
format_axis.XScale = 'log';
format_axis.YScale = 'log';
set(h_sp, format_axis);

title(h_sp,'Normalized resistive impedance');
ylabel(h_sp,'$\Re(Z^{\mathrm{H}}) =R^{\mathrm{H}}_{\mathrm{S}}/R_{\mathrm{S}}$');
xlabel(h_sp,'Normalized frequency $\Omega$');
text(h_sp, 2e-3, 0.99, '$G=0$');
% plot_arrow(0.2,0.95,0.2,1.08,'headwidth', 0.625,'headheight', 0.25);
text(h_sp, 1.5e-2, 50, '$G=0.01$');
text(h_sp, 1.5e-1, 5, '$G=0.1$');
text(h_sp, 1.5e0,  1.5, '$G=1$');
text(h_sp, 1.5e1,  1.2, '$G=10$');
% plot_arrow(2.4,1.26,2.4,1.13,'headwidth', 0.625,'headheight', 0.25);

%
h_sp = subplot(1,2,2);
Z_cap = squeeze(Ceff(:,:,end))./kron(ones(size(G)),Omega');
loglog(h_sp, Omega, Z_cap);
axis tight;hold on;box on;
h_l = loglog(h_sp, Omega, Omega.^-1,'k--');

set([h_sp.XLabel,h_sp.YLabel], format_axis_label);
set(h_sp.Title, format_title);

format_axis.Xlim = [1e-3,1e3];
format_axis.Ylim = [1e-3,1e3];
format_axis.XScale = 'log';
format_axis.YScale = 'log';

set(h_sp, format_axis);

title(h_sp,'Normalized capacitive impedance');
ylabel(h_sp,'$\Im(Z^{\mathrm{H}}) = (C_{\mathrm{DL}}/C^{\mathrm{H}}_{\mathrm{DL}})/ \Omega$');
xlabel(h_sp,'Normalized frequency $\Omega$');
text(h_sp,0.002,700, '$G=0$');
text(h_sp,0.0012, 7, '$G=0.01$');
text(h_sp,0.012, 0.7, '$G=0.1$');
text(h_sp,0.12, 0.07, '$G=1$');
text(h_sp,1.2, 0.007, '$G=10$');

set(h_f,format_figure,'Position',[0,0, 1500 750]);
set(setdiff(findobj(h_f, 'Type','line'),h_l), format_line,'Color','k');
set(findobj(h_f, 'Type','text'), format_text, 'VerticalAlignment','Middle', 'HorizontalAlignment','Left');


figure_name = 'EffectiveRC';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(:,101:1400,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%%

slope = diff(log10(Z_cap))./diff(log10(kron(ones(size(G)),Omega')));
d_Omega = (Omega(1:end-1) + Omega(2:end))/2;

[max_slope_G0, ind_max] = max(slope(:,1));

h_f = figure;
h_a = gca;
semilogx(h_a, d_Omega,slope)

hold on;box on;

format_axis.XTick = 10.^(-5:5);
format_axis.Xlim = [1e-5,1e5];
format_axis.Ylim = [-1.05,1.05];
format_axis.XScale = 'log';
format_axis.YScale = 'linear';

set(h_a, format_axis);
set([h_a.XLabel,h_a.YLabel], format_axis_label);
set(h_a.Title, format_title);
set(h_f,format_figure,'Position',[0,0, 1200 800]);

title(h_a,'log-log slope of normalized capacitive impedance');
ylabel(h_a,'$\mathrm{d} (\mathrm{lg}(\Im(Z^{\mathrm{H}})) / \mathrm{d} (\mathrm{lg}(\Omega))$');
xlabel(h_a,'Normalized frequency $\Omega$');
text(h_a,1e-4,-0.95, '$G=0$');
text(h_a,1e-2, 0.2, '$G=0.01$','Rotation',-80);
text(h_a,1e-1, 0.2, '$G=0.1$','Rotation',-80);
text(h_a,1e0, 0.2, '$G=1$','Rotation',-80);
text(h_a,1.2e1, 0.2, '$G=10$','Rotation',-80);


set(setdiff(findobj(h_f, 'Type','line'),h_l), format_line, 'Color','k');
set(findobj(h_f, 'Type','text'), format_text, 'VerticalAlignment','Middle', 'HorizontalAlignment','Left');

plot(h_a,[d_Omega(1),d_Omega(ind_max)],max_slope_G0*[1,1],'k--');
h_t = text(h_a,1e-4,max_slope_G0, 'max slope $G=0$');
set(h_t, format_text, 'VerticalAlignment','Bottom', 'HorizontalAlignment','Left');

figure_name = 'loglogslope_C';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(:,51:1150,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%%
rmpath(fullfile('..','LegendreBasisMatrix'))
