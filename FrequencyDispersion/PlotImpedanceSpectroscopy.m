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
ZH = 1./squeeze(BH0(:,:,end));
ZH(1,1) = inf/1i;


h_f = figure;

h_sp(1) = subplot(2,1,1);

loglog(h_sp(1), Omega, abs(ZH));
axis tight;hold on;box on;

set([h_sp(1).XLabel,h_sp(1).YLabel], format_axis_label);
set(h_sp(1).Title, format_title);

format_axis.Xlim = [1e-3,1e3];
format_axis.Ylim = [5e-1,1.5e3];
format_axis.YTick = [1 10 100 1000];
set(h_sp(1), format_axis);

ylabel(h_sp(1),{'Impedance magnitude $|Z^{\rm{H}}|/R_{\rm{S}}$'});
title(h_sp(1),'Impedance spectroscopy of disk electrode')


text(h_sp(1), 1.1e-3,1100 , '$G=0$');
text(h_sp(1), 1.1e-3,140, '$G=0.01$');
text(h_sp(1), 1.1e-3,14, '$G=0.1$');
text(h_sp(1), 1.1e-3,2.6, '$G=1$');
text(h_sp(1), 1.1e-3,1.4, '$G=10$');
h_l = plot(h_sp(1),  Omega, ones(size(Omega)),'k--');


h_sp(2) = subplot(2,1,2);
semilogx(h_sp(2), Omega, angle(ZH).*180/pi);

axis tight;hold on;box on;

set([h_sp(2).XLabel,h_sp(2).YLabel], format_axis_label);

format_axis.Xlim = [1e-3,1e3];
format_axis.Ylim = [-90,0];
format_axis.YTick = [-90 -60 -30 0];
for ii = 1 : length(format_axis.YTick)
    format_axis.YTickLabels{ii} = num2str(format_axis.YTick(ii),'$$%d ^{\\circ}$$');
end
set(h_sp(2), format_axis);

ylabel(h_sp(2), {'Impedance phase $\angle Z^{\rm{H}}$'});
xlabel(h_sp(2), 'Normalized frequency $\Omega$');

text(h_sp(2), 1.1e-3,-85, '$G=0$');
text(h_sp(2), 1.1e-3,-25, '$G=0.01$');
text(h_sp(2), 1.1e-2,-20, '$G=0.1$');
text(h_sp(2), 1.1e-1,-10, '$G=1$');
text(h_sp(2), 1.1e-0, -5, '$G=10$');


set(h_f,format_figure,'Position',[0,0, 1200 800]);
set(setdiff(findobj(h_f, 'Type','line'),h_l), format_line,'Color','k','LineWidth',1.5);

%%

opt = optimset('fminsearch');
opt = optimset(opt,'MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-12,'TolFun',1e-12);


for ii = 1 : length(G)
    fprintf('Fitting spectrum for G = %3.3g. \n', G(ii));
    fprintf('\tFitting with CDL.\n');
    para_init = [1 , 1 , 1/(G(ii) + eps)];
    [para, Jmin, ExitFlag] = fminsearch(@(x) JError3(x,Omega',ZH(:,ii)), para_init,opt);
    Fit.R_S = para(1);
    Fit.C_DL= para(2);
    Fit.R_CT = para(3);
    fprintf('\t\tFitted G: %1.4g. Fitting error: %1.3e.\n',Fit.R_S/Fit.R_CT , Jmin);
    fprintf('\t\tFitted R_S: %3.5f. Fitted R_CT:  %1.4e. Fitted C_DL: %1.4f. \n',Fit.R_S, Fit.R_CT, Fit.C_DL);
%     if ii == 1
%         Zfit = Fit.R_S + 1./(1i * Omega * Fit.C_DL);
%     else
        Zfit = Fit.R_S + Fit.R_CT./(1 + 1i * Omega * Fit.C_DL * Fit.R_CT);
%     end
    loglog(h_sp(1), Omega, abs(Zfit), '--r');
    semilogx(h_sp(2), Omega, angle(Zfit).*180/pi, '--r');
    
    fprintf('\tFitting with CPE.\n');
    para_init = [1 , 1 , 1, 1/(G(ii) + eps)];
    [para, Jmin, ExitFlag] = fminsearch(@(x) JError4(x,Omega',ZH(:,ii)), para_init,opt);
    Fit.R_S = para(1);
    Fit.Q = para(2);
    Fit.n = para(3);
    Fit.R_CT = para(4);
    fprintf('\t\tFitted G: %1.4g. Fitting error: %1.3e.\n', Fit.R_S/Fit.R_CT, Jmin);
    fprintf('\t\tFitted R_S: %1.4f. Fitted R_CT:  %1.4e. Fitted Q: %1.4f. Fitted n: %1.4f. \n\n',Fit.R_S, Fit.R_CT, Fit.C_DL, Fit.n);
    
    Zfit = Fit.R_S + Fit.R_CT./(1 + (1i * Omega).^Fit.n * Fit.Q * Fit.R_CT);

    loglog(h_sp(1), Omega, abs(Zfit), '--b');
    semilogx(h_sp(2), Omega, angle(Zfit).*180/pi, '--b');
end

plot(h_sp(2),10*[1,2],-40*[1,1], '-k','LineWidth',1.5);
plot(h_sp(2),10*[1,2],-50*[1,1], '--r');
plot(h_sp(2),10*[1,2],-60*[1,1], '--b');

text(h_sp(2),10*3,-40,'Solution');
text(h_sp(2),10*3,-50,'Fit with $C_{\mathrm{DL}}$');
text(h_sp(2),10*3,-60,'Fit with $\mathrm{CPE}$');

set(findobj(h_f, 'Type','text'), format_text, 'VerticalAlignment','Middle', 'HorizontalAlignment','Left');


figure_name = 'ImpedanceSpectrum';
saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(    imind,cm,fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%%
rmpath(fullfile('..','LegendreBasisMatrix'))

%%
function J = JError3(x, omega, Z)   %R_p+(CDL //Rp) model, relative error
% x(1) = R_s  x(2) = C_DL  x(3) = R_p
Z_fit = x(1) + x(3)./(1 + 1i * omega * x(2) * x(3));

J = sqrt( nanmean( ( abs(Z_fit - Z)./abs(Z) ).^2 ) );

if any(x <= 0)
    J = inf;
end

end

function J = JError4(x, omega, Z)   %R_p+(CPE //Rp) model, relative error
% x(1) = R_s  x(2) = Q  x(3) = n	x(4) = R_p
% ZCPE = ((j*w)^n * Q)^(-1)
Z_fit = x(1) + x(4)./(1 + (1i * omega) .^  x(3) * x(2) * x(4));

J = sqrt( nanmean( ( abs(Z_fit - Z)./abs(Z) ).^2 ) );

if any(x <= 0) || x (3) > 1
    J = inf;
end

end
