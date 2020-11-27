close all;
addpath(fullfile('..','LegendreBasisMatrix'))
load('figure_format.mat');
if ~exist('N','var')|| ~exist('M0','var')|| ~exist('A0','var')
    T_A = tic;
    fprintf('Loading matrices. ');
    load('MatrixData.mat','N','M0','A0');
    fprintf('Time: %s.\n\n',datestr(seconds(toc(T_A)),'MM:SS.FFF'));
end


M1 = M0(indP(1):indP(N),indP(1):indP(N));
A1 =  A0(indP(1):indP(N),indP(1):indP(N));
a0  = A0(indP(0):indP(N),indP(0));
a1  = A0(indP(1):indP(N),indP(0));
%%
log_range = 2;
num_per_dec = 12;

G = [0, logspace(-log_range,log_range, log_range*2*num_per_dec+1 )];

b0SS = zeros(N + 1, length(G));
b1SS_alt = zeros(N, length(G));
VSS = zeros(size(G));

for ii = 1 : length(G)
    b0SS(:,ii)      = (G(ii)*(A0 - 2*(a0*a0')) + M0 ) \ a0;
    b1SS_alt(:,ii) = (G(ii)*(A1 - 2*(a1*a1')) + M1 ) \ a1;
    VSS(ii) = 1./G(ii) + 2 * a0' * b0SS(:,ii);
    
end

per_err_BSS = abs([ones(size(G));b1SS_alt] - b0SS)./abs(b0SS) * 100;
fprintf('Maximum and median percent error of BSS: %1.3g%%, %1.3g%%.\n', max(per_err_BSS(:)), median(per_err_BSS(:)));

save('I_SS.mat','b0SS','VSS','G');

%%
h_f = figure;
h_a = axes;

plot(h_a, G, VSS);
hold on;box on;

format_axis.XTick = 10.^(-2:2);
format_axis.Xlim = [1e-2,1e2];
format_axis.Ylim = [1,max(VSS)];
format_axis.XScale = 'log';
format_axis.YScale = 'log';

set([h_a.XLabel,h_a.YLabel], format_axis_label);
set(h_a.Title, format_title);
set(h_a, format_axis);
ylabel(h_a, 'Normalized steady state voltage $V^{\mathrm{SS}}/V_{0}$');
xlabel(h_a, 'Normalized Faradaic conductance $G$');

set(h_f,format_figure,'Position',[0,0, 1200 600]);
set(findobj(h_f, 'Type','line'), format_line,'Color','k');

figure_name = 'VSS';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(:,51:1150,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');


%%
% Voltage Step Scaling Factor
h_f = figure;
h_a = axes;
plot(h_a, G, VSS.^-1);
box on;hold on;

format_axis.Ylim = [max(VSS)^-1,1];

set([h_a.XLabel,h_a.YLabel], format_axis_label);
set(h_a.Title, format_title);
set(h_a, format_axis);
ylabel(h_a, 'Scaling factor $K$');
xlabel(h_a, 'Normalized Faradaic conductance $G$');


set(h_f,format_figure,'Position',[0,0, 1200 600]);
set(findobj(h_f, 'Type','line'), format_line,'Color','k');

figure_name = 'VoltageStepScalingFactor';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(:,51:1150,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%%
rmpath(fullfile('..','LegendreBasisMatrix'))
