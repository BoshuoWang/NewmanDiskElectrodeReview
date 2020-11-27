close all;
addpath(fullfile('..','LegendreBasisMatrix'));

load('figure_format.mat');
if ~exist('N','var') ||  ~exist('N_vec0','var') || ~exist('M0','var')  || ~exist('A0','var') 
    T_A = tic;
    fprintf('Loading matrices. ');
    load('MatrixData.mat', 'N', 'N_vec0', 'M0', 'A0');
    fprintf('Time: %s.\n\n',datestr(seconds(toc(T_A)),'MM:SS.FFF'));
end


a0  = A0(indP(0):indP(N),indP(0));

%%
No_N = round(2.^(-7:0) * N);
Lambda = cell(size(No_N));
b0TZ = cell(size(No_N));


for kk = 1 : length(No_N)
    T_Lambda = tic;
    fprintf('Calculating eigenvalues of transient response with N = %d.\n', No_N(kk));
    ind = indP(0:No_N(kk));
    A0_N = A0(ind,ind);
    M0_N = M0(ind,ind);
    b0TZ{kk} = ones(No_N(kk)+1);
    
    AE = M0_N / ( A0_N );
    Lambda{kk} = sort(eig(AE,'vector'));
    
    ind = indP(1:No_N(kk));
    a1_N = a0(ind);
    A1_N = A0(ind,ind);
    M1_N = M0(ind,ind);
    for ii = 0: No_N(kk)
        b0TZ{kk}(ind, indP(ii)) = -(A1_N - M1_N / Lambda{kk}(indP(ii)) ) \ a1_N ;
    end
    
    fprintf('Time: %s.\n',datestr(seconds(toc(T_Lambda)),'MM:SS.FFF'));

end
save('V_TZ.mat','No_N','Lambda','b0TZ');


%%
per_acc = nan(size(No_N));
N_acc = nan(size(No_N));

h_f = figure;
h_sp = subplot(2,1,1);

hold on; box on;

for ii = 1 : length(No_N)
    if ii ~= length(No_N)
        ind_inacc = abs( Lambda{ii} ./ (Lambda{ii+1}( indP(0 : No_N(ii)))) - 1 )> 1e-4;   %0.001%
        per_acc(ii) = 1 - sum(ind_inacc)/indP(No_N(ii));
    else
        per_acc(ii) = nanmean(per_acc);
        ind_inacc = ( (1:No_N(ii)) > floor(per_acc(ii)*indP(No_N(ii))));
    end
    N_acc(ii) = indP(No_N(ii)) - sum(ind_inacc);
    
    plot(h_sp, N_vec0(ind_inacc), Lambda{ii}(ind_inacc),'-rx','MarkerIndices',(length(Lambda{ii}(ind_inacc)):-8:1),'Markersize',10);
end

plot(h_sp, N_vec0(1:N_acc(end)),Lambda{end}(1:N_acc(end)),'-b.','MarkerIndices',1:8:N_acc(end),'Markersize',15);


stats = regstats(Lambda{end}( indP( 0 : N_acc(end)) ), N_vec0(indP( 0 : N_acc(end)) ),'linear');

Interc = stats.beta(1);
Slope = stats.beta(2);  
Lambda_lin = Interc + Slope * N_vec0;
res_r2 = 1-stats.rsquare;
expon = floor(log10(res_r2));
signif = res_r2/10^expon;

save('V_TZ.mat', 'Interc', 'Slope', 'N_acc', '-append');

h_l = plot(h_sp, N_vec0 ,Lambda_lin,'--','Color',[1,1,1]*0.1);

Yscale = 100*ceil(Lambda_lin(end)/100);
% scale_p = 2^ceil(log2(Lambda_lin(end)));

plot(h_sp,No_N(end)*0.85,Yscale*0.5,'.b','Markersize',15);
plot(h_sp,No_N(end)*0.85,Yscale*0.4,'xr','Markersize',10);
plot(h_sp,No_N(end)*[0.825,0.875],Yscale*0.3*[1,1],'k--','Markersize',10);

text(h_sp,No_N(end)*0.98,Yscale*0.5,'Accurate','HorizontalAlignment','Right');
text(h_sp,No_N(end)*0.98,Yscale*0.4,'Inaccurate','HorizontalAlignment','Right');
text(h_sp,No_N(end)*0.98,Yscale*0.3,'Linear fit','HorizontalAlignment','Right');

text(h_sp, No_N(end)*0.98,Yscale*0.2,{['$\Lambda^{(i)} = ',num2str(Interc,'%1.5f'),'+',num2str(Slope/(pi/2)^2,'%1.5f'),'\times (\pi/2)^2 i$'],...
                                     ['$R^{2} = 1-',num2str(signif,'%1.3g'),'\times10^{',num2str(expon,'%d'),'}$']},...
                                     'HorizontalAlignment','Right');

title(h_sp,'Eigenvalues of the voltage step response');
ylabel(h_sp,'$\Lambda^{(i)}$');
xlabel(h_sp,'Order $i$');
format_axis.XLim = [0,No_N(end)];
format_axis.XTick = [0,No_N];
format_axis.YLim = [0,Yscale];
format_axis.XScale = 'lin';
format_axis.YScale = 'lin';

set(h_sp,format_axis);
format_axis.XTickLabelMode = 'manual';
format_axis.XTickLabels = h_sp.XTickLabels;
format_axis.XTickLabels(2:3) = {' '};

set(h_sp,format_axis)

set([h_sp.XLabel,h_sp.YLabel], format_axis_label);
set(h_sp.Title, format_title);

h_sp = subplot(2,1,2);
hold on; box on;

equspr = unique(round(2.^(0:0.2:10)));

for ii = 1 : length(No_N)
    [~,mk_ind] = intersect(N_vec0(indP(N_acc(ii):No_N(ii))),equspr);
    mk_ind = unique([mk_ind',No_N(ii)-N_acc(ii)+1]);

    plot(h_sp, N_vec0(indP(N_acc(ii):No_N(ii))), 1./Lambda{ii}(indP(N_acc(ii):No_N(ii))),'-rx','MarkerIndices',(mk_ind),'Markersize',10);
end
[~,mk_ind] = intersect(N_vec0(indP(0:N_acc(end))),equspr);
mk_ind = [mk_ind', N_acc(end)];
plot(h_sp, N_vec0(indP(0:N_acc(end))),1./Lambda{end}(indP(0:N_acc(end))),'-b.','MarkerIndices',mk_ind,'Markersize',15);

h_l = plot(h_sp, N_vec0 ,1./Lambda_lin,'--','Color',[1,1,1]*0.1);

% plot(h_sp,No_N(end)*0.4,1*0.5,'.b','Markersize',15);
% plot(h_sp,No_N(end)*0.4,1*0.5^2,'xr','Markersize',10);
% 
% text(h_sp,No_N(end)*0.98,1*0.5,'Accurate','HorizontalAlignment','Right');
% text(h_sp,No_N(end)*0.98,1*0.5^2,'Inaccurate','HorizontalAlignment','Right');


ylabel(h_sp,'$(\Lambda^{(i)})^-1$');
xlabel(h_sp,'Order $i$');
format_axis.XLim = [1/(2)^(1/4),No_N(end)*(2)^(1/4)];
format_axis.XTick = [1,No_N];
if isfield(format_axis,'XTickLabels')
    format_axis = rmfield(format_axis,'XTickLabels');
end
format_axis.XTickLabelMode = 'auto';
format_axis.XMinorTick ='off';
format_axis.YLim = [1e-7,1];
format_axis.XScale = 'log';
format_axis.YScale = 'log';

set(h_sp, format_axis);

set([h_sp.XLabel,h_sp.YLabel], format_axis_label);
set(h_sp.Title, format_title);

%%
set(h_f,format_figure,'Position',[0,0, 1200 1000]);
set(setdiff(findobj(h_f, 'Type','line'),h_l), format_line,'LineWidth',1);
set(findobj(h_f, 'Type','text'), format_text, 'VerticalAlignment','Middle');


figure_name = 'lambda';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(26:975,51:1150,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

%%
rmpath(fullfile('..','LegendreBasisMatrix'))
