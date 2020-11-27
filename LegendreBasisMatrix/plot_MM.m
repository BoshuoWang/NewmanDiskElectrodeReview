close all;
load('figure_format.mat');
%%
if ~exist('N_basis','var') || ~exist('xi','var') || ~exist('eta','var') || ~exist('U','var')
    T_M = tic;
    fprintf('Loading Basis Functions. ');
    load('BasisFunctions.mat', 'N_basis', 'eta', 'xi', 'U' );
    fprintf('Time: %s.\n\n',datestr(seconds(toc(T_M)),'MM:SS.FFF'));
        
end

N_plot = 32;
m = 1001;
xiM_num = [0,logspace(-4,4,m)]';
MM_num = zeros(length(xiM_num), N_plot);

for n = 0 : N_plot
    fac = factor(U{indP(n)},[eta, xi]);
    
    if n ~= 0
        MM_num(:,indP(n)) = double(vpa(fac(1))) .* double(vpa( subs( fac(2), eta, 1) )) .* ...
            double(vpa( subs( fac(3), xi,  xiM_num ) )) ;
    else
        MM_num(:,indP(n)) = double(vpa(fac(1))) .* double(vpa( subs( fac(2), xi, xiM_num) ));
    end
    ind = find(MM_num(:,indP(n)) < 1e-10 , 1);
    MM_num(ind:end,indP(n)) = 0;
end


%% MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
h_f = figure;
h_a = axes;

loglog(h_a, xiM_num,MM_num(:,1:N_plot),'k' );

hold on;
text(h_a, 8,0.1, '$$M_0$$');
text(h_a, 2.8,0.01, '$$M_2$$');
text(h_a, 2.2,1e-3, '$$M_4$$');
text(h_a, 1.5,1e-4, '$$M_6$$','BackgroundColor','w');
text(h_a, 1.74,2e-6, '$$M_8$$','BackgroundColor','w');
text(h_a, 0.9,2e-6, '$$\mathbf{\cdots}$$');
text(h_a, 0.1,2e-6, '$$M_{32}$$');


format_axis.Xlim = 10.^[-4,4];
format_axis.Ylim = [1e-6,1.0];
ylabel('$$M_{2n}(\xi)$$');
xlabel('$$\xi$$');
title('``Radial" component of the general solution');

set(h_a, format_axis, 'PlotBoxAspectRatio', [4,3,1]);
set([h_a.XLabel,h_a.YLabel], format_axis_label);
set(h_a.Title, format_title);

format_text.FontSize = 14;
set(h_f,format_figure,'Position',[50,50,1200,840]);
set(findobj(h_f, 'Type','line'), format_line);
set(findobj(h_f, 'Type','text'), format_text);

%%
figure_name = 'MM';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(:,51:1100,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');
