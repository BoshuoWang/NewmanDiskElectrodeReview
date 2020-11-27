close all;
load('figure_format.mat');
%%
if ~exist('N','var') || ~exist('N_vec0','var') || ~exist('A0','var') 
    T_A = tic;
    fprintf('Loading A0 matrix. ');
    load('MatrixData.mat','A0','N','N_vec0');
    fprintf('Time: %s.\n\n',datestr(seconds(toc(T_A)),'MM:SS.FFF'));
end

%%
h_f = figure;
set(h_f,format_figure);
set(h_f,'Position',[50,50,1200,800]);
h_a = axes;
surf(h_a, N_vec0,N_vec0,log10(abs(A0))-0.02,'EdgeColor','interp','FaceColor','interp')
hold on;

[C,h_contour]=contour3(h_a, N_vec0(1:N/4),N_vec0(1:N/4),log10(abs(A0(1:N/4,1:N/4))),(-3:-1));
set(h_contour, format_contour);
%  clabel(C,h_contour,'FontSize',13,'Interpreter','latex','LabelSpacing',450);

% [C,h_contour]=contour(N_vec,N_vec,abs(AH),10.^[-6,-6]);
% set(h_contour,format_contour);
 
[C,h_contour]=contour3(h_a, N_vec0,N_vec0,log10(abs(A0)),(-10:-4));
set(h_contour, format_contour);
% clabel(C,h_contour,'FontSize',13,'Interpreter','latex');


format_axis.xlim = [min(N_vec0),max(N_vec0)];
format_axis.ylim = format_axis.xlim;
format_axis.zlim = [-10,0];

set([h_a.XLabel,h_a.YLabel,h_a.ZLabel], format_axis_label);
ylabel(h_a, '$$n$$');
xlabel(h_a, '$$m$$');
zlabel(h_a, '$$|a_{n,m}|$$');
% set(h_a.Title, format_title);
% title('$$|$${\boldmath$${A}$$}$$_{0}(n,m)|$$')

set(h_a, format_axis, 'View', [15 40]);
set(h_a,'XTick',[0,N/8:N/8:N],'YTick',[0,N/8:N/8:N],'YDir','Reverse')

caxis(h_a,[-10,0]);
colormap(h_a,parula(1024));

% h_cb = colorbar(h_a, 'EastOutside');
% set(h_cb, format_color_bar)
% for ii = 1 : length(h_cb.Ticks)
%     h_cb.TickLabels{ii} = sprintf('$$10^{%d}$$',h_cb.Ticks(ii));    
% end
% h_cb.Position(1) = h_cb.Position(1) + 0.05;
% h_a.Position(3) = h_a.Position(3) - 0.05;

h_a.ZTick = (-10:0);
for ii = 1 : length(h_a.ZTick)
    h_a.ZTickLabel{ii} = sprintf('$$10^{%d}$$',h_a.ZTick(ii));    
end
 
% set(H_cb.Title, 'String', '$$|$${\boldmath$${A}$$}$$_{0}(n,m)|$$', format_color_bar_title);

%%
figure_name = 'A0';
% saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
im = frame2im(getframe(h_f));
imwrite(im(26:775,26:1125,:),fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');
