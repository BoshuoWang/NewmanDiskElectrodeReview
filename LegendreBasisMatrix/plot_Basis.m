%%
close all;
load('figure_format.mat');
%%
T_M = tic;
fprintf('Loading Basis Functions. ');
load('BasisFunctions_U0J0_num.mat','N_basis','N_plot','r_num','J0_num','U0_num');
fprintf('Time: %s.\n\n',datestr(seconds(toc(T_M)),'MM:SS.FFF'));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_f = figure;
h_sp = subplot(20,2,[1,2]);
title({'Potential and current density of basis functions on electrode surface (dimensionless and normalized)'});
set(h_sp.Title, format_title);
set(h_sp,format_blank_axis);

%
h_sp = subplot(20,2,[3,39]);
plot(h_sp,r_num, U0_num(:,indP(1:N_plot)), 'k');

hold on;
box on;
h_U0_0 = plot(h_sp,r_num, U0_num(:,indP(0)), 'k');
plot(h_sp,[0,1],[0,0],'k:');

text(h_sp,0.5,0.95, '$$n = 0$$');
text(h_sp,0.95,-0.45,'$$n = 1$$', 'HorizontalAlignment','Right');
text(h_sp,0.5,-0.45, '$$\leftarrow$$');
text(h_sp,0.2,-0.45, '$$n = 10$$', 'HorizontalAlignment','Right');


xlabel(h_sp,'Radial position $r/r_{0}$');
ylabel(h_sp,{'Potential $P_{2n}(\eta(r,0^{+}))$'});

format_axis.Xlim = [0,1];
format_axis.Ylim = [-1.05,1.05];
 

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel], format_axis_label);


%
h_sp = subplot(20,2,[4,40]);
plot(h_sp,r_num(1:end-1), J0_num(:,indP(1:N_plot))/(4/pi), 'k');     % normalization factor: 4*B_0/pi=4/pi
hold on;
box on;
h_J0_0 = plot(h_sp,r_num(1:end-1), J0_num(:,indP(0))/(4/pi), 'k');     % normalization factor: 4*B_0/pi=4/pi
plot(h_sp,[0,1],[0,0],'k:');

text(h_sp,0.01, -0.8, '$$n = 0$$', 'BackgroundColor', 'w', 'Margin', 1);
text(h_sp,0.03, 7, '$$\uparrow$$', 'BackgroundColor', 'w', 'Margin', 1);
text(h_sp,0.01, 17, '$$n = 10$$');

xlabel(h_sp,'Radial position $r/r_{0}$');
ylabel(h_sp,{'Current density  $P_{2n}(\eta(r,0^{+}))\cdot  M^{\prime}_{2n}(\xi(r,0^{+})) \pi /(4\eta)$'});

format_axis.Xlim = [0,1];
format_axis.Ylim = [-20,20];

set(h_sp, format_axis);
set([h_sp.XLabel,h_sp.YLabel], format_axis_label);

set(h_f,format_figure);
set(findobj(h_f, 'Type','line'), format_line);
set(findobj(h_f, 'Type','text'), format_text);

set([h_U0_0, h_J0_0], 'LineWidth', 2.5);

figure_name = 'U0J0';
saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(    imind,cm,fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_f = figure;
h_a = axes;
box on;hold on;

plot(h_a,r_num(1:end-1), J0_num(:,indP(0))/(4/pi), 'k');     % normalization factor: 4*B_0/pi=4/pi
% plot([0,1],[0,0],'k:');
plot(h_a,r_num,ones(size(r_num)),'--k')

axis equal
format_axis.Xlim = [0,1];
format_axis.Ylim = [0,2];
set(h_a, format_axis,'YTick',[0:0.5:format_axis.Ylim(end)]);
set([h_a.XLabel,h_a.YLabel], format_axis_label);

xlabel(h_a,'Radial position $r/r_{0}$');
ylabel(h_a,'Current density $J_0^{\mathrm{P}}(r)/\overline{J_0}$');


set(h_f,format_figure, 'Position', [50,50,600,800]);
set(findobj(h_f, 'Type','line'), format_line);
set(findobj(h_f, 'Type','text'), format_text, 'VerticalAlignment','Middle', 'HorizontalAlignment','Center');


figure_name = 'JP';
saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(    imind,cm,fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');


%%
T_M = tic;
fprintf('Loading Basis Functions. ');
load('BasisFunctions_U_num.mat','N_plot','R_num','Z_num','U_num','r_num_end','z_num_end');
fprintf('Time: %s.\n\n',datestr(seconds(toc(T_M)),'MM:SS.FFF'));


%%
figure_length_x = 1920;
figure_length_y = 640;
format_figure.Position = [0, 0, figure_length_x, figure_length_y];

x_length = r_num_end * 2;
x_gap = 1;
x_N_gaps = [1,2,2,3,3,3,4,4,4,4];

x_total = x_gap * 5 + sum(unique(x_length));
x_norm_unit = 1/x_total;
x_norm = x_norm_unit * figure_length_x;

x_gap = x_norm_unit * x_gap;
x_length = x_norm_unit * x_length;


y_length = z_num_end;
y_gap = 1;
y_N_gaps = [1,1:2,1:3,1:4];
y_title_space = 0.5;
y_title = 0.25;
y_cb_space = 1;
y_cb = 0.5;

y_total = y_gap * 2 + y_length(indP(0)) + y_title_space + y_cb_space;
y_norm_unit = 1/y_total;
y_norm = y_norm_unit * figure_length_y;

y_gap = y_norm_unit * y_gap;
y_title_space = y_norm_unit * y_title_space;
y_title = y_norm_unit * y_title;
y_cb_space = y_norm_unit * y_cb_space;
y_cb = y_norm_unit * y_cb;
y_length = y_norm_unit * y_length;


axes_pos = cell(size(N_plot));
for ii = 1: length(N_plot) 
    axes_pos{ii} = [x_gap * x_N_gaps(ii) + sum(unique(x_length(1:ii))) - x_length(ii),...
                    1 - y_title_space - (y_gap + y_length(ii) )* y_N_gaps(ii) , ...
                    x_length(ii), y_length(ii)];
end

c_lvl = (-1:0.002:1);

%%
h_f = figure;
h_sp_0 = axes(h_f,'Position',[x_gap, 1 - (y_title_space + y_length(indP(0)) + y_gap + y_cb), sum(unique(x_length)) + x_gap * 3, y_length(indP(0)) + y_gap + y_title]);
set(h_f,format_figure);
set(h_sp_0.Title, format_title);
set(h_sp_0,format_blank_axis);
set([h_sp_0.XLabel,h_sp_0.YLabel], format_axis_label);
title(h_sp_0, {'Potential field distribution (dimensionless and normalized)'},'Interpreter','latex');
xlabel(h_sp_0, '$$r/r_{0}$$','Color','k');    
% ylabel(h_sp_0, '$$z/a$$','Color','k');
caxis(h_sp_0,[-1,1]);
colormap(h_f,bluewhitered(512));
h_cb = colorbar(h_sp_0, 'SouthOutside');
set(h_cb, format_color_bar, 'TickLength',0.001);
h_cb.Position(2) = 1 - (y_title_space + y_length(indP(0)) + y_gap *  2 + y_cb);
%%
for ii = 1 : length(N_plot) 
    h_sp = axes(h_f,'Position',axes_pos{ii});
    hold on;box on;
    
    contourf(h_sp,  R_num{ii}, Z_num{ii}, U_num{ii}(:,:),c_lvl,'LineStyle','none');
    contourf(h_sp, -R_num{ii}, Z_num{ii}, U_num{ii}(:,:),c_lvl,'LineStyle','none');
    [C,h_contour]=contour(h_sp, R_num{ii}, Z_num{ii}, U_num{ii}(:,:),[0,0]);
    set(h_contour, format_contour , 'LineStyle', '--','LineWidth', 0.5);
    [C,h_contour]=contour(h_sp, -R_num{ii}, Z_num{ii}, U_num{ii}(:,:),[0,0]);
    set(h_contour, format_contour, 'LineStyle', '--','LineWidth', 0.5);
    if ii == 1
        [C,h_contour]=contour(h_sp, R_num{ii}, Z_num{ii}, U_num{ii}(:,:),[0.1: 0.1 : 0.6, 0.8]);
        set(h_contour, format_contour , 'LineStyle', '-','LineWidth', 0.5);
        clabel(C,h_contour,'FontSize',12,'Interpreter','latex','LabelSpacing',450);
        [C,h_contour]=contour(h_sp, R_num{ii}, Z_num{ii}, U_num{ii}(:,:),[0.7, 0.9]);
        set(h_contour, format_contour, 'LineStyle', '-','LineWidth', 0.5);
        [C,h_contour]=contour(h_sp, -R_num{ii}, Z_num{ii}, U_num{ii}(:,:),[0.1: 0.1: 0.6, 0.8]);
        set(h_contour, format_contour, 'LineStyle', '-','LineWidth', 0.5);
        clabel(C,h_contour,'FontSize',12,'Interpreter','latex','LabelSpacing',450);
        [C,h_contour]=contour(h_sp, -R_num{ii}, Z_num{ii}, U_num{ii}(:,:),[0.7, 0.9]);
        set(h_contour, format_contour, 'LineStyle', '-','LineWidth', 0.5);
        ylabel(h_sp, '$$z/a$$');
        set(h_sp.YLabel, format_axis_label);
    end
    h_sp.XTick = sort(unique([h_sp.XTick,-1,1]));
    
    plot(h_sp, [-1,1],[0,0],'k-','LineWidth',1)
    axis(h_sp, 'equal');
        
    caxis(h_sp, [-1,1]);
    colormap(h_sp, bluewhitered(512));

    set(h_sp.Title, format_title,'FontSize',14);
    title(h_sp, ['$$n = ',num2str(N_plot(ii)),'$$']);
    
    format_axis.Xlim = [-r_num_end(ii),r_num_end(ii)];
    format_axis.Ylim = [0,z_num_end(ii)];
    
    view(h_sp, 2);
    set(h_sp, format_axis, 'FontSize',14);
    drawnow;
end


figure_name = 'U';
saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
[imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
imwrite(    imind,cm,fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');
