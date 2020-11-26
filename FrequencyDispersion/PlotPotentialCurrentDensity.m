%%
close all;
addpath(fullfile('..','LegendreBasisMatrix'));
%%
load('figure_format.mat');

if (  ~exist('G','var')  || ~exist('Omega','var') || ~exist('bH','var') || ~exist('BH0','var') ...
   || ~exist('Reff','var') || ~exist('Ceff','var') || ~exist('Reff_inf','var') || ~exist('Ceff_inf','var')  )
    load('freq_disp_data.mat','G','Omega','bH', 'BH0','Reff','Ceff','Reff_inf','Ceff_inf');
    fprintf('Loading frequency dispesion results. \n');
end

%%
load('BasisFunctions_U0J0_num.mat','N_basis','N_plot','r_num','J0_num','U0_num');

Omega_plot = 10.^(-3:0.5:3);
[~,ind_omega_plot] = intersect(Omega, Omega_plot);
Omega_plot_lines = kron(Omega_plot,ones(size(r_num))); 
r_plot_lines = repmat(r_num,size(Omega_plot));


PhiH0 = zeros(length(r_num), length(Omega), length(G));
JH0 = zeros(length(r_num)-1, length(Omega), length(G));

format_axis.XLim = 10.^[-3,3];
format_axis.XTick = 10.^(-3:3);
format_axis.XScale = 'log';
format_axis.XDir = 'reverse';
format_axis.YLim = [0,1];
format_axis.YTick = 0:0.2:1;
   
for ii = 1 : length(G)
    for jj = 1 : length(Omega)
        PhiH0(:,jj,ii) = U0_num * bH{end}(indP(0:N_basis), jj, ii);
        JH0(:,jj,ii) = J0_num * bH{end}(indP(0:N_basis), jj ,ii) / bH{end}(indP(0),jj,ii) * pi/4;
    end 
    h_f = figure;
    h_sp = subplot(20,2,[1,2]);
    title(sprintf('Potential and current density for $ G = %g$', G(ii)));
    set(h_sp.Title, format_title);
    set(h_sp,format_blank_axis);

    %
    h_sp = subplot(20,2,[3,39]);
    surf(Omega, r_num, real(squeeze(PhiH0(:,:,ii))),'LineStyle' ,'none','FaceAlpha',0.5);
    hold on;box on;
    plot3(Omega_plot_lines, r_plot_lines, real(squeeze(PhiH0(:,ind_omega_plot,ii)))); 
    
    
    format_axis.ZLim = [0,1.05];
    
    set(h_sp, format_axis);
    set([h_sp.XLabel,h_sp.YLabel,h_sp.ZLabel], format_axis_label);
    caxis(h_sp, [0,1])
    view(h_sp, [75, 25])
    xlabel({'Frequency  $\Omega$'});
    ylabel('Radial position $r/r_{0}$');
    zlabel({'Potential $\varphi^{\mathrm{H}}_{0}(r)/V^{\mathrm{H}}$'});
    
    %
    h_sp = subplot(20,2,[4,40]);
    h_s = surf(Omega, r_num(1:end-1), real(squeeze(JH0(:,:,ii))),'LineStyle' ,'none','FaceAlpha',0.5);     % normalization factor: 4*B_0/pi=4/pi
    hold on; box on;
    plot3(Omega_plot_lines(1:end-1,:), r_plot_lines(1:end-1,:), real(squeeze(JH0(:,ind_omega_plot,ii)))); 
    
    format_axis.ZLim = [0,4];
   
    set(h_sp, format_axis);
    set([h_sp.XLabel,h_sp.YLabel,h_sp.ZLabel], format_axis_label);
    caxis(h_sp, [0.5,4])
    view(h_sp, [75,25])
    
    xlabel({'Frequency  $\Omega$'});
    ylabel('Radial position $r/r_{0}$');
    zlabel({'Current density  $J^{\mathrm{H}}_{0}(r)/\overline{J_0^{\mathrm{H}}}$'});

    
    set(h_f,format_figure,'Position',[0,0,1800,900]);
    set(findobj(h_f, 'Type','line'), format_line, 'LineWidth', 1.5,'Color','k');
    set(findobj(h_f, 'Type','text'), format_text);
    
    figure_name = sprintf('U0J0_G%g',G(ii));
    saveas(h_f,fullfile('Figures',[figure_name,'.fig']));
    [imind,cm] = rgb2ind(frame2im(getframe(h_f)),256);
    imwrite(    imind,cm,fullfile('Figures',[figure_name,'.tif']),'tif','WriteMode','overwrite', 'Resolution',500,'Compression','none');

end


%%
rmpath(fullfile('..','LegendreBasisMatrix'))
