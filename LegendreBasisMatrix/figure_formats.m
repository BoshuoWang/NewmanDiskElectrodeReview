close all;
clearvars 'format_*'
%% Figure format
figure_length_x = 1500;         % Default
figure_length_y = 750;          % Default
format_figure.Position = [50,50,figure_length_x,figure_length_y];
format_figure.Color = 'w';      % background color of figure window


%% Axis format
format_axis.FontSize = 16;
format_axis.TickLabelInterpreter = 'latex';
format_axis.Color = 'w'; % background color of the plot area
format_axis.TickDir = 'out';
format_axis.LineWidth = 1;
format_axis.Box = 'on';

format_axis_label.FontSize = 16;
format_axis_label.Interpreter = 'latex';
format_axis_label.Color = 'k';

format_title.FontSize = 18;
format_title.Interpreter = 'latex';
format_title.Color = 'k';

format_text.FontSize = 15;
format_text.Interpreter = 'latex';


%% Blank placeholder axis
format_blank_axis = format_axis;
format_blank_axis.XColor = 'w';
format_blank_axis.YColor = 'w';
format_blank_axis.ZColor = 'w';
format_blank_axis.TickDir = 'in';
format_blank_axis.XTick = [];
format_blank_axis.YTick = [];
format_blank_axis.ZTick = [];

format_color_bar_title.FontSize = 16;
format_color_bar_title.Interpreter = 'latex';
format_color_bar_title.Color = 'k';
format_color_bar_title.FontWeight = 'Normal';
format_color_bar_title.Units = 'Normalized';
format_color_bar_title.HorizontalAlignment = 'Center';
format_color_bar_title.VerticalAlignment = 'Middle';
% 
format_color_bar.FontSize = 14;
format_color_bar.TickLabelInterpreter = 'latex';
format_color_bar.TickDirection = 'out';
format_color_bar.LineWidth = 1;


%% Contour lines and labels
format_contour.LineStyle = '-';
format_contour.LineWidth = 1.5;
format_contour.LineColor = [0 0 0];

format_line.LineWidth = 1.5;

save('figure_format.mat','-regexp','format_*');
