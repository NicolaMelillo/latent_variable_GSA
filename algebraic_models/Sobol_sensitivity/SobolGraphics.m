function [ f1 ] = SobolGraphics( input )
%%

%%%
% Nicola Melillo, University of Pavia.
% 28 October 2018
%%%

%% input structure


%%% outplut_SA: output SobolSensitivity
%               {S_mean,S_std,S_05prctile,S_95prctile;...
%                T_mean,T_std,T_05prctile,T_95prctile;...
%                PE_mean, PE_95prctile, PE_05prctile, PE_std};
%%% param_extracted: cell array containing the parameters names
%%% title_bar: title of the figure
%%% fs: figure size -> 0 standard | 1 widescreen
%%% print: print figure -> 0 no | 1 yes
%%% format: string containing the figure format, e.g. '-dpdf', '-dpng',
%           '-djpg', '-depsc' (vectorial). See MATALB documentation for
%           figures formats.
%           if print is equal to 0 this entry can be omitted
%%% print_dir: directory in which we want to print the figure
%              if print is equal to 0 this entry can be omitted
%%% resolution: figure resolution, e.g. '-r150'
%               if print is equal to 0 this entry can be omitted


output_SA = input.output_SA;
param_extracted = input.param_extracted;
title_bar = input.title_bar;
fs = input.fs;
print_x = input.print;
ylimit = input.ylimit;

%% organize input parameters...

interesting_param = 1:length(param_extracted);

if fs == 0
    figure_size = [0.2884 0.3086 0.4217 0.6680];
elseif fs == 1
    figure_size = [0 0 1 1];
end

S_mean = output_SA{1,1};
S_inf = output_SA{1,3};
S_sup = output_SA{1,4};
T_mean = output_SA{2,1};
T_inf = output_SA{2,3};
T_sup = output_SA{2,4};

%% Figure 1 - main and total effect

matrix_toplot = [S_mean(interesting_param),T_mean(interesting_param)];
matrix_prc{1} = [S_inf(interesting_param) S_sup(interesting_param)];
matrix_prc{2} = [T_inf(interesting_param) T_sup(interesting_param)];

f1 = figure();
set(f1,'units','normalized','outerposition',figure_size);
b = bar(matrix_toplot);
b(1).FaceColor = [112 216 231]/255;
b(1).EdgeColor = [0,0,0]/255;
b(1).LineWidth = 1;
b(2).FaceColor = [25,134,150]/255;
b(2).EdgeColor = [0,0,0]/255;
b(2).LineWidth = 1;

set(gca,'XTick',1:length(interesting_param),'FontSize',16)
set(gca,'xticklabel',param_extracted(interesting_param),'FontSize',16)
set(gca,'XTickLabelRotation',45)
set(gca,'TickLabelInterpreter','none')
title(title_bar,'FontSize',20)

pause(0.1); %pause allows the figure to be created
lb = numel(b);

for ib = 1:lb
    matrix_prc_ib = matrix_prc{ib};
    val_inf = matrix_toplot(:,ib) - matrix_prc_ib(:,1);
    val_sup = - matrix_toplot(:,ib) + matrix_prc_ib(:,2);
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = b(ib).XData+b(ib).XOffset;
    hold on
    errorbar(xData',matrix_toplot(:,ib),val_inf, val_sup, 'k.','LineWidth',1);
end

if ~isempty(ylimit)
    ylim(ylimit)
end

leg = legend('Main Effect','Total Effect','location','best');
set(leg,'FontSize',17);

pause(2)

%% print

if print_x == 1
    
    format = input.format;
    print_dir = input.print_dir;
    resolution = input.resolution;
    name_print1 = title_bar;
    
    printpdf( f1, name_print1, print_dir, format, resolution )
end

end

