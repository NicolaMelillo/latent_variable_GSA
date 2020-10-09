function [  ] = printpdf( h, name_print, print_dir, format, resolution )

%%%
% Version 25 October 2018
% by Nicola Melillo
% BMS lab - University of Pavia
% Copyright 2018
%%%

if ~exist(print_dir, 'dir')
    mkdir(print_dir)
end

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
name = strcat(print_dir,'/',name_print);
print(h,name,format,resolution)

end