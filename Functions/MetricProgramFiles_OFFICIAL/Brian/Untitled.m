function [header values] test(data)

header = {'test'};
numcells = numel(data.cellmaskdata);
values = zeros(numcells,1);
values(numcells) = 1/0;




end