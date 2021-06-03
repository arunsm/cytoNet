function [cellSizeMetrics] = generateCellSizePlot(cellinfo_extended, pixelwidth, separations, fileName, ... 
    outputDirName)
%cellinfo_extended : object as generated by generateNPCgraph
%pixelwidth : width of a pixel, in um
%separations : arbitrary list of area separators, in um
%fileName : name of file being processed (originally supplied to
%generateNPCgraph)
%outputDirName : name of output directory where image is to be saved, as
%used in engine.m

%Generate cell size histogram

figure('visible', 'off');
set(gcf, 'color', 'w');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 2]);

avgArea_pixel = cellinfo_extended.avgArea; 
% convert from pixels^2 to um^2 based on pixelwidth (um/pixel side) : 
avgArea_um = avgArea_pixel .* (pixelwidth^2); 
% convert from um^2 area to diameter :
% diameter = 2 * (sqrt(area/pi))
avgDia_um = ((avgArea_um./pi()).^(0.5)) * 2; 
edges = unique([0, separations, max(cellinfo_extended.avgArea)]); 
h = histogram(avgDia_um, edges); 
%don't go all the way up to the top or bottom edge, in case they're way off
%in the distance (ie, [0 20 25 200] will squish the middle edges in the
%graphed output. 
offset = min(diff(edges)); %use the smallest distance between edges ... 
xticks_toset = edges; 
%... as the leading and trailing bin size
xticks_toset(1) = xticks_toset(2) - offset; 
xticks_toset(end) = xticks_toset(end-1) + offset; 
% xticks = xticks_toset; 
%actually, that's not a great idea. Let's do it this way : 
xticks = xticks_toset(2:end-1); 
xlim([min(xticks_toset) max(xticks_toset)]); 

output = 'Cell size metrics :'; 
for i=1:length(h.Values)
    tmpstring = sprintf('%i cells between %i and %i um in diameter.', ... 
        h.Values(i), h.BinEdges(i), h.BinEdges(i+1)); 
    output = [output sprintf('\n') tmpstring]; 
end
cellSizeMetrics = output; 

% copied from the generatePlots code
[filepath,name,ext] = fileparts(fileName);
print(sprintf('%s%s%s_cellareas', outputDirName, filesep, name),'-dtiffn')
end
