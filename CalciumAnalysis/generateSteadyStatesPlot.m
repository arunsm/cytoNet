function [] = generateSteadyStatesPlot(cellinfo_extended, threshold, ...
    minFrames, fileName, outputDirName, scaleValues)
%cellinfo_extended : object as generated by generateNPCgraph
%threshold : from 0 to 1, the percentage of the maximum intensity of all
%cells that is categorized as a "high calcium value"
%minFrames : minimum number of consecutive frames needed to categorize as a
%"steady state"
%fileName : name of file being processed (originally supplied to
%generateNPCgraph)
%outputDirName : name of output directory where image is to be saved, as
%used in engine.m
%scaleValues : logical -- should we scale all the stacked graph's axes to
%the max and min values seen in the entire set of intensity values?

%Generate steady state graphs for all cells with frameMin consecutive
% frames above the specified percentage of max intensity.

signalMax = max(cellinfo_extended.rawIntensity(:)); 
signalThresh = threshold * signalMax; 
%one row per cell ... any signals above threshold?
logicalThresh = (cellinfo_extended.rawIntensity >= signalThresh); 
%sum across columns to find cells with any threshold-exceeding frames
idxAbove = find(sum(logicalThresh, 2)); 
signalsAbove = cellinfo_extended.rawIntensity(idxAbove, :); 
%now generate a subset where the number of consecutive signals above
%threshold meets or exceeds our criterion
%Make a copy of signalsAbove, and drop out any items that don't have runs
%of minFrames length. 
% I'm going to be sneaky and leverage the image processing toolbox, per 
% https://www.mathworks.com/matlabcentral/answers/230702-how-to-find-consecutive-values-above-a-certain-threshold
runsToShow = zeros(size(signalsAbove)); 
for i=1:length(idxAbove)
    hasRun = false; 
    binaryVector = signalsAbove(i,:); 
    [labeledVector, numRegions] = bwlabel(binaryVector); 
    potentialRun = cellinfo_extended.rawIntensity(idxAbove(i),:);
    measurements = regionprops(labeledVector, ... 
        cellinfo_extended.rawIntensity(idxAbove(i),:), 'Area', 'PixelValues');
    for k = 1 : numRegions
      if measurements(k).Area >= minFrames
        % Area (length) is 3 or greater, so store the values.
        % ca{k} = measurements(k).PixelValues;
        hasRun = true; 
      end
    end
    if hasRun
        %insert into the set of runs we want to show
        runsToShow(i,:) = potentialRun'; 
    end
end
%and remove any extra rows
runsToShow = runsToShow(any(runsToShow,2),:);

%now generate the actual plots

%for stackedplot, we need data in columns
runsToShow = runsToShow'; 
frames = 1:size(runsToShow, 1); 
%and in case we need to scale, set up our scaling factors
intensityFlat = runsToShow(:); %flatten so we can get overall min/max
intensityMin = min(intensityFlat); 
intensityMax = max(intensityFlat); 
% Show cell numbers as column names
labels = cellstr(num2str(idxAbove)); 
%build plots in groups of n
n = 10
currentPlot = 1; %because someone already used 'plot_idx', cough cough
for i=1:10:size(runsToShow,2)
    figure('visible', 'off');
    set(gcf, 'color', 'w');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [6 2]);
    s = stackedplot(frames, runsToShow(:,i:min(i+n-1, size(runsToShow,2))), ... 
        'DisplayLabels', labels(i:min(i+n-1, size(labels, 1)))); 
    if(scaleValues)
        % set common y limits for our stacked graphs
        for j=1:length(s.AxesProperties)
            s.AxesProperties(j).YLimits = [intensityMin intensityMax]; 
        end
    end
    % copied from the generatePlots code
    [filepath,name,plot_idx] = fileparts(fileName);
    %only add a "_2" to the end if we have more than one graph
    outputFileName = sprintf('%s%s%s_highSteadyState', outputDirName, filesep, name); 
    if currentPlot > 1
       outputFileName = strcat(outputFileName, '_', num2str(currentPlot)); 
    end
    print(outputFileName,'-dtiffn')
    currentPlot = currentPlot + 1; 
end
end
