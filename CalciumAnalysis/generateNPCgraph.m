%Generates a graph overlaid on the original input image based on the
%correlation coefficients of the mean pixel intensity in the regions of
%interest found by the SegmentNPC function.

%Input: filename-image file you wish to generate a graph for; should
%contain multiple image files.

%Input: timetotal-total timeframe image set covers

%Input: timestop-timepoint you wish to stop reading in data. 

function [errorReport, cellinfo, framestop, bestsig1, bestsig2, badsig1, badsig2, adj, corrs, ...
    cutoff, ccov, bw, image, cellinfo_extended] = generateNPCgraph(filename, timetotal, timestop)

close all
errorReport = [];

%Segment median image, find ROI
info = imfinfo(filename);
framestop = round(numel(info)*(timestop/timetotal));
[image, bw] = segmentNPC(filename);
CC = bwconncomp(bw);
numcells = CC.NumObjects;
fprintf('Image Segmented \n')

if numcells == 0
    fn = getFileName(filename);
    errorReport = makeErrorStruct(['no cells found in file ', fn], 1);
    cellinfo = [];
    framestop = [];
    bestsig1 = [];
    bestsig2 = [];
    badsig1 = [];
    badsig2 = [];
    adj = [];
    corrs = [];
    cutoff = [];
    ccov = [];
    bw = [];
    image = [];
    cellinfo_extended = [];
    return;
end

%Gather info on ROI
% ZM : cellinfo(n,o,p) is the information about the nth cell, with p at o=1,2
% being the x-y location and p at o=3 holding data about the intensity at that
% location. The index of p reflects the frame, so for a 15-layer TIF, we will
% have p = 1:15.
% This means that p only has different data for o=3, since for o=1, p=the x
% coordinate (for each of the frames, it is the same) and for o=2, p=the y
% coordinate (again, the same for each frame). 
% At this point, we are looking at the segmentNPC output, which is
% image-specific, rather than frame-specific.
cellinfo = zeros(numcells, 3,framestop);

% ZM : Additionally output a struct array with derived information.
% Fields are : 
%       avgArea (in pixels; average area over all frames for each cell)
cellinfo_extended = struct; 
cellinfo_extended.totalCellCount = numcells;
cellinfo_extended.rawIntensity = zeros(numcells, framestop); %just a placeholder now

stats = regionprops(CC, 'Centroid');
for k = 1:CC.NumObjects
        current = stats(k);
        location = current.Centroid;
        cellinfo(k, 1, :) = location(1);
        cellinfo(k, 2, :) = location(2);
end

% Here, we start looking at each frame. For cellinfo(n,o,p), we now load
% information about the nth cell: 1:p frames into the o=3 dimension.

% ZM : to get mean areas for each cell in cellinfo_extended, we have to loop
% through all the frames and get the area for the cell in that frame. 
areas = zeros(size(cellinfo, 1), framestop); 
parfor j=1:framestop
    im = toGray(imread(filename, j));
    %Subtract background intensity
    Fback = mean2(im(~bw));
    im = im - Fback;
    %Gather Pixel intensities and area for each cell
    stats = regionprops(CC, im, 'PixelValues', 'Area');
    %Store mean intensity in cellinfo
    intensities = zeros(size(cellinfo(:,3,j)));
    % ZM : build a similar array of current-frame areas
    currentAreas = intensities; 
% ZM : CC.NumObjects is already loaded into numcells.
    for k = 1:CC.NumObjects
        current = stats(k);
        intensities(k) = sum(current.PixelValues) / current.Area;
        currentAreas(k) = current.Area; 
    end
% The idea here is that for cellinfo(n,o,p), o=3 holds the intensities in p
% for each of the frames (from 1:p) for the nth cell. 
    cellinfo(:,3,j) = intensities;
    areas(:, j) = currentAreas; 
end

%Process data
for j = 1:numcells
    %Subtract linear trend from mean intensity data
    cellinfo(j,3,:) = detrend(squeeze(cellinfo(j,3,:))) + cellinfo(j,3,1);
    %ZM : Save a copy of the intensity data before processing
    cellinfo_extended.rawIntensity(j,:) = cellinfo(j,3,:); 
    %Create moving baseline vector
    F0 = squeeze(cellinfo(j,3,:));
    region = round(length(F0)/10);
    for k = 1:length(F0)
        if k <= region
            F0(k) = min(F0(1:k+region));
        elseif k >= length(F0) - region
            F0(k) = min(F0(k-region:length(F0)));
        else 
            F0(k) = min(F0(k-region:k+region));
        end
    end
    %Subtract baseline from intensity, then divide by baseline
    % ZM - this is the delta F / F0
    tmpF = (squeeze(cellinfo(j,3,:))-F0)./F0;
    cellinfo(j,3,:) = tmpF; 
end
fprintf('Data Processed \n')

%Filter out inactive cells to avoid false positives and measure activity
totalactivity = 0;
delete = pi*ones(size(cellinfo(1,3,:))); %marker vector
for j = 1:numcells
    signal = squeeze(cellinfo(j,3,:));
    spikes = findpeaks(signal, 'MinPeakProminence', abs(mean(signal))); %Spike is peak that drops by the mean of the signal
    if length(spikes) < 1
        cellinfo(j,3,:) = delete;
        numcells = numcells - 1;
    end
    totalactivity = totalactivity + length(spikes);
end
%Delete data from inactive cells
j = 1;
j_allcells = 1; % ZM index to track which of the allcells we are on
while size(cellinfo,1) > numcells
    if cellinfo(j,3,:) == delete
        %ZM : thanks, I hate this syntax. Anyway, it deletes the jth
        %indexed entry entirely.
        cellinfo(j,:,:) = [];
        %ZM : I am implementing the same delete-inactive-entries method with
        %all of our other cell infos.
        areas(j,:) = [];
        cellinfo_extended.rawIntensity(j,:) = [];
    else
        j = j + 1;
    end
    j_allcells = j_allcells + 1; 
end
%average area by getting row-wise mean of areas (active cells only, now)
cellinfo_extended.avgArea = mean(areas,2); 
cellinfo_extended.centroid = cellinfo(:,1:2,1); 
fprintf('Inactive Cells Removed \n')


% %ZM : For sanity checking, let's show our image and overlay centroids of
% %each cell. 
% % Make a copy of cellinfo with just the active cells
% imshow(image,'InitialMagnification',200); hold on; 
% scatter(cellinfo(:,1,1), cellinfo(:,2,1), '+')

fprintf('Data Extracted \n')

%Generate adjacency matrix based on cross-correlation of pixel intensities
%over time, and keep track of two signal pairs to show good and
adj = zeros(numcells, numcells);
corrs = [];
[ccov, cutoff] = generateCutoff(cellinfo);
fprintf('Cutoff Generated\n')
bestsig1 = 1;
bestsig2 = 1;
badsig1 = 1;
badsig2 = 1;
for j = 1:numcells
    correlations = zeros(length(adj), 1);
    parfor k = j:numcells
%         current= max(xcorr(cellinfo(j,3,:), cellinfo(k,3,:), 'coeff'));
          current= max(xcorr(squeeze(cellinfo(j,3,:)), squeeze(cellinfo(k,3,:)), 'coeff'));
%         current= max(xcov(cellinfo(j,3,:), cellinfo(k,3,:), 'coeff'));
        if j == k
            current = 0;
        end
        correlations(k) = current;
    end
    adj(j,:) = correlations;
end
for j = 1:numcells
    for k = j:numcells
        adj(k,j) = adj(j,k);
        corrs = [corrs adj(j,k)];
        if adj(j,k) > adj(bestsig1, bestsig2)
            bestsig1 = j;
            bestsig2 = k;
        end
        if adj(j,k) > adj(badsig1, badsig2) & adj(j,k) < (cutoff / 3);
            badsig1 = j;
            badsig2 = k;
        end
    end
end
fprintf('Adjacency Matrix Created \n')

end

%Generates a cuttoff correlation value by scrambling each data set, finding
%all correlation coefficients, and then using the 99th percentile of the
%set of coefficients as the cutoff value
function [ccov, cutoff] = generateCutoff(cellinfo)
s = size(cellinfo);

%Scrambles data, preserving total activity 
for j = 1:s(1)
    t = randi(s(3));
    front = cellinfo(j,3,1:t-1);
    back = cellinfo(j,3,t:s(3));
    index = s(3) - t;
    cellinfo(j, 3, 1:index + 1) = back;
    cellinfo(j, 3, index + 2:s(3)) = front;
end
ccov = [];
sigdevs = std(cellinfo(:,3,:), 0 , 3);
sigmeans = mean(cellinfo(:,3,:), 3);
%Generate vector of Correlation Coefficients
parfor j = 1:s(1)
    for k = j+1:s(1)
        if j ~= k
%             ccov = [ccov max(xcorr(cellinfo(j,3,:), cellinfo(k,3,:), 'coeff'))];
%             ccov = [ccov max(xcov(cellinfo(j,3,:), cellinfo(k,3,:), 'coeff'))];
            ccov = [ccov max(xcorr(squeeze(cellinfo(j,3,:)), squeeze(cellinfo(k,3,:)), 'coeff'))];
        end
    end
end
%Makes cutoff 99th percentile of the coefficients
cutoff = prctile(ccov, 99);
end

%Appends a tag to a filename after removing file extension from filename
function nm = outputFileName(fileName, tag)
k = strfind(fileName, '.');
if isempty(k)
    nm = [fileName, tag];
else
    nm = [fileName(1:(k(end)-1)), tag];
end
end

