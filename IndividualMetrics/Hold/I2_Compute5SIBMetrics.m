function [header values] = Compute5SIBMetrics(data)

fprintf('Compute5SIBMetrics: Distorted contours; see testSIB.m\n');
% fprintf('Compute5SIBMetrics: Multiple nuclei\n');
% fprintf('Compute5SIBMetrics: Contour band averages include zeros from outside the band\n');
% fprintf('Compute5SIBMetrics: Remove nucleus stain from this metric\n');

metrics = {'MeanLocationAboveAvg', 'MeanLocation', 'MaxLocation', 'MaxIntensity', 'Slope'};
nummetrics = numel(metrics);
staindata = data.staindata;
numstains = numel(staindata);
stainNames = fieldnames(data.stainkey);
numStainNames = numel(stainNames);
% Make sure that there at least as many stain names as stain images
assert(numstains <= numStainNames, sprintf('Compute5SIBMetrics: %d stain images but only %d stain names', numstains, numStainNames));

% Omit metrics for nucleus stain
nucleusIndex = data.stainkey.nucleus;
if (nucleusIndex > 0) & (nucleusIndex <= numstains)
    switch nucleusIndex
        case 1
            stainIndicesNoNucleus = 2:numstains;
        case numStainNames
            stainIndicesNoNucleus = 1:numstains-1;
        otherwise
            stainIndicesNoNucleus = [1:nucleusIndex-1, nucleusIndex+1:numstains];
    end
else
    stainIndicesNoNucleus = 1:numstains;
end
stainNames = stainNames(stainIndicesNoNucleus);
numStainsNoNucleus = numel(stainIndicesNoNucleus);

header = {};
for i = 1:nummetrics
    h = strcat(metrics{i}, {' '}, stainNames)';
    for j = 1:numel(h)
        header{end+1} = h{j};
    end
end

stainImage = data.staindata;
nuMask = data.nucleusmaskdata;
cyMask = data.cytomaskdata;
numOfCell = numel(nuMask);
numOfStain = numel(stainImage);
set_stain = 1:numOfStain;
aDiv = 100;
rDiv = 10;
[row, col] = size(nuMask{1});

values = zeros(numOfCell, numel(header));

SIB = ComputeSIB(stainImage,nuMask,cyMask,numOfCell,numOfStain,aDiv,rDiv,set_stain,row,col, data.cellmaskdata);
location = linspace(100/rDiv,100,rDiv);
loc80 = .80;
loc20 = .20;
idx80 = round(loc80*rDiv);
idx20 = round(loc20*rDiv);
    
for idx_cell = 1:numOfCell
    for idx_stain = 1:numel(stainIndicesNoNucleus);
        intensity = SIB.intensityavg{idx_cell}(stainIndicesNoNucleus(idx_stain),:);

        idx_location = intensity >= mean(intensity);
        MeanLocationAboveAvg = mean(intensity(idx_location).*location(idx_location))/max(intensity);
        MeanLocation = mean(location.*intensity)/max(intensity);
        MaxLocation = find(intensity == max(intensity),1);
        MaxIntensity = max(intensity);
        Slope = (intensity(idx80)-intensity(idx20))/(loc80-loc20);
        
        values(idx_cell, idx_stain + (0 * numStainsNoNucleus)) = MeanLocationAboveAvg;
        values(idx_cell, idx_stain + (1 * numStainsNoNucleus)) = MeanLocation;
        values(idx_cell, idx_stain + (2 * numStainsNoNucleus)) = MaxLocation;
        values(idx_cell, idx_stain + (3 * numStainsNoNucleus)) = MaxIntensity;
        values(idx_cell, idx_stain + (4 * numStainsNoNucleus)) = Slope;
    end
end
    
end

function SIB = ComputeSIB(stainImage,nuMask,cyMask,numOfCell,numOfStain,aDiv,rDiv,set_stain,row,col, cellMask)

theta = linspace(eps,pi-(eps),round(aDiv/2));
nTheta = length(theta);
maxfractionaldelta = 0;
for idx_cell = 1:numOfCell
    tic
    fprintf('Processing cell #%g:\n',idx_cell)
    
%    cyMask{idx_cell} = cellMask{idx_cell} & ~nuMask{idx_cell};
    cm = (cellMask{idx_cell} & ~nuMask{idx_cell}) | bwperim(nuMask{idx_cell});


    r = (cm ~= cyMask{idx_cell});
    rSum = sum(double(r(:)));
    fprintf('Cell: %d,   %d different pixels\n', idx_cell, rSum);
 
    
    % If there are multiple nuclei, use the largest
    cc = bwconncomp(nuMask{idx_cell});
    largestArea = 0;
    pixelIdxList = [];
    for i = 1:cc.NumObjects
        pil = cc.PixelIdxList{i};
        area = numel(pil);
        if area > largestArea
            largestArea = area;
            pixelIdxList = pil;
        end
    end
    mask = logical(zeros(cc.ImageSize));
    mask(pixelIdxList) = logical(1);
    
    % find the centroid of the nucleus
    statTemp = regionprops(mask,'Centroid');  
    centroid = statTemp.Centroid;    
    xc = round(centroid(1));
    yc = round(centroid(2));
    
    % init
    Xinner = zeros(nTheta,2);
    Yinner = zeros(nTheta,2);
    Xouter = zeros(nTheta,2);
    Youter = zeros(nTheta,2);
    for idx_angle = 1:nTheta
%         X2 = 1:col;
        X = [1 col];
        Y = tan(theta(idx_angle))*(X-xc)+yc;
%         Y2 = tan(theta(idx_angle))*(X2-xc)+yc;
        % find the intersection with the nuclear mask
        [inXnu inYnu] = ext_lineinmask(mask,X,Y);
        if theta(idx_angle)<pi/2
            Xinner(idx_angle,1) = (inXnu(2));
            Xinner(idx_angle,2) = (inXnu(end-1));
            Yinner(idx_angle,1) = (inYnu(2));
            Yinner(idx_angle,2) = (inYnu(end-1));
        else
            Xinner(idx_angle,2) = (inXnu(2));
            Xinner(idx_angle,1) = (inXnu(end-1));
            Yinner(idx_angle,2) = (inYnu(2));
            Yinner(idx_angle,1) = (inYnu(end-1));
        end
        
        % find the intersection with the cytoplasm mask
        [inXcy inYcy] = ext_lineinmask(cyMask{idx_cell},X,Y);
        
        if theta(idx_angle)<pi/2
            
            Xouter(idx_angle,1) = (inXcy(2));
            Xouter(idx_angle,2) = (inXcy(end-1));
            Youter(idx_angle,1) = (inYcy(2));
            Youter(idx_angle,2) = (inYcy(end-1));
        else
            
            Xouter(idx_angle,2) = (inXcy(2));
            Xouter(idx_angle,1) = (inXcy(end-1));
            Youter(idx_angle,2) = (inYcy(2));
            Youter(idx_angle,1) = (inYcy(end-1));
        end
    end
    
    % round the intersection
    Xinner = round(Xinner(:));
    Yinner = round(Yinner(:));
    Xouter = round(Xouter(:));
    Youter = round(Youter(:));
    
    % get the mask band specified by rDiv+1 (number of radial divisions)
    Xset = ext_ndlinspace(Xinner,Xouter,rDiv+1);
    Yset = ext_ndlinspace(Yinner,Youter,rDiv+1);
    masknu = poly2mask(Xset(:,1),Yset(:,1),row,col);
    
    mb = cell(numOfCell,rDiv);
    mb{idx_cell,1} = logical(poly2mask(Xset(:,2),Yset(:,2),row,col)-masknu);
    for idx_band = 2:rDiv
        temp = poly2mask(Xset(:,idx_band),Yset(:,idx_band),row,col);
        mb{idx_cell,idx_band}= logical(poly2mask(Xset(:,idx_band+1),Yset(:,idx_band+1),row,col)-temp);
    end
    
    % get the average stain intensity in a band for every stain
    intensityAvg = zeros(numOfStain,rDiv);
    for idx_stain = set_stain
        for idx_band = 1:rDiv
            bandMask = mb{idx_cell,idx_band};
            % Some masks are not completely connected
            cc = bwconncomp(bandMask);
            numBandPixels = 0;
            for i = 1:cc.NumObjects
                numBandPixels = numBandPixels + numel(cc.PixelIdxList{i});
            end
%             intensity = (mb{idx_cell,idx_band}).*(cyMask{idx_cell}).*(stainImage{idx_stain});
z = bandMask & cyMask{idx_cell};
comp = z ~= bandMask;
diff = sum(comp(:));
if diff ~= 0
%    fprintf('Image: %d Stain: %d Band: %d diff=%d\n', idx_cell, idx_stain, idx_band, diff);
end



            intensity = bandMask.*(cyMask{idx_cell}).*(stainImage{idx_stain});
%             intensityAvg(idx_stain,idx_band) = mean2(intensity);
            intensityAvg(idx_stain,idx_band) = sum(intensity(:)) / numBandPixels;
        end
    end
    SIB.intensityavg{idx_cell} = intensityAvg;
    
    time = toc;
    fprintf('Elapsed time is %3.1f seconds.\n\n',time)
end

end

