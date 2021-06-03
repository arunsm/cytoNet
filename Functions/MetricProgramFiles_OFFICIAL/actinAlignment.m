%actinAlignment() takes an image file of a single cell.
% The image should contain only the data from the actin-channel
%Returned is... statstics on spread - NEED TO DECIDE

% actinAlignment fuction now returns individual pixel angles and 
% suggested number of bins for histogramming
% 1/17/2013 Byron
%function [FitResults LowestError] = actinAlignment(imagefile)
function [actinThetas numBins] = actinAlignment(imagefile)
info = imfinfo(imagefile);
I = imreadG16(imagefile,info,1);
I = im2double(I);

H1 = fspecial('sobel'); %returns the matrix that is the sobel filter
                       % (horizontal, the transpose is the vertical)
Gh = imfilter(I,H1); %apply horizontal sobel filter
Gv = imfilter(I,H1'); %apply vertical sobel filter

%Calculate angles (theta) and edge intensities (R)
theta = atan2(Gv,Gh).*(180/pi); %REVERSED ratio from paper, returns degrees
R = (Gh.^2 + Gv.^2).^0.5;


%force interval of theta to be [-90,90] intead of [-180,180] originally
thetaAbove90Bool = theta > 90;
thetaBelowNeg90Bool = theta < -90;
thetaCorrectBool = (theta <= 90) & (theta >= -90);
theta = thetaCorrectBool.*theta + thetaAbove90Bool.*((thetaAbove90Bool.*theta)-180) + thetaBelowNeg90Bool.*((thetaBelowNeg90Bool.*theta)+180);


% Change anything "near" -90 to 90 (what is "near" still needs to be
% determined)
thetaNearNeg90 = theta==-90;
thetaNotNearNeg90 = ~thetaNearNeg90;
theta = thetaNotNearNeg90.*theta + thetaNearNeg90.*90;


%Use a mask created from thresholding R (edge intensities)
%to extract the actin's angles from all the angles 
%(edges assumed to be the actin's edges)
thresh = graythresh(R) + 0.1;
actinMask = im2bw(R,thresh);


% Code added to create a color image of fiber angles
% 1/16/2013 Byron
% 1. Select colormap and store as an array
colormap(jet);
cmap = colormap;
% 2. Convert theta to colormap indices
numColors = size(colormap, 1);
thetaCM = min(floor((((theta + 90) / 180) * numColors) + 1), numColors);
% 3. Create 3D array of RGB values from colormap
thetaZero = zeros(size(theta));
thetaRGB = cat(3, thetaZero, thetaZero, thetaZero);
for i = 1:size(theta, 1)
    for j = 1:size(theta, 2)
        clr = cmap(thetaCM(i, j), :);
        thetaRGB(i, j, :) = clr;
    end
end
% 4. Convert mask to RGB
maskRGB = cat(3, actinMask, actinMask, actinMask);
% 5. Apply mask to thetaRGB
thetaRGB = thetaRGB .* maskRGB;
%close all;
%figure, imshow(thetaRGB);
%colorbar;
%figure;
% End of code added to create a color image of fiber angles


%extract the angles of the actin edges using the mask
count = sum(sum(actinMask));            
[nrows ncols] = size(theta);
actinThetas = zeros(1,count); 
index = 1;
for r = 1:nrows
    for c = 1:ncols
        if(actinMask(r,c) == 1)
            actinThetas(index) = theta(r,c);
            index = index + 1;
        end
    end
end



avg = mean(actinThetas);
stddev = std(actinThetas);
%kurt = kurtosis(actinThetas)
%skew = skewness(actinThetas);


% Histogram created using "Scott's normal reference rule" for bin width

ncolsActin = size(actinThetas,2);
h = (3.5*stddev)/(ncolsActin.^(1/3)); %h is bin width
k = (max(actinThetas)-min(actinThetas))/h; %k is the number of bins
numBins = round(k);
% figure(8); hist(actinThetas, round(k)) %Plot the histogram
% title('Histogram of Angle Data')
% xlabel('Bin Center (degrees)')
% ylabel('Number of Pixels')
% figure();
[hist_bin_size binCenters] = hist(actinThetas,round(k));
binSizesNormalized = hist_bin_size/size(actinThetas,2);
% numBins = size(hist_centers,2);
%figure(15); bar(hist_centers,binSizesNormalized); %"normalized hist"


% New Metric does not use curves found by PeakFit software
% 1/17/2013  Byron
%[FitResults LowestError] = automatePeakFit(binCenters,binSizesNormalized);
%FitResults = adjustAnglesOutOfRange(FitResults);

%figure(15); bar(binCenters,binSizesNormalized); %"normalized hist"

% % Images to see what the code is doing

% figure(2);
% T = image(abs(theta));
% title('Abs(theta)')
% map = colormap;

% figure(3);
% imshow(I);
% title('OriginalImage')
% 
% figure(4);
% imshow(R);
% title('Edge Intensities (R)')
% 
% figure(5);
% image(double(actinMask).*theta);
% colormap(map); %gets colormap from figure(2) where theta is imaged
% title('Mask Applied to Angles')
% 
% figure(6);
% imshow(actinMask)
% title('Actin Mask (R-Mask)')

end

function FitResults = adjustAnglesOutOfRange(FitResults)
% Adjust angles outside of [-90 90] degrees range. Only interested in Peak Location (column 'C'=3)

peakColumn =2;
peakLoc = FitResults(:,peakColumn); %column of all peak Locations

for i =1:size(peakLoc,1)
    while(abs(peakLoc(i))>90)
        if(peakLoc(i)>90)
            peakLoc(i) = peakLoc(i)-180;
        else if(peakLoc(i)<-90)
                peakLoc(i) = peakLoc(i)+180;
            end
        end
    end
end

FitResults(:,peakColumn) = peakLoc;

end