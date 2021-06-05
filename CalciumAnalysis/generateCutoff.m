%Generates a cuttoff correlation value by scrambling each data set, finding
%all correlation coefficients, and then using the 99th percentile of the
%set of coefficients as the cutoff value
function cutoff = generateCutoff(cellinfo, percentile, M)
scrambled_cellinfo = cellinfo;
s = size(scrambled_cellinfo);

%Scrambles data, preserving total activity 
for j = 1:s(1)
    t = randi(s(3));
    front = scrambled_cellinfo(j,3,1:t-1);
    back = scrambled_cellinfo(j,3,t:s(3));
    index = s(3) - t;
    scrambled_cellinfo(j, 3, 1:index + 1) = back;
    scrambled_cellinfo(j, 3, index + 2:s(3)) = front;
end
crossCorrelation = [];
%sigdevs = std(cellinfo(:,3,:), 0 , 3);
%sigmeans = mean(cellinfo(:,3,:), 3);
%Generate vector of Correlation Coefficients
parfor j = 1:s(1)
    for k = j+1:s(1)
        if j ~= k
            crossCorrelation = [crossCorrelation max(xcov(squeeze(scrambled_cellinfo(j,3,:)), squeeze(scrambled_cellinfo(k,3,:)), M, 'coeff'))];
        end
    end
end
%Makes cutoff 99th percentile of the coefficients
cutoff = prctile(crossCorrelation, percentile);
end