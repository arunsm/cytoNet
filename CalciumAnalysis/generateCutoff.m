% Generates a cuttoff correlation value by scrambling each data set, finding
% all correlation coefficients, and then using the 99th percentile of the
% set of coefficients as the cutoff value
function cutoff = generateCutoff(timeSeries, cutoffPercentile, M)
scrambledTimeSeries = timeSeries;
s = size(scrambledTimeSeries);
nCells = s(2);
nFrames = s(1);

% scramble data, while preserving total activity
for j = 1:nCells
    t = randi(nFrames);
    front = scrambledTimeSeries(1:t-1, j);
    back = scrambledTimeSeries(t:nFrames, j);
    index = nFrames - t;
    scrambledTimeSeries(1:(index+1), j) = back;
    scrambledTimeSeries((index+2):nFrames, j) = front;
end

% generate vector of correlation coefficients
crossCorrelation = [];
parfor j = 1:nCells
    scrambledTimeSeries1 = scrambledTimeSeries(:, j);
    for k = (j+1):nCells
        scrambledTimeSeries2 = scrambledTimeSeries(:, k);
        if j ~= k
            crossCorrelation = [crossCorrelation max(xcov(scrambledTimeSeries1, scrambledTimeSeries2, M, 'coeff'))];
        end
    end
end

% set cutoff
cutoff = prctile(crossCorrelation, cutoffPercentile);
end