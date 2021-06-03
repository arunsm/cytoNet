function testengine()
dir('C:\Users\ebionic\Dropbox\MetricProgramFiles_OFFICIAL/TestData')
cellMaskFileNames={'C:\Users\ebionic\Dropbox\MetricProgramFiles_OFFICIAL/TestData/testCellMasks.tif'};
nucleusMaskFileNames={'C:\Users\ebionic\Dropbox\MetricProgramFiles_OFFICIAL/TestData/testNucleusMasks.tif'};
stainFileNames={'C:\Users\ebionic\Dropbox\MetricProgramFiles_OFFICIAL/TestData/testStains.tif'};

outputFileName='../testresults.csv';

stainOrder = {'actin', 'microtubule', 'nucleus'};
%metrics='../test.m';
metrics='IndividualMetrics';

fprintf('Calling engine...\n');
engine(cellMaskFileNames, nucleusMaskFileNames, stainFileNames, stainOrder, metrics, outputFileName);
end

