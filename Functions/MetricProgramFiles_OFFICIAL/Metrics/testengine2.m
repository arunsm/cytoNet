function testengine2()
cellMaskFileNames={'../TestData/testCellMasks.tif'};
nucleusMaskFileNames={'../TestData/testNucleusMasks.tif'};
stainFileNames={'../TestData/testStains.tif'};

outputFileName='../testresults.csv';

stainOrder = {'actin', 'microtubule', 'nucleus'};
%metrics='../test.m';
metrics='../Orientation.m';

fprintf('Calling engine...\n');
engine(cellMaskFileNames, nucleusMaskFileNames, stainFileNames, stainOrder, metrics, outputFileName);
end

