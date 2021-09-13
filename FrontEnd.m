
function [errorReport, imageProcessingParams] = FrontEnd(inputPath, outputPath, segmentationStruct, Mask2GraphStruct, sourceCodePath)

tic
addpath(genpath([sourceCodePath, filesep, 'Functions']));

if ~exist(outputPath, 'file')
    mkdir(outputPath);
end

% calling engine
[errorReport, imageProcessingParams] = Engine(inputPath, outputPath, segmentationStruct, Mask2GraphStruct);
writeLog('[FrontEnd] Done');
toc
end
