
function [errorReport, imageProcessingParams] = FrontEnd(InputPath, OutputPath, segmentationStruct, Mask2GraphStruct, sourceCodePath)

tic
addpath(genpath([sourceCodePath, filesep, 'Functions']));

if ~exist(OutputPath, 'file')
    mkdir(OutputPath);
end

% calling engine
[errorReport, imageProcessingParams] = Engine(InputPath, OutputPath, segmentationStruct, Mask2GraphStruct);
writeLog('[FrontEnd] Done');
toc
end
