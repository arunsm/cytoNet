% inputDirName to have two sub-directories - MASKS and IMAGES; each will
% have matching file base names to indicate corresponding image and mask
% combinations

function statusArr = calciumSignalling(inputDirName, outputDirName, sourceCodePath)

s = makeErrorStruct('', '');
errorReport = repmat(s, 0);

% add source code directory to path
addpath(genpath(sourceCodePath));

statusArr = repmat(makeErrorStruct('', ''), 0);

dirImages = dir(strcat(inputDirName, filesep, 'IMAGES'));
maxDigits = countDigits(numel(dirImages));
folderNum = 1;
for i = 1:numel(dirImages)
    if dirImages(i).isdir
        continue;
    end
    
    % enter image files
    fileName = dirImages(i).name;
    imageFilePath = strcat(inputDirName, filesep, 'IMAGES', filesep, fileName);
    status = validateImageFile(strcat(inputDirName, filesep, 'IMAGES'), fileName);
    fileNameBase = createBaseName(fileName);
    
    % search for mask file with same base name as image file
    dirMasks = dir(strcat(inputDirName, filesep, 'MASKS', filesep, fileNameBase, '.*'));
    if numel(dirMasks) == 1
        maskPath = [inputDirName, filesep, 'MASKS', filesep, dirMasks.name];
        maskAvailable = 1;
        try
            mask = imread(maskPath);
        catch
            errorReport(end+1) = makeErrorStruct(['mask file ', dirMasks.name, ' is not a readable image file; using segmented maximum intensity image as mask'], 2);
            maskAvailable = 0;
        end
    elseif ~islogical(mask)
        errorReport(end+1) = makeErrorStruct(['mask file ', dirMasks.name, ' is not a binary image; using segmented maximum intensity image as mask'], 2);
    else
        maskAvailable = 0;
    end
    
    if ~isempty(status)
        statusArr(end+1) = status;
    else
        outputSubDirName = [outputDirName, filesep, pad(folderNum, maxDigits)];
        folderNum = folderNum + 1;
        [status, ~] = mkdir(outputSubDirName);
        if status ~= 1
            error('[calciumSignalling] Unable to create directory %s', outputSubDirName);
        end

        writeLog(sprintf('[calciumSignalling] calling calciumEngine; fileName=%s\n', fileName));
        
        if maskAvailable
            errorReport = calciumEngine(errorReport, imageFilePath, outputSubDirName, maskPath);
        else
            errorReport = calciumEngine(errorReport, imageFilePath, outputSubDirName);
        end
        if ~isempty(errorReport)
            statusArr(end+1) = errorReport;
        end
    end
    writeLog('[calciumSignalling] done');
end
end

function numDigits = countDigits(n)
numDigits = 1;
while n >= 10
    n = n / 10;
    numDigits = numDigits + 1;
end
end

function s = pad(n,len)
s = sprintf('%d', n);
while numel(s) < len
    s = ['0', s];
end
end


function status = validateImageFile(inputDirName, fileName)
status = [];
filePath = [inputDirName, filesep, fileName];
try
    info = imfinfo(filePath);
catch e
    status = makeErrorStruct(sprintf('%s is not a readable image file', fileName), 2);
    return;
end
if numel(info) < 3
    status = makeErrorStruct(['file ', fileName, ' does not contain at least 3 images'], 1);
    return;
end
width = info(1).Width;
height = info(1).Height;
for i = 2:numel(info);
    if width ~= info(i).Width || height ~= info(i).Height
        status = makeErrorStruct(sprintf('%s does not contain images of a constant size', fileName), 2);
        return;
    end
end
end
