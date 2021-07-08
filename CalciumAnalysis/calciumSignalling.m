function statusArr = calciumSignalling(inputDirName, outputDirName, sourceCodePath)
%addpath(genpath(sourceCodePath)); % add source code directory to path?
addpath(genpath([sourceCodePath filesep 'Functions']));

statusArr = repmat(makeErrorStruct('', ''), 0);

d = dir(inputDirName);
maxDigits = countDigits(numel(d));
folderNum = 1;
for i = 1:numel(d)
    if d(i).isdir
        continue;
    end
    fileName = d(i).name;
    status = validateImageFile(inputDirName, fileName);
    if ~isempty(status)
        statusArr(end+1) = status;
    else
        outputSubDirName = [outputDirName, filesep, pad(folderNum, maxDigits)];
        folderNum = folderNum + 1;
        [status,msg] = mkdir(outputSubDirName);
        if status ~= 1
            error('[calciumSignalling] Unable to create directory %s', outputSubDirName);
        end
        fid = fopen([outputDirName, filesep, 'log.txt'], 'a');
        fprintf(fid, '[calciumSignalling] Calling calciumEngine; fileName=%s\n', fileName);
        fclose(fid);
        errorReport = calciumEngine([inputDirName, filesep, fileName], outputSubDirName); % AM - calling calciumEngine; send optional mask here?
        if ~isempty(errorReport)
            statusArr(end+1) = errorReport;
        end
    end
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
