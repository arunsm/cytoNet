function fileName = getFileName(path)
k = strfind(path, filesep);
if isempty(k)
    fileName = path;
else
    fileName = path((k(end) + 1):end);
end
end
