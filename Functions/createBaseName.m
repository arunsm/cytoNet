function baseName = createBaseName(imageName)
    k = strfind(imageName, '.');
    
    if isempty(k)
        baseName = imageName;
    else
        baseName = imageName(1:(k(end)-1));
    end
end