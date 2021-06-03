function writeresults(header, table, outfilename, badmetric)
separator = ',';
header_string = '';
for i = 1:length(header)
    if i == 1
        header_string = header{i};
    else
        header_string = [header_string, separator, header{i}];
    end
end

[pathstr name ext] = fileparts(outfilename);


% Create the directory if it does not exist
if ~strcmp(pathstr, '') & ~exist(pathstr, 'dir')
    mkdir(pathstr);
end


%write the string to a file

fid = fopen(outfilename,'w');
if nargin > 3
    badmetric_string = '';
    for i = 1:length(badmetric)
        if i == 1
            badmetric_string = sprintf('%d', badmetric(i))
        else
            badmetric_string = [badmetric_string, separator, sprintf('%d', badmetric(i))];
        end
    end
    badmetric_string
    badmetric
fprintf(fid,'%s\r\n',badmetric_string);    
end
fprintf(fid,'%s\r\n',header_string);
fclose(fid);

dlmwrite(outfilename, table, '-append', 'delimiter', separator);


fprintf('Wrote file %s\n', outfilename);
end
            
        
        
        
