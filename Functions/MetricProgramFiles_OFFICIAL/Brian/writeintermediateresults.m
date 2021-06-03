function writeintermediateresults(mf, mr, filename)
[r c] = size(mr);
separator = ',';
header_string = '';
for i = 1:r
    fid = fopen(filename, 'a');
    
    lo = 0;
    k = 0;
    
    for j = 1:c
        if isempty(mr{i, j})continue;end
        header = mr{i, j}.header;
        values = mr{i, j}.values;
        [l, k] = size(values);
        if l > lo
            lo = l;
        end
        fprintf(fid, '%d, %d', l+1, k);
        fprintf(fid, '\n')
        iscell(mr{i, j}. values);
        if iscell(header)
            for k = 1:length(header)
                fprintf(fid, '%s,', header{k});
            end
        else
            fprintf(fid, '%s,', header);
        end
        fprintf(fid, '\n');
        if isempty(mr{i, j}) continue; end
        [a b] = size(mr{i, j}.values)
        for p = 1:a
            for m = 1:b
                value = mr{i, j}.values(p, m)
                fprintf(fid, '%f, ', value);
            end
            fprintf(fid, '\n')
        end
    end
    
end
fclose('all')
end
