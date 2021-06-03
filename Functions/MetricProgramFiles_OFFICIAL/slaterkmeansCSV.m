% Assume that the analysis will be for cluster sizes ranging from k1 to k2
% and that there are m metrics and n cells.  The upper table of cell
% distances to cluster centers will have n+1 rows and sum(k1, k2)+1
% columns. The lower table of cluster center coordinates will have
% sum(k1, k2)+1 rows and m+1 columns.  Assume that there are r blank rows
% between the upper and lower tables.  The entire data foorprint will be
% n+1+r+sum(k1, k2)+1 rows and max(sum(k1, k2)+1, m+1) columns.

function slaterkmeans(filename)
% Number of repetitions for kmeans clustering
repetitions = 10000;
% Range of k values for kmeans clustering
firstk = 2;
lastk = 10;
% Number of empty rows between upper and lower tables
tablegap = 5;

% Create output file name
dotIndices = strfind(filename, '.');
suffix = 'kmeans.csv';
if length(dotIndices) == 0
    outfilename = strcat(filename, suffix);
else
    nameEnd = dotIndices(end) - 1;
    outfilename = strcat(filename(1:nameEnd), suffix);
end
if exist(outfilename, 'file')
    delete(outfilename);
    fprintf('Deleted file %s\n', outfilename);
end

[data varnames casenames] = tblread(filename, '\t');
numcells = length(casenames);
nummetrics = length(varnames);

sumk = sum(firstk:lastk);
numtablerows = numcells + 1 + tablegap + sumk + 1;
numtablecols = max(sumk+1, nummetrics+1);
table = cell(numtablerows, numtablecols);

firstdatarow = 2;
col = 1;

% Write cell names in first column
table{1, col} = 'Cell Name';
cellnamerow = firstdatarow;
for i = 1:numcells
    table{callnamerow, col} = strtrim(casenames(i, :));
    cellnamerow = cellnamerow + 1;
end

% Write 'Center Coordinates' label below cell names to indicate where
% cluster center coordinates appear
clusterheaderrow = numcells + tablegap;
table{clusterheadrow, col} = 'Center Coordinates';


% Write individual metric name labels on same row as 'Center Coordinates'
% label
metricnames = cell(1,nummetrics);
metriccol = col + 1;
lastmetriccol = colstr;
for i = 1:nummetrics
    table{clusterheadrow, metriccol} = strtrim(varnames(i, :));
    metriccol = metriccol + 1;
end


for k = firstk:lastk
    fprintf('Starting kmeans clustering (k=%d)...\n', k);
    [IDX C sumd D] = kmeans(data, k, 'Replicates', repetitions);
    [sil, h] = silhouette(data, IDX);
    avg = mean(sil);
    col = col + 1;
    % Write K (number of clusters) label on first row 
    table{1, col} = sprintf('k=%d cluster; avg sil=%d', k, avg);
    % Write cluster number of each cell in a single column
    for i = 1:numel(IDX)
        table{firstdatarow + i - 1, col} = IDX(i);
    end
    % Use k columns to write distance of each cell to each cluster center 
    for i = 1:k
        col = col + 1;
        % Write label
        table{1, col} = sprintf('Dist to cluster %d', i);
        % Write distance to i-th cluster center in a single column
        for j = 1:size(D, 1)
            table{firstdatarow + i - 1, col} = D(j, i);
        end
    end
    % Use k rows below the individual cell data table to list the
    % coordinates of each of the k cluster centers
    clusterdatarow = clusterheaderrow + 1;
    for i = 1:k
        table{clusterdatarow, 1} = sprintf('k=%d', k);
        for j = 1:nummetrics
           table{clusterdatarow, j + 1} = C(i, j);
        end
        clusterdatarow = clusterdatarow + 1;
    end
        
end



end


function s2 = inc(s)
numchars = length(s);
s2 = char(1:numchars);
carry = 1;
for i = numchars:-1:1
    c = char(s(i) + carry);
    if (c > 'Z')
        c = 'A';
        carry = 1;
    else
        carry = 0;
    end
    s2(i) = c;
end
if carry == 1
    s2 = strcat('A', s2);
end
end
