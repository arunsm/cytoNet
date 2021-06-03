% Returns the number of columns in header/data 
function numcols = validate(header, data, numrows)
numcols = -1;
id = 'validate:invalid_metric_output';

% Check that header is either a string or a cell array
if ~ischar(header) & ~iscell(header)
    throw(MException(id, 'Header is neither a string nor a cell array of strings'));
end
[hr hc] = size(header);

% Check that header has only one row
if hr ~= 1
    throw(MException(id, 'Header has %d rows', hr));
end

% If header is a cell array, check that all elements are strings
if iscell(header)
    for i = 1:hc
        if ~ischar(header{i})
            throw(MException(id, 'Header contains non-string values'));
        end
    end
else
    % If header is a single string, then there should be only 1 column of
    % data
    hc = 1;
end

% Check that data contains numbers
if ~isfloat(data) || sum(isinf(data(:))) ~= 0 || sum(isnan(data(:))) ~= 0
    throw(MException(id, 'Data is not an array of numbers or contains NaN or inf values'));
end
[dr dc] = size(data);


% Check that header and data have the same number of columns
if dc ~= hc
    throw(MException(id, 'Header (%d) and data (%d) do not have the same number of columns', hc, dc));
end

% Check that data has the expected number of rows
if dr ~= numrows
    throw(MException(id, 'Data table has %d rows; expected %d rows', dr, numrows));
end

numcols = dc;
end
                
            