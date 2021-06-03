%
% 1. Metrics are modularized in separate files in a
%    special directory.
% 2. Multiple metrics may be in a single file.
% 3. Metrics return a table of cells.  The top row
%    is a header containing strings naming the Metrics.
%    Each column of data conatins metric results in
%    the order that cell data was passed to the metric.
% 4. Metric function arguments are cell masks and images
% 5. Each metric function must always return a header and
%    values table with the same number of columns
% 6. Stain images contain double values ranging from 0 to 1.


%Need image magnification level as a parameter!!!!!
%Product of magnification and diameter of field of view is constant
%d1m1 = d2m2

%function engine(cellmaskfilenames, nucleusmaskfilenames, cytomaskfilenames, stainfilenames, stainorder, metricsdir, outfilename)

function engine(cellmaskfilenames, nucleusmaskfilenames, stainfilenames, stainorder, metricsdir, outfilename)

%fprintf('Begin engine...\n');

fprintf('Begin stain key construction...\n');
% The stain key indicates the indices of individual stains
% stainkey = struct('nucleus', 0, 'actin', 0, 'microtubule', 0);
% for i = 1:numel(stainorder)
%     switch stainorder{i};
%         case 'nucleus'
%             stainkey.nucleus = i;
%         case 'actin'
%             stainkey.actin = i;
%         case 'microtubule'
%             stainkey.microtubule = i;
%         otherwise
%                 error(sprintf('Unexpected stain name: %s', stainorder{i}));
%     end
% end

% Force fieldnames to occur in the order of their indices
for i = 1:numel(stainorder)
    stainkey.(stainorder{i}) = i;
end


% The data struct holds input data for each metric function
data = struct('cellmaskfilename', '', 'nucleusmaskfilename', '', 'stainfilename', '', 'stainkey', stainkey, 'cellmaskdata', 0, 'nucleusmaskdata', 0, 'staindata', 0);

numcellmaskfiles = numel(cellmaskfilenames);
numnucleusmaskfiles = numel(nucleusmaskfilenames);
numstainfiles = numel(stainfilenames);
assert(numcellmaskfiles == numstainfiles, 'The number of cell mask files, %d, is not the same as the number of stain files, %d\n', numcellmaskfiles, numstainfiles);

fprintf('Begin retrieving metrics...\n');
metricfunctions = getmetricfunctions(metricsdir);
numfunctions = numel(metricfunctions);
badmetric = false(numfunctions, 1); 
% Store metric function output in metric results
metricresults = cell(numcellmaskfiles, numfunctions);

fprintf('Begin metric loop...\n');
% For each image/mask file, run the metrics.  Each file is opened only
% once, and although the contents are buffered, only the contents of one
% file are buffered at any time.
for  i = 1:numcellmaskfiles
    cellmaskfname = cellmaskfilenames{i};
    fprintf('Begin reading cell mask file %s\n', cellmaskfname);
    cellmaskfinfo = imfinfo(cellmaskfname);
    numcellmasks = numel(cellmaskfinfo);
    cellmaskdata = cell(numcellmasks, 1);
    for j = 1:numcellmasks
        cellmaskdata{j, 1} = logical(imread(cellmaskfname, 'Info', cellmaskfinfo, 'Index',j));
    end
    
    nucleusmaskfname = nucleusmaskfilenames{i};
    fprintf('Begin reading nucleus mask file %s\n', nucleusmaskfname);
    nucleusmaskfinfo = imfinfo(nucleusmaskfname);
    numnucleusmasks = numel(nucleusmaskfinfo);
    nucleusmaskdata = cell(numnucleusmasks, 1);
    for j = 1:numnucleusmasks
        nucleusmaskdata{j, 1} = logical(imread(nucleusmaskfname, 'Info', nucleusmaskfinfo, 'Index',j));
    end
    
    %    cytomaskfname = cytomaskfilenames{i};
    %    fprintf('Begin reading cytoplasm mask file %s\n', cytomaskfname);
    %    cytomaskfinfo = imfinfo(cytomaskfname);
    %    numcytomasks = numel(cytomaskfinfo);
    %    cytomaskdata = cell(numcytomasks, 1);
    %    for j = 1:numcytomasks
    %        cytomaskdata{j} = logical(imread(cytomaskfname, 'Info', cytomaskfinfo,% 'Index',j));
    %    end
    
    stainfname = stainfilenames{i};
    fprintf('Begin reading stain file %s\n', stainfname);
    stainfinfo = imfinfo(stainfname);
    numstains = numel(stainfinfo);
    staindata = cell(numstains, 1);
    for j = 1:numstains
        staindata{j} = uimread(stainfname, j);
    end
    
    data.cellmaskfilename = cellmaskfname;
    %    data.cytomaskfilename = cytomaskfname;
    data.nucleusmaskfilename = nucleusmaskfname;
    data.stainfilename = stainfname;
    data.cellmaskdata = cellmaskdata;
    data.nucleusmaskdata = nucleusmaskdata;
    %    data.cytomaskdata = cytomaskdata;
    data.staindata = staindata;
    data.stainkey = stainkey;
    for j = 1:numfunctions
        metric = metricfunctions{j};
        % Don't use metrics that crashed earlier in the loop
        if badmetric(j) || metric.skip continue; end
        f = metric.func;
        fprintf('Invoking function %s...\n', metric.name);
        success = false;
        try
            [header values] = f(data);
            nummetrics = validate(header, values, numcellmasks);
            if metric.number < 0
                metricfunctions{j}.number = nummetrics;
            else
                if nummetrics ~= metricfunctions{j}.number
                    n1 = nummetrics;
                    n2 = metricfunctions{j}.number;
                    throw(Mexception('engine:inconsistent_metric', 'metric returned  %d and %d column results', n1, n2));
                end
            end
            success = true;
        catch e
            % Don't use metric.skip = 1.  It changes a copy of
            % metricfunctions{j}
            metricfunctions{j}.skip = 1;
            badmetric(j) = true
            fprintf('Caught exception in function %s: %s  %s\n', metric.name, e.identifier, e.message);
            stackSize = size(e.stack, 1);
            for s = 1:stackSize
                frame = e.stack(s);
                fprintf('file: %s  function: %s  line: %d\n', frame.file, frame.name, frame.line);
            end
            data
            fprintf('Cell mask file: %s\n', data.cellmaskfilename);
            fprintf('Nucleus mask file: %s\n', data.nucleusmaskfilename);
            fprintf('Stain file: %s\n', data.stainfilename);
            fprintf('\n');
        end
        if success
            % Create the struct without header in case header is an array
            % of cells.
            metricresults{i, j} = struct('header', 0, 'values', values);
            metricresults{i, j}.header = header;
        n = sprintf('%d%s%d', i,'_',j);
        newfilename = newoutfilename(outfilename, n);
        [header values] = assembledata(metricresults, metricfunctions);
        writeintermediateresults(metricfunctions, metricresults, newfilename);
%         writeresults(header, values, newfilename, badmetric);
        if i~= 1 || j~=1
        delete(lastfile);
        end
        lastfile = newfilename;

        end
    end
end

[header data] = assembledata(metricresults, metricfunctions);
writeresults(header, data, outfilename);
delete(newoutfilename(outfilename, n));
end
%cd(startingLocation);




