function metricfunctions = getmetricfunctions(metricsdir)

% If metricsdir is a single file, then get its function
if exist(metricsdir, 'file') & ~exist(metricsdir, 'dir')
   [pathstr name ext] = fileparts(metricsdir);
   addpath(pathstr);
   metricfunctions{1} = struct('name', name, 'func', str2func(name), 'skip', 0, 'number', -1);
else
   % Otherwise get the functions in the .m files
   addpath(metricsdir);
   dirContents = dir(strcat(metricsdir, '/*.m'));
   dirSize = size(dirContents, 1);
   metricfunctions = cell(dirSize, 1);
   for i = 1:dirSize
      [pathstr name ext] = fileparts(dirContents(i).name);
      fileName = dirContents(i).name;
      functionName = fileName(1:end-2);
      metricfunctions{i} = struct('name', name, 'func', str2func(name), 'skip', 0, 'number', -1);
   end
end
    

