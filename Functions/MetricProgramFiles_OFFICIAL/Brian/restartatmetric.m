function restartatmetric (filename)
fnum = 1;
i = 1;
for i = 1:9
        num = sprintf('%d', i);
        a = strfind(filename, num);
        if isequal(a, []) == 0
            a = a(1);
        end
end

b =strfind(filename, '_');
c = strfind(filename, '.');
filenum = sprintf('%s', filename(a:b-1));
if isequal(c, [])
    metricnum = sprintf('%s', filename(b+1:length(filename)));
else
    metricnum = sprintf('%s', filename(b+1:c-1));
end
fprintf('Resuming at file %s and metric %s. \n', filenum, metricnum);
end