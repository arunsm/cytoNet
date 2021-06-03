function s = newoutfilename(filename, n)
m = findstr('.', filename);
if length(m) == 0
    s = sprintf('%s', filename, n);
else
    pos = m(length(m));
    s = sprintf('%s', filename(1:pos-1), n ,filename(pos:length(filename)));
end
end



