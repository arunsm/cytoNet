function writeLog(s)
fid = fopen('log.txt', 'a');
fprintf(fid, '%s\n', s);
fclose(fid);
end
