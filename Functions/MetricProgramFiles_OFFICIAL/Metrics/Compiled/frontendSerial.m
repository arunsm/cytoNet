function frontendSerial(inputfiledir, outfile)
% frontend(inputfiledir, outfile)

cellmasksuffix = '-MASKS.tif';
nucleusmasksuffix = '-NUMasks.tif';
stainsuffix = '-IMG.tif';
stainorder = {'microtubule', 'actin', 'nucleus'};


cellmaskd = dir(strcat(inputfiledir, '/*', cellmasksuffix));
nucleusmaskd = dir(strcat(inputfiledir, '/*', nucleusmasksuffix));
staind = dir(strcat(inputfiledir, '/*', stainsuffix));

numcellmasks = numel(cellmaskd);
numnucleusmasks = numel(nucleusmaskd);
numstains = numel(staind);

assert(numcellmasks == numnucleusmasks, sprintf('There are %d cell mask files, but %d nucleus mask files', numcellmasks, numnucleusmasks));

assert(numcellmasks == numstains, sprintf('There are %d cell mask files, but %d stain files', numcellmasks, numstains));


cellmaskfiles = cell(numcellmasks, 1);
nucleusmaskfiles = cell(numcellmasks, 1);
stainfiles = cell(numcellmasks, 1);

cellmasksuffixlength = numel(cellmasksuffix);
nucleusmasksuffixlength = numel(nucleusmasksuffix);
stainsuffixlength = numel(stainsuffix);
for i = 1:numcellmasks
    cellmaskname = cellmaskd(i).name;
    cellmaskprefix = cellmaskname(1:end-cellmasksuffixlength);
    nucleusmaskname = nucleusmaskd(i).name;
    nucleusmaskprefix = nucleusmaskname(1:end-nucleusmasksuffixlength);
    stainname = staind(i).name;
    stainprefix = stainname(1:end-stainsuffixlength);
    assert(strcmp(cellmaskprefix, nucleusmaskprefix),...
	   sprintf('Mismatched file names: %s %s', cellmaskname, nucleusmaskname));
    assert(strcmp(cellmaskprefix, stainprefix),...
	   sprintf('Mismatched file names: %s %s', cellmaskname, stainname));
    cellmaskfiles{i} = strcat(inputfiledir, '/', cellmaskname);
    nucleusmaskfiles{i} = strcat(inputfiledir, '/', nucleusmaskname);
    stainfiles{i} = strcat(inputfiledir, '/', stainname);
end


engineSerial(cellmaskfiles, nucleusmaskfiles, stainfiles, stainorder, outfile); 


end
