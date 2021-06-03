function frontend1(inputfiledir, metricsdir, outfile)
% frontend1(inputfiledir, metricsdir, outfile)

cellmasksuffix = '-LAYEREDMASKS.tif';
nucleusmasksuffix = '-LAYEREDMASKS.tif';
stainsuffix = '-IMG.tif';
stainorder = {'microtubule', 'actin', 'nucleus'};

%%% NEW - Arun %%%
% cellmaskd, nucleusmaskd and staind have been set to read only one image
% at a time, instead of all images

% cellmaskd = dir(strcat(inputfiledir, filesep, '*', cellmasksuffix));
% nucleusmaskd = dir(strcat(inputfiledir, filesep, '*', nucleusmasksuffix));
% staind = dir(strcat(inputfiledir, filesep, '*', stainsuffix));

cellmaskd = dir(strcat(inputfiledir, cellmasksuffix));
nucleusmaskd = dir(strcat(inputfiledir, nucleusmasksuffix));
staind = dir(strcat(inputfiledir, stainsuffix));

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
   
   %%% NEW - Arun %%%
%     cellmaskfiles{i} = strcat(inputfiledir, filesep, cellmaskname);
%     nucleusmaskfiles{i} = strcat(inputfiledir, filesep, nucleusmaskname);
%     stainfiles{i} = strcat(inputfiledir, filesep, stainname);
    cellmaskfiles{i} = strcat(inputfiledir, '-LAYEREDMASKS.tif');
    nucleusmaskfiles{i} = strcat(inputfiledir, '-LAYEREDMASKS.tif');
    stainfiles{i} = strcat(inputfiledir, '-IMG.tif');
end


engine(cellmaskfiles, nucleusmaskfiles, stainfiles, stainorder, metricsdir, outfile); 


end
