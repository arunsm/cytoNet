function [] = createLayeredMasks(MaskFileName)
    Masks = imread(MaskFileName);
    layeredMaskFileName = strrep(MaskFileName, '-MASKS.tif', '-LAYEREDMASKS.tif');
    CC = bwconncomp(Masks);
    
    for i = 1:CC.NumObjects
        mask_i = false(CC.ImageSize(1), CC.ImageSize(2));
        mask_i(CC.PixelIdxList{i}) = true;
        imwrite(mask_i, layeredMaskFileName, 'WriteMode', 'Append');
    end
end