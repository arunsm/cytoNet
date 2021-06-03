function BorderLength = sharedPixels(allMasksExtended)

nNodes = size(allMasksExtended, 1);

BorderLength = zeros((nNodes*(nNodes - 1))/2, 1);

ctr = 1;

for i = 1:nNodes
    
    mask1 = allMasksExtended{i};
    
    for j = i+1:nNodes
        if i == j
            continue;
        end
        
        mask2 = allMasksExtended{j};
        
        % find overlapping pixels in the extended masks
        border = mask1 & mask2;
        BorderLength(ctr) = sum(border(:));
        
        ctr = ctr + 1;
    end
end

end