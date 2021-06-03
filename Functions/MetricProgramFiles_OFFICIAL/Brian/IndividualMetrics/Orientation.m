function [header values] = Orientation(data)

fieldnames(data);
names = fieldnames(data.stainkey);

masks = data.cellmaskdata;
nummasks = numel(masks);
values = zeros(nummasks, 1);
m = data.cellmaskdata{1};

count2 = 1;
nucleus = data.stainkey.nucleus;
for i = 1:length(data.staindata)
    if i == nucleus
        continue
    end
    H1 = fspecial('sobel'); %returns the matrix that is the sobel filter
    % (horizontal, the transpose is the vertical)
    I = data.staindata{i};
    Gh = imfilter(I,H1); %apply horizontal sobel filter
    Gv = imfilter(I,H1'); %apply vertical sobel filter
    
    %Calculate angles (theta) and edge intensities (R)
    theta = atan2(Gv,Gh).*(180/pi); %REVERSED ratio from paper, returns degrees
    neg = find(theta<0);
    theta(neg) = theta(neg) + 360;
    great = find(theta>=180);
    theta(great) = theta(great) - 180;
    
     
    R = (Gh.^2 + Gv.^2).^0.5;
    thresh = graythresh(R) + 0.1;
    edgeMask = im2bw(R,thresh);
    
    d = 10; % Sets the range of angles
    
    for j = 1:nummasks
        m = data.cellmaskdata{j};
        stats = regionprops(m, 'Orientation');
        cellOrientation = stats.Orientation;

        cc = bwconncomp(edgeMask & m);
        count = 0;
        edge = 0;
        range = cellOrientation;
        for k = 1:length(cc.PixelIdxList)
            idx = cc.PixelIdxList{k};
            theta2 = theta(idx);
            for l = 1:length(theta2)
                if abs(range - theta2(l)) >= 90
                    diff = 180 - abs(range - theta2(l));
                else
                    diff = abs(range - theta2(l));
                end
                
                if diff <= d;
                    count = count + 1;
                end
                
            end
            edge = edge + length(idx);
        end
        
        fraction = count/edge;

        head = names{i, 1};
        header{1, count2} = sprintf('%s alignment',head);
        values(j, count2) = fraction;
        
    end
  count2 = count2 + 1;
end

end


