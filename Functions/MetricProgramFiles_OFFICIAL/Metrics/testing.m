    data.cellmaskdata = true(8);

    H1 = fspecial('sobel'); %returns the matrix that is the sobel filter
    % (horizontal, the transpose is the vertical)
    I = A;
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
    
    d = 10;
    
        m = data.cellmaskdata;
        stats = regionprops(m, 'Orientation');
        cellOrientation = stats.Orientation;
        disp(cellOrientation);
        
        cc = bwconncomp(edgeMask & m);
        count = 0;
        edge = 0;
        range = cellOrientation + d;
        for k = 1:length(cc.PixelIdxList)
            idx = cc.PixelIdxList{k};
            theta2 = theta(idx);
            for l = 1:length(theta2)
                if abs(range - theta2(l)) >= 0
                    diff = 180 - abs(range - theta2(l));
                else
                    diff = abs(range - theta2(l));
                end
                
                if diff <= cellOrientation + d;
                    count = count + 1;
                end
                
            end
            edge = edge + length(idx);
        end
        
        fraction = count/edge;
