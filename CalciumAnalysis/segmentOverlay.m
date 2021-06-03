%Input: Name of file to be segmented.
%segmentOverlay takes a grayscale image and segments it the segmentNPC method.
%It then overlays the boundaries of the segmentation on the original
%grayscale.
function thresholded = segmentOverlay(bw, med)
edges = bwmorph(bw,'remove');
i = med;
R = i;
G = i;
B = i;
B(bw) = 1;
B(edges) = 0;
G(edges) = 0;
R(edges) = 0;
thresholded = cat(3, R, G, B);
end
