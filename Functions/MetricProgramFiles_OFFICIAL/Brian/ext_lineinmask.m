function [inX inY outX outY polyX polyY] = ext_lineinmask(varargin)
%LINEINMASK determines portions of the line in and out of a mask
%
%   [inX inY] = LINEINMASK(XS,YS,BW,XV,YV) computes the segments of the
%       line specified by verticies XV and YV that are inside the mask BW.
%       Points XS and YS specify the limits (or pixel values) of the
%       coordinates of the matrix BW.
%
%   [inX inY] = LINEINMASK(BW,XV,YV) assumes XS=1:N and YS=1:M where
%       [M,N]=SIZE(BW).
%
%   [inX inY] = LINEINMASK(...,METHOD) Uses the given method to create the
%       polygon from the mask.
%           Inner | Outer | {Exact}
%
%   [inX inY outX outY] = LINEINMASK(...) outputs the segments of the line
%       both in and out of the mask.
%
%   [inX inY outX outY polyX polyY] = LINEINMASK(...) outputs the polygon
%       describing the mask.
%
%   Example:
%       I = imread('rice.png');
%       level = graythresh(I);
%       bw = im2bw(I,level);
%       bw = bwareaopen(bw, 50);
%       [a b] = size(bw);
%       xv = b*rand(10,1);
%       yv = a*rand(10,1);
%       [inX inY outX outY] = lineinmask(bw,xv,yv,'exact');
%       figure
%       imagesc(bw)
%       colormap(bone)
%       hold on;
%       h1 = plot(inX,inY,'r');
%       h2 = plot(outX,outY,'b');
%       legend([h1 h2],'Inside Mask','Outside Mask');
% 
% By J Sullivan, August 2011

% Parse Inputs
[xs ys bw xv yv method] = parseInputs(varargin{:});

% Turn the mask into a polygon
c = getthepoly(xs,ys,bw,method);

% Initialize the various lines
x1 = xv(1);
y1 = yv(1);
x2 = NaN;
y2 = NaN;
flag = false;

% Parse the line segments
for ii = 2:length(xv)
    % Check to see if the line intersects the polygon mask
    pnts = InterX(xv(ii-1:ii),yv(ii-1:ii),c(1,:),c(2,:));
    
    % If it doesn't intersect
    if isempty(pnts);
        if flag
            x2 = [x2 xv(ii)];
            y2 = [y2 yv(ii)];
        else
            x1 = [x1 xv(ii)];
            y1 = [y1 yv(ii)];
        end
    else
        % Put them in the correct order
        if xv(ii) ~= xv(ii-1)
            [~, inds] = sort([xv(ii-1); pnts(:,1); xv(ii)]);
        else
            [~, inds] = sort([yv(ii-1); pnts(:,1); yv(ii)]);
        end
        if inds(1) ~= 1;
            inds = flipud(inds);
        end
        pnts = pnts(inds(2:end-1)-1,:);
        
        % Add the cross points and NaNs
        for jj = 1:size(pnts,1)
            if flag
                x2 = [x2 pnts(jj,1) NaN];
                y2 = [y2 pnts(jj,2) NaN];
                x1 = [x1 pnts(jj,1)];
                y1 = [y1 pnts(jj,2)];
            else
                x1 = [x1 pnts(jj,1) NaN];
                y1 = [y1 pnts(jj,2) NaN];
                x2 = [x2 pnts(jj,1)];
                y2 = [y2 pnts(jj,2)];
            end
            flag = ~flag;
        end
        
        % Add the end point
        if flag
            x2 = [x2 xv(ii)];
            y2 = [y2 yv(ii)];
        else
            x1 = [x1 xv(ii)];
            y1 = [y1 yv(ii)];
        end
    end
end

if nargout >= 5
    polyX = c(1,:);
    polyY = c(2,:);
end

% Store it as in/out
if interp2(xs,ys,bw,xv(1),yv(1),'nearest',0)
    inX = x1;
    inY = y1;
    outX = x2;
    outY = y2;
else
    inX = x2;
    inY = y2;
    outX = x1;
    outY = y1;
end
end


%Parse Inputs
function [xs ys bw xv yv method] = parseInputs(varargin)
if nargin < 5
    bw = varargin{1};
    xv = varargin{2};
    yv = varargin{3};
    if nargin == 4;
        method = varargin{4};
    else
        method = 'exact';
    end
    [a b] = size(bw);
    xs = 0:b+1;
    ys = 0:a+1;
elseif nargin >= 5
    xs = varargin{1};
    ys = varargin{2};
    bw = varargin{3};
    xv = varargin{4};
    yv = varargin{5};
    if nargin == 6;
        method = varargin{6};
    else
        method = 'exact';
    end
    [a b] = size(bw);
    if length(xs) == 2 && length(xs) ~= size(bw,2)
        dx = diff(xs)/b;
        dy = diff(ys)/a;
        xs = linspace(xs(1)-dx,xs(2)+dx,b+2);
        ys = linspace(ys(1)-dy,ys(2)+dy,a+2);
    end
end
bw = [false(a+2,1) [false(1,b); bw; false(1,b)] false(a+2,1)];
bw = 1*logical(bw);
end

%Get the polygon
function c = getthepoly(xs,ys,bw,method)
switch lower(method(1))
    case 'e'
        c = contourc(xs,ys,double(bw),[0.5 0.5]);
    case 'o'
        c = contourc(xs,ys,double(bw),[0 0]);
    case 'i'
        c = contourc(xs,ys,double(bw),[1 1]);
end
lc = length(c);
ind = 1;

while ind < lc
    indNext = c(2,ind)+ind+1;
    c(:,ind) = NaN;
    ind = indNext;
end
end

function [xi yi] = InterX(x1,y1,x2,y2)
%INTERX Intersection of curves
%   P = INTERX(x1,y1,x2,y2) returns the intersection points of two curves.
%   The curves can be either closed or open and are described by the points
%   (x1,y1) or (x2,y2), where each variable is a 1D vector of points. The
%   intersection of groups of curves (e.g. contour lines, multiply
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L1 = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(x1,y1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(x1,y1,x1,y1) returns all the points of the curve
%   together with any self-intersection points.
%
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX(x1,y1,x2,y2);
%       plot(x1,y1,x2,y2,P(:,1),P(:,2),'ro')

%   Author : NS
%   Tweaked: J Sullivan
%   Version: 3.1, 21 Jan. 2011

%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point.
%   Each factor of the 'C' arrays is essentially a matrix containing
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.

%...Argument checks and assignment of L2
error(nargchk(2,4,nargin));
if nargin == 2,
    x2 = x1;
    y2 = y1;
    hF = @lt;   %...Avoid the inclusion of common points
else
    hF = @le;
end

%...Preliminary stuff
x1 = x1(:);
y1 = y1(:);
x2 = x2(:)';
y2 = y2(:)';

dx1 = diff(x1); dy1 = diff(y1);
dx2 = diff(x2); dy2 = diff(y2);

%...Determine 'signed distances'
S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);

C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

%...Obtain the segments where an intersection is expected
[i,j] = find(C1 & C2);
if isempty(i)
    P = zeros(2,0);
    if nargout > 1
        xi = P(:,1);
        yi = P(:,2);
    else
        xi = P;
    end
    return;
end;

%...Transpose and prepare for output
i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0

%...Solve system of eqs to get the common points
P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
    dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows');

if nargout > 1
    xi = P(:,1);
    yi = P(:,2);
else
    xi = P;
end

    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end