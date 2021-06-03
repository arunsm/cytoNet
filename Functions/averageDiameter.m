function d2 = averageDiameter(d1, dRest)
    
    d2 = zeros(numel(dRest), 1);
    for i = 1:numel(dRest)
        d2(i) = (d1 + dRest(i))/2;
    end

end