function c = CorrCoeff(sig1,sig2,std1, std2, mn1, mn2)
N1 = length(sig1);
crosscov = 0;
for i = 1:N1
    current = 0;
    for j = (i):N1
        current = current + (sig1(j - i + 1) - mn1)*(sig2(j) - mn2);
    end
    current = current / (N1 - 1);
    crosscov = max(crosscov, abs(current));
end
for i = 1:N1
    current = 0;
    for j = (i):N1
        current = current + (sig2(j - i + 1) - mn2)*(sig1(j) - mn1);
    end
    current = current / (N1 - 1);
    crosscov = max(crosscov, abs(current));
end
c = crosscov / (std1*std2);
end