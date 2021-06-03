function totalVariableSize = returnTotalVariableSize(a)
    VariableSize = cat(1, a.bytes);
    totalVariableSize = sum(VariableSize)/10^9; % in gigabytes
end