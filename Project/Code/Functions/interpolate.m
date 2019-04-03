function x = interpolate(func,value,grid)
    na = length(grid);
    idxLow = max(sum(value>grid),1);
    if idxLow == na
        idxLow = idxLow - 1;
    end
    idxHigh = idxLow + 1;
    while grid(idxHigh) == grid(idxLow)
        if idxHigh < na
            idxHigh = idxHigh + 1;
        else
            idxLow = idxLow - 1;
        end
    end
    
    x = func(idxLow) + (value-grid(idxLow))* (func(idxHigh)-func(idxLow))/(grid(idxHigh)-grid(idxLow));
end