function newlimits = equalaxes(currenlimits)

lo = min(currenlimits(1:2:end));
hi = max(currenlimits(2:2:end));
newlimits = [lo,hi,lo,hi];

end