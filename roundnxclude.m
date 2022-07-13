function  [finalperc,finalcurrent] = roundnxclude(optcurrent,minperc)

totalcurrent = sum(optcurrent);
percurrent = optcurrent*100./totalcurrent;
    
finalcurrent = round(totalcurrent,1);
smallperc = sum(percurrent < minperc);

% Exclude contacts < options.optimizer.general.minperc
while smallperc >0
    [mincurrent,contind] = nanmin(percurrent);
    if mincurrent < minperc
        percurrent(contind) = nan;
    end
    percurrent = percurrent*100/nansum(percurrent);
    smallperc = sum(percurrent < minperc);
end
percurrent(isnan(percurrent)) = 0;

% Round to 100%
persum = sum(floor(percurrent));
[~,rankedmods] = sort(mod(percurrent,1),'descend');
j = 1;
while persum <100
    percurrent(rankedmods(j)) = percurrent(rankedmods(j))+1;
    persum = sum(floor(percurrent));
    j = j+1;
end
finalperc = floor(percurrent);

