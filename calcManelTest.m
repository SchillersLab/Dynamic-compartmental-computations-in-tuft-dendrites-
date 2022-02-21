function calcManelTest(activityMatrix, isPearson, distanceMatrix)
    if isPearson
        activityMatrix = 1 - abs(activityMatrix);
    end
    
    [manstat, p1, index, np, r] = mantelexact(activityMatrix, distanceMatrix);
    
end