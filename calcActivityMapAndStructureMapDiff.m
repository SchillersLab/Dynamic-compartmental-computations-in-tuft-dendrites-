function diffMatrix = calcActivityMapAndStructureMapDiff(roiActivityDistanceMatrix, roiTreeDistanceMatrix, isPearsonActivity)
    if isPearsonActivity
        roiActivityDistanceMatrix = 1 - abs(roiActivityDistanceMatrix);
    else
        roiActivityDistanceMatrix = roiActivityDistanceMatrix ./ max(roiActivityDistanceMatrix, [], 'all');
    end
    
    norm_structure = (roiTreeDistanceMatrix ./ max(roiTreeDistanceMatrix, [], 'all'));
    
    diffMatrix = abs(roiActivityDistanceMatrix - norm_structure);     
end