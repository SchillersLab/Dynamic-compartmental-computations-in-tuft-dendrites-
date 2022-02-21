function matrixForSP = getMatrixForDistanceSP(gRoi)
    matrixForSP = zeros(size(gRoi.Nodes, 1));
    
    for nid = 1:size(gRoi.Nodes, 1)
        for nid_sec = 1:nid
            [~, d] = shortestpath(gRoi, gRoi.Nodes.ID(nid), gRoi.Nodes.ID(nid_sec));
            matrixForSP(nid, nid_sec) = d;
        end
    end
end