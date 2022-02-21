function matrixForHS = getMatrixForDistanceHS_NoFix(gRoi, rootNodeID)
    matrixForHS = zeros(size(gRoi.Nodes, 1));
    
    for nid = 1:size(gRoi.Nodes, 1)
        [p, d] = shortestpath(gRoi, gRoi.Nodes.ID(nid), rootNodeID);
        
        w_sum = 0;
        prevIndex = 0;
        for indexP = 1:length(p)
            current = find(gRoi.Nodes.ID == p(indexP));
            
            if (current == nid)
                matrixForHS(nid, current) = 1;
                prevIndex = current;
            else                                               
                index_e = findedge(gRoi,p(indexP),gRoi.Nodes.ID(prevIndex));
                w_sum = w_sum + gRoi.Edges.Weight(index_e);
                
                matrixForHS(nid, current) = w_sum;
                
                prevIndex = current;
            end
        end
    end
end
