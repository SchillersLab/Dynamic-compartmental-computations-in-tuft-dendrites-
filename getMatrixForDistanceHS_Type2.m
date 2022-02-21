function matrixForHS = getMatrixForDistanceHS_Type2(gRoi, rootNodeID)
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
                
                if contains(gRoi.Nodes.Name(nid), 'roi') && contains(gRoi.Nodes.Name(current), 'roi') && gRoi.Nodes.Depth(nid, 1) == gRoi.Nodes.Depth(current, 1) 
                    nextP = find(gRoi.Nodes.ID == p(indexP + 1));
                    if contains(gRoi.Nodes.Name(current), 'roi') && ~contains(gRoi.Nodes.Name(nextP), 'roi')
                        matrixForHS(nid, current) = w_sum;
                    else
                        matrixForHS(nid, current) = 0;
                    end
                else
                    matrixForHS(nid, current) = w_sum;
                end
                
                prevIndex = current;
            end
        end
    end
end

