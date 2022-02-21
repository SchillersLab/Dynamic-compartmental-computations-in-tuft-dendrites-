function [matrixForHS, indexTemp, matrixForSP] = getMatrixForDistanceHS_Type3(gRoi, rootNodeID)    
    tempId = [];
    MainNodeTemp = zeros(1, size(gRoi.Nodes, 1));
    
    for nid = 1:size(gRoi.Nodes, 1)
        [p, d] = shortestpath(gRoi, gRoi.Nodes.ID(nid), rootNodeID);

        for indexP = 1:length(p)
            current = find(gRoi.Nodes.ID == p(indexP));
            
            if nid ~= current && contains(gRoi.Nodes.Name(nid), 'roi') && contains(gRoi.Nodes.Name(current), 'roi') && gRoi.Nodes.Depth(nid, 1) == gRoi.Nodes.Depth(current, 1) 
                nextP = find(gRoi.Nodes.ID == p(indexP +  1));
                if contains(gRoi.Nodes.Name(current), 'roi') && ~contains(gRoi.Nodes.Name(nextP), 'roi')
                    MainNodeTemp(current) = nextP;
                    tempId(end+1) = current;
                end
            end           
        end       
    end
    
    tempU = unique(tempId);
    indexTemp = length(tempU);
    
    for tempI = 1:length(tempU)
        index_e = findedge(gRoi,tempU(tempI),MainNodeTemp(tempU(tempI)));
        w_old = gRoi.Edges.Weight(index_e);
        gRoi = rmedge(gRoi,index_e);
        gRoi = addnode(gRoi, sprintf('temp%02d', tempI));
        
        current = find(strcmp(gRoi.Nodes.Name,sprintf('temp%02d', tempI)));
            
        gRoi.Nodes.ID(current) = current;
        gRoi = addedge(gRoi,tempU(tempI),current, 1);
        gRoi = addedge(gRoi,MainNodeTemp(tempU(tempI)),current, w_old);        
    end
    
    matrixForHS = getMatrixForDistanceHS(gRoi, rootNodeID);
    matrixForSP = getMatrixForDistanceSP(gRoi);
end



