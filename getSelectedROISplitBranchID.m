function selectedROISplit = getSelectedROISplitBranchID(gRoi, depthForSplit, selectedROISplit, selectedROIName, rootNode)
    if depthForSplit == -1
        mainTreeSplite_ID = [];
        alreadyPast = [];
        mainTreeSplite_ID = getFirstNodeWithRoi(rootNode, gRoi, mainTreeSplite_ID, alreadyPast);
        mainTreeSplite = gRoi.Nodes(mainTreeSplite_ID, :);
    else
        mainTreeSplite = gRoi.Nodes(gRoi.Nodes.Depth(: ,1) == depthForSplit, :);
    end

    alreadyPast = [];
    for index = 1:size(mainTreeSplite, 1)
        [selectedROISplit,alreadyPast]  = setMainTreeSplit(gRoi, selectedROISplit,  mainTreeSplite(index, :).ID, selectedROIName, mainTreeSplite(index, :).ID, alreadyPast);
    end 
end

function mainTreeSplite = getFirstNodeWithRoi(rootNode, gRoi, mainTreeSplite, alreadyPast)
    alreadyPast(end + 1) = rootNode;
    
    if contains(gRoi.Nodes(rootNode, :).Name, 'roi')
        mainTreeSplite(end+1) = rootNode;
    else
        nid = neighbors(gRoi,rootNode);
        sumROI = 0;
        sumAll = 0;
        
        for index = 1:length(nid)        
            if (sum (alreadyPast == nid(index)) == 0)
                if contains(gRoi.Nodes(nid(index), :).Name, 'roi')
                  sumROI = sumROI + 1;  
                end
                
                sumAll = sumAll + 1;
            end
        end
        
        if sumROI > 0
            mainTreeSplite(end+1) = rootNode;
        else
            for index = 1:length(nid)        
                if (sum (alreadyPast == nid(index)) == 0)
                    mainTreeSplite = getFirstNodeWithRoi(nid(index), gRoi, mainTreeSplite, alreadyPast);
                end
            end
        end
    end 
end

function [selectedROISplit, alreadyPast] = setMainTreeSplit(gRoi, selectedROISplit, currentID, selectedROIName, rootID, alreadyPast)
    
    alreadyPast(end + 1) = currentID;
    indexSelectedROI = find(contains(selectedROIName, gRoi.Nodes(currentID, :).Name));
    if ~isempty(indexSelectedROI) && selectedROISplit(indexSelectedROI(1)) == -1
        selectedROISplit(indexSelectedROI(1)) = rootID;
    end
    
    nid = neighbors(gRoi,currentID); 
    
    for index = 1:length(nid)        
        if (sum (alreadyPast == nid(index)) == 0) && (gRoi.Nodes(nid(index), :).Depth(1) >= gRoi.Nodes(currentID, :).Depth(1))
            [selectedROISplit, alreadyPast] = setMainTreeSplit(gRoi, selectedROISplit, nid(index), selectedROIName, rootID, alreadyPast);  
        end
    end
end