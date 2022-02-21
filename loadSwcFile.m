function [gRoi, rootNodeID, selectedROI] = loadSwcFile(neuronTreeFile, outputpath, doComboForCloseRoi)
    data = importdata(neuronTreeFile);
    rootNodeID = [];

    % create nodes
    edgeIndex = 1;
    lengthArrays = (size(data.textdata,1)) - 1;
    nodesLabel = cell(1, lengthArrays);
    nodesX = zeros(1, lengthArrays);
    nodesY = zeros(1, lengthArrays);
    nodesZ = zeros(1, lengthArrays);
    nodesType = zeros(1, lengthArrays);
    nodesID = zeros(1, lengthArrays);
    nodesDepth = zeros(lengthArrays, 2);
    
    for index = 2:(size(data.textdata,1))
        nodesLabel(index-1) = lower(data.textdata(index, 8));
        nodesX(index-1) = str2double(data.textdata(index, 3));
        nodesY(index-1) = str2double(data.textdata(index, 4));
        nodesZ(index-1) = str2double(data.textdata(index, 5));
        nodesType(index - 1) = str2double(data.textdata(index,2));
        nodesID(index - 1) = str2double(data.textdata(index,1));
        nodesDepth(index - 1, :) = [-1, -1];
        
        parent = data.textdata(strcmp(data.textdata(:, 1), data.textdata(index,6)), :);
        if  ~isempty(parent)
            edgeArray(edgeIndex, :) = [str2double(parent(1)), str2double(data.textdata(index,1))];
            edgeDistFromParent(edgeIndex) = norm([nodesX(index-1), nodesY(index-1), nodesZ(index-1)] - [str2double(parent(3)), str2double(parent(4)), str2double(parent(5))]);
            edgeIndex = edgeIndex + 1;
        else
            rootNodeID = nodesID(index - 1);
        end
        
        if contains(nodesLabel(index-1), 'roi')
            nodesColor(index-1, :) = [1, 0, 0];
        else
             nodesColor(index-1, :) = [0, 0, 1];
        end
        
        nodesLabel(index-1) = replace(nodesLabel(index-1), '_', '');
        nodesLabel(index-1) = replace(nodesLabel(index-1), '-', '');
    end
    
    EdgeTable = table(edgeArray, edgeDistFromParent','VariableNames',{'EndNodes', 'Weight'});
    NodeTable = table(nodesID', nodesLabel', nodesX', nodesY', nodesZ', nodesType', nodesDepth, nodesColor,'VariableNames',{'ID','Name', 'X', 'Y', 'Z', 'Type', 'Depth', 'ColorN'});
    gRoi = graph(EdgeTable, NodeTable);
    
    gRoi = setDepth(gRoi, rootNodeID, 0, 1);
      
%     selectedROINames = setSubGraph1Depth(gRoi, rootNodeID);
    selectedROI = gRoi.Nodes(contains(gRoi.Nodes.Name, 'roi'), :);
    
    if doComboForCloseRoi
        
        [gRoi, indexNodesList] = combRoiClose(gRoi, rootNodeID, rootNodeID, zeros(1, max(gRoi.Nodes.ID)));
        
        selectedROI = gRoi.Nodes(contains(gRoi.Nodes.Name, 'roi') & (~contains(gRoi.Nodes.Name, 'old')), :);
   
    end

    plotGraphWithROI(gRoi, [outputpath, '\GraphWithROI'], gRoi.Nodes.ColorN, {'Graph ROI"s'})
    
    plotTree(gRoi, outputpath);
end

function gRoi = setDepth(gRoi, rootNodeID, depthBefore, indexB)
    gRoi.Nodes(rootNodeID, :).Depth = [depthBefore, indexB];
    nid = neighbors(gRoi,rootNodeID);
    nid = nid(gRoi.Nodes(nid, :).Depth(: ,1) == -1);
    for index = 1:length(nid)
        if (length(nid) == 1)
            gRoi = setDepth(gRoi, nid(index), depthBefore, index);
        else
            gRoi = setDepth(gRoi, nid(index), depthBefore + 1, index);    
        end
    end
end

function [gRoi, indexNodesList] = combRoiClose(gRoi, fNodeID, sNodeID, indexNodesList)
    if contains(gRoi.Nodes.Name(fNodeID), 'roi') && contains(gRoi.Nodes.Name(sNodeID), 'roi')
        if gRoi.Nodes(fNodeID, :).Depth(1,1) == gRoi.Nodes(sNodeID, :).Depth(1,1)
            indexOut = findedge(gRoi,fNodeID, sNodeID);
            if gRoi.Edges.Weight(indexOut) <= 30
                gRoi.Nodes(fNodeID, :).Name = {[gRoi.Nodes.Name{fNodeID} '&' gRoi.Nodes.Name{sNodeID}]};
                gRoi.Nodes(fNodeID, :).X = (gRoi.Nodes(fNodeID, :).X + gRoi.Nodes(sNodeID, :).X) ./ 2;
                gRoi.Nodes(fNodeID, :).Y = (gRoi.Nodes(fNodeID, :).Y + gRoi.Nodes(sNodeID, :).Y) ./ 2;
                gRoi.Nodes(fNodeID, :).Z = (gRoi.Nodes(fNodeID, :).Z + gRoi.Nodes(sNodeID, :).Z) ./ 2;
                
                gRoi.Nodes(sNodeID, :).Name = {['old' gRoi.Nodes.Name{sNodeID}]};
                
                tOut = neighbors(gRoi,sNodeID);
                tOut(tOut == fNodeID) = [];
                
                gRoi = rmedge(gRoi,ones(1, length(tOut)) * sNodeID,tOut);
                gRoi = rmedge(gRoi, fNodeID, sNodeID);
                
                tOutSec = neighbors(gRoi,fNodeID);
                tOutSec(tOutSec == sNodeID) = [];
                gRoi = rmedge(gRoi,tOutSec, ones(1, length(tOutSec)) * fNodeID);
                
                tOut = [tOut, tOutSec];
                wOut = zeros(1, length(tOut));
                sOutNew = ones(1, length(tOut)) * fNodeID;
                
                for i = 1:length(tOut)
                    wOut(i) = norm([gRoi.Nodes.X(tOut(i)), gRoi.Nodes.Y(tOut(i)), gRoi.Nodes.Z(tOut(i))] -...
                        [gRoi.Nodes.X(fNodeID), gRoi.Nodes.Y(fNodeID), gRoi.Nodes.Z(fNodeID)]);
                end
                
                gRoi = addedge(gRoi,sOutNew,tOut,wOut);
%                 gRoi = rmnode(gRoi,sNodeID);
                indexNodesList(sNodeID) = 1;
                sNodeID = fNodeID;
            end
        end
    end
    
    indexNodesList(sNodeID) = 1;
    nid = neighbors(gRoi,sNodeID);
    for index = 1:length(nid)
        if indexNodesList(nid(index)) ~= 1
            [gRoi, indexNodesList] = combRoiClose(gRoi, sNodeID, nid(index), indexNodesList);
        end
    end
end

function selectedROINames = setSubGraph1Depth(gRoi, rootRoiNode)
    fDepth = gRoi.Nodes(gRoi.Nodes.Depth(: ,1) == 0, :);
    sDepth = gRoi.Nodes(gRoi.Nodes.Depth(: ,1) == 1, :);
    tDepth = gRoi.Nodes(gRoi.Nodes.Depth(: ,1) == 2, :);
    wDepth = gRoi.Nodes(gRoi.Nodes.Depth(: ,1) == 3, :);
    
    selectedROI = gRoi.Nodes(contains(gRoi.Nodes.Name, 'roi'), :);
    is_set = zeros(1, size(selectedROI, 1));
    selectedROIPath = cell(1, size(selectedROI, 1));
    
    selectedROINames = cell(1, size(selectedROI, 1));
    index_roi_name = 1;
    
    for roi_i = 1:size(selectedROI,1)
        selectedROIPath{roi_i} = shortestpath(gRoi, rootRoiNode, selectedROI(roi_i, :).ID);
    end
%     for fI = 1:size(fDepth, 1)
        for sI = 1:size(sDepth, 1)
            for tI = 1:size(tDepth, 1)               
                for wI = 1:size(wDepth, 1)
                    for roi_i = find(is_set == 0)
                        if sum(selectedROIPath{roi_i} == tDepth(tI, :).ID) == 1 && ...
                            sum(selectedROIPath{roi_i} == fDepth(end, :).ID) == 1 && ...
                            sum(selectedROIPath{roi_i} == sDepth(sI, :).ID) == 1 && ...
                            sum(selectedROIPath{roi_i} == wDepth(wI, :).ID) == 1
                            selectedROINames(index_roi_name) = selectedROI(roi_i, :).Name;
                            index_roi_name = index_roi_name + 1;
                            is_set(roi_i) = 1;
                        end
                    end
                end
                
                for roi_i = find(is_set == 0)
                    if sum(selectedROIPath{roi_i} == tDepth(tI, :).ID) == 1 && ...
                        sum(selectedROIPath{roi_i} == fDepth(end, :).ID) == 1 && ...
                        sum(selectedROIPath{roi_i} == sDepth(sI, :).ID) == 1
                        selectedROINames(index_roi_name) = selectedROI(roi_i, :).Name;
                        index_roi_name = index_roi_name + 1;
                        is_set(roi_i) = 1;
                    end
                end
            end
            
            for roi_i = find(is_set == 0)
                if sum(selectedROIPath{roi_i} == fDepth(end, :).ID) == 1 && ...
                    sum(selectedROIPath{roi_i} == sDepth(sI, :).ID) == 1
                    selectedROINames(index_roi_name) = selectedROI(roi_i, :).Name;
                    is_set(roi_i) = 1;
                    index_roi_name = index_roi_name + 1;
                end
            end
        end

        for roi_i = find(is_set == 0)
            if sum(selectedROIPath{roi_i} == fDepth(end, :).ID) == 1 
                selectedROINames(index_roi_name) = selectedROI(roi_i, :).Name;
                is_set(roi_i) = 1;
                index_roi_name = index_roi_name + 1;
            end
        end
%     end
end

function plotTree(gRoi, outputpath)
    figTree = figure;
    
    % Create axes
    axes1 = axes('Parent',figTree);
    hold(axes1,'on');

    % Create zlabel
    zlabel({'Z'});

    % Create ylabel
    ylabel({'Y'});

    % Create xlabel
    xlabel({'X'});

    % Create title
    title({'Neuron Tree'});

    view(axes1,[-39.6 18.8]);
    box(axes1,'on');
    grid(axes1,'on');
    
    recursivePlot(gRoi)
    
    [leg, legColor] = getSwcLegPlot();
    legend(leg, legColor);
    %TODO------------------------------
    % Save Tree Plot
    mysave(figTree, [outputpath, '\3DTreeWithROI']);
    %TODO------------------------------
end

function [leg, legColor] = getSwcLegPlot()
    leg = zeros(5, 1);
    leg(1) = plot(0,0, 'color', [0,0,0]);
    leg(2) = plot(0,0,'color', [1,0.6,0.6]);
    leg(3) = plot(0,0,'color', [0.6,0.6,1]);
    leg(4) = plot(0,0,'color', [0.6,1,1]);
    leg(5) = plot(0,0,'color', [0.6,1,0.6]);
    legColor = {'undefined', 'soma', 'axon', '(basal) dendrite', 'apical dendrite'};
end

function [color, name] = getSwcTypeColor(type)
    switch type
        case 0
            color = [0,0,0];
            name = 'undefined';
        case 1
            color = [1,0.6,0.6];
            name = 'soma';
        case 2
            color = [0.6,0.6,1];
            name = 'axon';
        case 3
            color = [0.6,1,1];
            name = '(basal) dendrite';
        case 4
            color = [0.6,1,0.6];
            name = 'apical dendrite';
        case 5
            color = [0,0,0];
            name = 'undefined';
        
    end
end

function recursivePlot(gRoi)   
  for nIndex = 1:size(gRoi.Nodes.Name, 1)
    if contains(gRoi.Nodes.Name{nIndex}, 'roi')
        scatter3(gRoi.Nodes.X(nIndex), gRoi.Nodes.Y(nIndex), gRoi.Nodes.Z(nIndex), '*', 'blue', 'HandleVisibility','off');
        text(gRoi.Nodes.X(nIndex), gRoi.Nodes.Y(nIndex), gRoi.Nodes.Z(nIndex), gRoi.Nodes.Name{nIndex}, 'color', 'blue')
    end
  end
    
  for index = 1:size(gRoi.Edges.EndNodes,1)
    fNode = gRoi.Nodes(findnode(gRoi,gRoi.Edges.EndNodes(index, 1)), :);
    sNode = gRoi.Nodes(findnode(gRoi,gRoi.Edges.EndNodes(index, 2)), :);
    
    [color, ~] = getSwcTypeColor(sNode.Type(1));
    plot3([fNode.X(1), sNode.X(1)], [fNode.Y(1), sNode.Y(1)], [fNode.Z(1), sNode.Z(1)], 'color', color, 'HandleVisibility','off'); 
  end
end