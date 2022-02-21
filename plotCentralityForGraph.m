function cent_weighted = plotCentralityForGraph(gRoi, roiActivityDistanceMatrix, selectedROI, outputpath, isFixNeeded, selectedROISplitDepth1)
    if isFixNeeded
        roiActivityDistanceMatrix = roiActivityDistanceMatrix + 1;
%         roiActivityDistanceMatrix = abs(roiActivityDistanceMatrix);
        roiActivityDistanceMatrix(isnan(roiActivityDistanceMatrix)) = 0;
    end
    
    [cent_weighted, cent_notweighted] = graph_centrality(roiActivityDistanceMatrix, selectedROISplitDepth1);
    
    centrality_measures = fieldnames(cent_notweighted);


    for name_i = 1:length(centrality_measures)
%         PlotGraphWithCentrality(cent_notweighted.(centrality_measures{name_i}), gRoi, selectedROI,...
%             fullfile(outputpath, centrality_measures{name_i}, 'noW'), centrality_measures{name_i}, 'noW');
%         
%         
%         PlotGraphWithCentrality(cent_weighted.(centrality_measures{name_i}), gRoi, selectedROI,...
%             fullfile(outputpath, centrality_measures{name_i}, 'W'), centrality_measures{name_i}, 'W');
%         
        t = table(selectedROI,selectedROISplitDepth1, cent_weighted.(centrality_measures{name_i}),...
            cent_notweighted.(centrality_measures{name_i}), 'VariableNames', {'ROI', 'Labels', 'CentralityW', 'CentralityNW'});
        
        writetable(t, fullfile(outputpath, [centrality_measures{name_i}, '.csv']));
    end     
    
end

function PlotGraphWithCentrality(cent, gRoi, selectedROI, outputpath, measure, type)
     nodesColor = zeros(length(gRoi.Nodes.Name),3);
    
    rgbColors = vals2colormap(cent, 'jet', [min(cent), max(cent)]);
    
    for index = 1:length(selectedROI)
        locRoi = find(strcmp(gRoi.Nodes.Name, selectedROI(index)));
        nodesColor(locRoi, :) = rgbColors(index,:);
    end
    
    labelsNames = gRoi.Nodes.Name;
    f_r = find(contains(gRoi.Nodes.Name', 'roi') == 1);
    
    for index_roi = f_r
        labelsNames(index_roi) = {sprintf('roi%d', sscanf(gRoi.Nodes.Name{index_roi}, 'roi%d'))};
    
    end
  
    
    figGraph = figure;
    plot(gRoi, 'NodeLabel', labelsNames, 'NodeFontWeight', 'bold', 'NodeColor', nodesColor);
    title({'Roi centrality', measure, type});
    colorbar;
    colormap('jet');
    
    if min(cent) == max(cent)
        caxis([min(cent), max(cent)+1]);
    else
        caxis([min(cent), max(cent)]);
    end
    
    figGraph.Position = [figGraph.Position(1), figGraph.Position(2), figGraph.Position(3) + 200, figGraph.Position(4) + 70]; 
    
    fileName2 = [outputpath, '\GraphWithROI_ColoredCentrality'];
    mysave(figGraph, fileName2);
   
end