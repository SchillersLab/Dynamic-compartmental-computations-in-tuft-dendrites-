function figGraph = plotGraphWithROI(gRoi, filename, ColorN, titlePG, markerSize)
    if nargin == 4
        markerSize = 3;
    end

    figGraph = figure;
    
    labelsNames = gRoi.Nodes.Name;
    f_r = find(contains(gRoi.Nodes.Name', 'roi') == 1);
    
    for index_roi = f_r
        labelsNames(index_roi) = {sprintf('roi%d', sscanf(gRoi.Nodes.Name{index_roi}, 'roi%d'))};
    
    end
  
    plot(gRoi, 'NodeColor', ColorN, 'NodeLabel', labelsNames, 'NodeFontWeight', 'bold', 'MarkerSize', markerSize);
    title(titlePG);
    figGraph.Position = [figGraph.Position(1), figGraph.Position(2), figGraph.Position(3) + 200, figGraph.Position(4) + 70]; 
    mysave(figGraph, filename);
    
end