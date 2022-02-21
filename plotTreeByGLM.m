function plotTreeByGLM(gRoi, outputpathCurr, selectedROI, cont, inds, roiNamesL, typesU, timeS, selectedROISplitDepth1, colorMatrix1)
    plotColorByMainCont(gRoi, outputpathCurr, selectedROI, cont, inds, roiNamesL, typesU, timeS);
    
    summaryPerRoi = table(cell(0, 1), zeros(0,1), cell(0,1), zeros(0, 1), zeros(0,1));
    summaryPerRoi.Properties.VariableNames = {'Roi_name', 'Main_Direction_value', 'Main_Direction_Name', 'Types_Pass_15_precent', 'SubTreeSplit'};
    for i_roi = 1:length(selectedROI)
        roiNumI=sscanf(selectedROI{i_roi}, 'roi%05d');
        locationI_glm = find(roiNumI == roiNamesL);
        
        contLocationI = find(inds{timeS} == locationI_glm);

        if ~isempty(contLocationI)
            [preferdDir_value, preferdDir_index] = max(cont{timeS}(contLocationI, :));
            summaryPerRoi.Roi_name(end + 1) = selectedROI(i_roi);
            summaryPerRoi.Main_Direction_value(end) = preferdDir_value;
            summaryPerRoi.Main_Direction_Name(end) = typesU(preferdDir_index);
            summaryPerRoi.Types_Pass_15_precent(end) = sum(cont{timeS}(contLocationI, :) > 0.15);
            summaryPerRoi.SubTreeSplit(end) = selectedROISplitDepth1(i_roi);
        end
    end
    
    writetable(summaryPerRoi, [outputpathCurr, '\TableSummarySelectivity.csv']);
    classesU = unique(summaryPerRoi.SubTreeSplit);
    avgC = zeros(1, length(classesU));
    stdC = zeros(1, length(classesU));
   
    colorMatR = zeros(length(classesU), 3);
    
    labels = {};
    for i = 1:length(classesU)
        avgC(i) = mean(summaryPerRoi.Main_Direction_value(summaryPerRoi.SubTreeSplit == classesU(i)), 'omitnan');
        stdC(i) =  std(summaryPerRoi.Main_Direction_value(summaryPerRoi.SubTreeSplit == classesU(i)), 'omitnan');
        colorMatR(i, :) = colorMatrix1(find(summaryPerRoi.SubTreeSplit == classesU(i), 1, 'first'),:);
        
        if classesU(i) == -1
            labels(i) = gRoi.Nodes.Name(1);
        else
            labels(i) = gRoi.Nodes.Name(classesU(i));
        end
    end
    
    f = figure;
    hold on;
    b= bar(1:length(classesU), avgC);    
    b.FaceColor = 'flat';
    b.CData = colorMatR;
   
    errorbar(1:length(classesU), avgC, stdC, 'Color', 'k');
    xticks(1:length(classesU));
    xticklabels(labels);
    xtickangle(90);
    ylabel('R2');
    title({'Mean Primary Glm Value within tree'});
    
    mysave(f, [outputpathCurr, '\SelectivityValuesBySubTree']);
    
    rgbColorsByValues = vals2colormap(summaryPerRoi.Main_Direction_value, 'jet', [min(summaryPerRoi.Main_Direction_value), max(summaryPerRoi.Main_Direction_value)]);
    rgbColorsByCount = vals2colormap(summaryPerRoi.Types_Pass_15_precent, 'jet', [min(summaryPerRoi.Types_Pass_15_precent), max(summaryPerRoi.Types_Pass_15_precent)]);
    
    nodesColorByV = zeros(length(gRoi.Nodes.Name),3);
    nodesColorByC = zeros(length(gRoi.Nodes.Name),3);
    
    markerSize = ones(size(gRoi.Nodes, 1), 1) * 1;
    for index = 1:size(summaryPerRoi, 1)
        locRoi = find(strcmp(gRoi.Nodes.Name, summaryPerRoi.Roi_name(index)));
        nodesColorByV(locRoi, :) = rgbColorsByValues(index,:);
        nodesColorByC(locRoi, :) = rgbColorsByCount(index,:);
        markerSize(locRoi) = 6;
    end
    
    figGraph = plotGraphWithROI(gRoi, [outputpathCurr, '\GraphWithROIByGlmMainDirValue'], nodesColorByV, {'Graph ROI"s Primary Glm Value'}, markerSize);
    colorbar;
    colormap('jet');
    
    if min(summaryPerRoi.Main_Direction_value) == max(summaryPerRoi.Main_Direction_value)
        caxis([min(summaryPerRoi.Main_Direction_value), max(summaryPerRoi.Main_Direction_value)+1]);
    else
        caxis([min(summaryPerRoi.Main_Direction_value), max(summaryPerRoi.Main_Direction_value)]);
    end
    mysave(figGraph, [outputpathCurr, '\GraphWithROIByGlmMainDirValueBar']);
  
    
    figGraph = plotGraphWithROI(gRoi, [outputpathCurr, '\GraphWithROIByGlmPass15PercentCount'], nodesColorByC, {'Graph ROI"s, Summary above 15% contribution'}, markerSize);  
    colorbar;
    colormap('jet');
    
    if min(summaryPerRoi.Types_Pass_15_precent) == max(summaryPerRoi.Types_Pass_15_precent)
        caxis([min(summaryPerRoi.Types_Pass_15_precent), max(summaryPerRoi.Types_Pass_15_precent)+1]);
    else
        caxis([min(summaryPerRoi.Types_Pass_15_precent), max(summaryPerRoi.Types_Pass_15_precent)]);
    end
    
    mysave(figGraph, [outputpathCurr, '\GraphWithROIByGlmMainDirPrecantBar']);
  
end


function plotColorByMainCont(gRoi, outputpathCurr, selectedROI, cont, inds, roiNamesL, typesU, timeS)
    contSum = nanmean(cont{timeS}, 1);
    [~, max3Index] = sort(contSum,'descend', 'MissingPlacement','last');
    max3Index = max3Index(1:3);
    
    maxTypesU = typesU(max3Index);
    
    roiColors = ones(size(gRoi.Nodes, 1), 3);
    markerSize = ones(size(gRoi.Nodes, 1), 1) * 1;
    for i_roi = 1:length(selectedROI)        
        roiNumI=sscanf(selectedROI{i_roi}, 'roi%05d');
        locationI_glm = find(roiNumI == roiNamesL);

        contLocationI = find(inds{timeS} == locationI_glm);

        if ~isempty(contLocationI) 
            contMax = cont{timeS}(contLocationI, max3Index);
            [~, sortI] = sort(contMax,'descend', 'MissingPlacement','last');  

            roiLocation = find(strcmp(gRoi.Nodes.Name, selectedROI{i_roi}));
            rgbC = getColorRGB(num2str(sortI));

            roiColors(roiLocation, :) = rgbC;

            markerSize(roiLocation) = 6;
        end       
    end
    
    plotGraphWithROI(gRoi, [outputpathCurr, '\GraphWithROIByGlm'], roiColors, {'Graph ROI"s'}, markerSize);
    
    hold on;
    p = perms([1,2,3]);
     leg = [];
    for i = 1:size(p, 1)       
        leg(i) = plot(0,0, 'Color', getColorRGB(num2str(p(i, :))), 'LineWidth', 1.5);
        legColor(i) = {[maxTypesU{p(i, 1)}, ', ', maxTypesU{p(i, 2)}, ', ', maxTypesU{p(i, 3)}]};
    end
    
    legend(leg, legColor);
    legend('Location', 'bestoutside');
    
    mysave(gca, [outputpathCurr, '\GraphWithROIByGlm_Leg']);
end
    

function rgbC = getColorRGB(index)
    switch index
        case '1  2  3'
            rgbC = [230, 25, 75] ./ 255;
        case '2  1  3'
            rgbC = [60, 180, 75] ./ 255;
        case '3  1  2'
            rgbC = [0, 130, 200] ./ 255;
        case '1  3  2'
            rgbC = [145, 30, 180] ./ 255;
        case '3  2  1'
            rgbC = [245, 130, 48] ./ 255;
        case '2  3  1'
            rgbC = [128, 128, 128] ./ 255;
    end
end