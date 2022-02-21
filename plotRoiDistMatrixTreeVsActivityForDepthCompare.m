function plotRoiDistMatrixTreeVsActivityForDepthCompare(gRoi, outputpath, mainTreeBranchROI, roiTreeDistanceMatrix, roiActivityDistanceMatrix, roiActivityDistanceFunction, roiActivityPeakSize)
    %     Plot ROI Activity VS Tree Distance
    fig = figure;
    hold on;
    % Create ylabel
    ylabel({'Calcium Event Distance'});

    % Create xlabel
    xlabel({'Dendritic distance'});

    title({'ROI Activity VS Tree Distance'});
%     leg = zeros(3, 1);
    leg = [];
%     legColor = []; 
    classesM = unique(mainTreeBranchROI.Depth(:, 1));
    
    
    for index = 1:length(classesM)
        [classesColor(index, :), legColor{index}] = getDepthColor(classesM(index));
        leg(end+1) = plot(0,0, 'color', classesColor(index, :), 'LineWidth', 2.5);
    end
    
    classesColor(end + 1, :) = [1,0.2,0.2];
    legColor(end + 1) = {'Between Depths'};
    leg(end+1) = plot(0,0, 'color', classesColor(end, :), 'LineWidth', 2.5);
      
    corrIndexMatrix = 1;   
    corrIndexMatrixInsideMainBranch = 1;
    for index = 1: size(roiTreeDistanceMatrix, 2)
        for secIndex = (index + 1): size(roiTreeDistanceMatrix, 2)
            if (mainTreeBranchROI.Depth(index, 1) == mainTreeBranchROI.Depth(secIndex, 1))
                cur_depth = find(classesM == mainTreeBranchROI.Depth(index, 1));
                color = classesColor(cur_depth, :);
                corrMatrixForROIInsideTheMainBranch(corrIndexMatrixInsideMainBranch, :) = [roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex)];
                corrIndexMatrixInsideMainBranch = corrIndexMatrixInsideMainBranch + 1;                        
            else
                color = classesColor(end, :);              
            end
            
            hold on;
            corrMatrixForROI(corrIndexMatrix, :) = [roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex)];
            corrIndexMatrix = corrIndexMatrix + 1;
            scatter(roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex), 'filled', 'MarkerFaceColor', color);
        end
    end

    legend(leg, legColor);
    legend('Location','northwestoutside')
    
%     ylim([0,1]);
       
    mysave(fig, [outputpath, '\ActivityDistVSDendriticDistForROIByDepth_' roiActivityDistanceFunction ,'_eventsSize', roiActivityPeakSize, '_numofTreeDepth', num2str(length(classesM))]);
    
%     coutI = 1;
    nodesColor = zeros(length(gRoi.Nodes.Name),3);
    
    for i = 1:size(mainTreeBranchROI.Depth, 1)
        depth_calc = mainTreeBranchROI.Depth(i, 1);
        locRoi = find(strcmp(gRoi.Nodes.Name, mainTreeBranchROI.Name(i)));
        nodesColor(locRoi, :) = classesColor(depth_calc == classesM, :);    
    end

    plotGraphWithROI(gRoi, [outputpath, '\GraphWithROIDepth_' num2str(length(classesM))], nodesColor, {'Graph Depth'})
    
    snapnow;
    fclose('all');
    close all;
    clear resultsT;
end

function [color, depthName] = getDepthColor(depthNum)
    switch depthNum
        case 1
            color = [0.2,1,0.6];
            depthName = 'Depth 1';
        case 2
            color = [0.2,1,1];
            depthName = 'Depth 2';
        case 3
            color = [0.2,0.6,1];
            depthName = 'Depth 3';
        case 4
            color = [0.2,0.2,1];
            depthName = 'Depth 4';
        case 5
            color = [0.6,0.2,1];
            depthName = 'Depth 5';
        case 6
            color = [1,0.2,1];
            depthName = 'Depth 6';           
        case 7
            color = [1,0.2,0.6];
            depthName = 'Depth 7';
        case 8
            color = [1,0.6,0.2];
            depthName = 'Depth 8';
        otherwise
            color = [0,0,0];
            depthName = 'Depth above 8';
    end
end