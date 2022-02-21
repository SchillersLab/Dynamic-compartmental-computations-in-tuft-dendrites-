function [pictureNames, resultsForGLM] = plotROIDistMatrixTreeVSActivity(event_count, gRoi, outputpath, firstBranchROI,mainTreeBranchROI, roiTreeDistanceMatrix, roiActivityDistanceMatrix, do_corrtest, roiActivityDistanceFunction, roiActivityPeakSize, selectedRoi, index_apical, distType, clusterType)
    classesM = unique(mainTreeBranchROI);
    classesF = unique(firstBranchROI);
    
    if all(firstBranchROI == mainTreeBranchROI)
        isMainSubTree = true;
    else
        isMainSubTree = false;
    end
    
    nodesColor = zeros(length(gRoi.Nodes.Name),3);
    tempC_2 = zeros(length(classesF(classesF ~= -1)), 3);
    
    ma_c = getTreeColor('main', -1 , true);
  
    tempC_2(:, 1) = ma_c(:, 1);
    tempC_2(:, 2) = ma_c(:, 2);
    tempC_2(:, 3) = ma_c(:, 3);
     
    for clr = 1:length(classesM)
        curTree = find(mainTreeBranchROI == classesM(clr) & classesM(clr) ~= -1);
        
        for i = 1:length(curTree)
            if mainTreeBranchROI(curTree(i)) == classesM(clr) && classesM(clr) ~= -1
                locRoi = find(strcmp(gRoi.Nodes.Name, selectedRoi(curTree(i))));
                nodesColor(locRoi, :) = getTreeColor('within', (clr), isMainSubTree, length(classesF(classesF ~= -1)));
%                 coutI = coutI + 1;
            end
        end
        
        if (classesM(clr) ~= -1) && ~isMainSubTree
            fID = find(classesF == firstBranchROI(curTree(1)));
        
            nodesColor(classesM(clr), :) = getTreeColor('within', fID, true);
        end
    end
    
    
    nodesColor(classesF(classesF ~= -1), :) = tempC_2;
     
    titlePG = {'Number of subtree ', num2str(length(classesM))};
    fileName2 = [outputpath, '\GraphWithROI_' num2str(length(classesM))];
    plotGraphWithROI(gRoi, fileName2, nodesColor, titlePG)
 
    
    %     Plot ROI Activity VS Tree Distance
    fig = figure;
    hold on;
    subplot(9,1,2:8);
    sgtitle({'ROI Activity VS Tree Distance'});
    
    % Create ylabel
%     xlabel({'Calcium Event ' distType, ['(', num2str(event_count),')']});
    ylabel({'Calcium Event ' distType, ['(', num2str(event_count),')']});

    % Create xlabel
%     ylabel({'Dendritic distance'});
    xlabel({'Dendritic distance'});

    
%     leg = zeros(3, 1);
    leg = [];
%     legColor = []; 
    
    
    index_classes = 1;   
    classesColorName = {};
    resultsForGLM = {};
    
    for index = 1: size(roiTreeDistanceMatrix, 2)
        for secIndex = (index + 1): size(roiTreeDistanceMatrix, 2)
            
            fMainName = 'ND';
            secMainName = 'ND';
            
            fFirstName = 'ND';
            secFirstName = 'ND';
            
            if mainTreeBranchROI(index) ~= -1
                fMainName = gRoi.Nodes(mainTreeBranchROI(index),:).Name{1};
            end
            
            if mainTreeBranchROI(secIndex) ~= -1
                secMainName = gRoi.Nodes(mainTreeBranchROI(secIndex),:).Name{1};
            end
            
            if firstBranchROI(index) ~= -1
                fFirstName = gRoi.Nodes(firstBranchROI(index),:).Name{1};
            end
            
            if firstBranchROI(secIndex) ~= -1
                secFirstName = gRoi.Nodes(firstBranchROI(secIndex),:).Name{1};
            end
            
            if strcmp(fFirstName, 'ND') || strcmp(secFirstName, 'ND')
                color = getTreeColor('ND', -1, true);
                colorName = {'NotInDepth'};
            elseif ~strcmp(fFirstName, secFirstName)
                colorName = {'BetweenMainDepth'};
                color = getTreeColor('main',-1, true);
            elseif strcmp(fMainName, 'ND') || strcmp(secMainName, 'ND')
                color = getTreeColor('ND', -1, true);
                colorName = {'NotInDepth'};
            elseif ~strcmp(fMainName, secMainName)
                colorName = {['BetweenSecondDepth_', fFirstName]};
                color = getTreeColor('within', find(classesF == firstBranchROI(index)), true);
            else
                color = getTreeColor('within', find(classesM == mainTreeBranchROI(index)), isMainSubTree, length(classesF));
                colorName = {fMainName};
            end
               
            hold on;
            
            if contains(colorName, 'BetweenMainDepth') && ~isMainSubTree
                continue;
            end              
            
            if sum((index_apical == index) | (index_apical == secIndex)) > 0
%                 scat = scatter(roiActivityDistanceMatrix(index, secIndex, 1), roiTreeDistanceMatrix(index, secIndex), '*', 'MarkerEdgeColor', color);
                scat = scatter(roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex, 1), '*', 'MarkerEdgeColor', color);
            else
%                 scat = scatter(roiActivityDistanceMatrix(index, secIndex, 1), roiTreeDistanceMatrix(index, secIndex), 20, 'filled', 'MarkerFaceColor', color);
                scat = scatter(roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex, 1), 20, 'filled', 'MarkerFaceColor', color);           
            end
            
            scat.set('UserData', [selectedRoi{index} 'x' selectedRoi{secIndex}])

            colorName = colorName{1};
            colorName = replace(colorName,{'&','-', '_'}, 'x'); 
                
            if isempty(classesColorName) || ...
                    (sum(strcmp(classesColorName, colorName)) == 0)
                classesColorName(index_classes) = {colorName};
                classesColor(index_classes) = {color};
                
                leg(index_classes) = plot(0,0, 'color', color, 'LineWidth', 2.5);
                legColor(index_classes) = {colorName};
                
                resultsT.(colorName) = [];
                index_classes = index_classes + 1;
            end
            
            resultsForGLM(end + 1, :) = {([selectedRoi{index} 'x' selectedRoi{secIndex}]),...
                roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex, 1), roiActivityDistanceMatrix(index, secIndex, 2), mainTreeBranchROI(index), ...
                mainTreeBranchROI(secIndex), (mainTreeBranchROI(secIndex) == mainTreeBranchROI(index) && mainTreeBranchROI(index) ~= -1), roiActivityPeakSize};  
            
            resultsT.(colorName)(end + 1, :) = [roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex, 1), roiActivityDistanceMatrix(index, secIndex, 2)];
        end
    end
 
    if ~exist('resultsT', 'var')
        pictureNames = {'', ''};
        resultsForGLM = cell(1, 5);
        fclose('all');
        close all;
        return;
    end
    
    f_name = fieldnames(resultsT);
    indexT = 1;
    txtSubT = '';
    
    betweenSecLocation  = (contains(f_name, 'BetweenSecondDepth'));
    resultsSplit.Between = [];
    resultsSplit.within = [];
    
    for t = 1:length(f_name)
        if isempty(resultsT.(f_name{t}))
            continue;
        end
        
        locationColor = strcmp(classesColorName, f_name{t});
%         mdl = fitlm(resultsT.(f_name{t})(:,2), resultsT.(f_name{t})(:,1));
        

        if size(resultsT.(f_name{t}), 1) > 1
            mdl = fitglm(resultsT.(f_name{t})(:,1), resultsT.(f_name{t})(:,2));

            if isMainSubTree
                yfit = predict(mdl, resultsT.(f_name{t})(:,1));
        %         plot(resultsT.(f_name{t})(:,2), yfit, 'Color', classesColor{locationColor}, 'LineStyle', '--');
                plot(resultsT.(f_name{t})(:,1), yfit, 'Color', classesColor{locationColor}, 'LineStyle', '--');

                txtSubT = strcat(txtSubT, sprintf(' {\\color[rgb]{%f,%f,%f}R^2=%0.4f}\n\r',classesColor{locationColor}(1),...
                    classesColor{locationColor}(2), classesColor{locationColor}(3), mdl.Rsquared.Ordinary));
            end

            r2(indexT, 1) = {mdl.Rsquared.Ordinary};
            r2(indexT, 2) = f_name(t);

            indexT = indexT + 1;
        end
        
        writetable(array2table(resultsT.(f_name{t})),fullfile(outputpath, ['AVSDForROI_' num2str(t) '_Depth', num2str(length(classesM)), '_', f_name{t}, '.csv']));
        
        if ~isMainSubTree
            if betweenSecLocation(t) == 1
                resultsSplit.Between((end+1):(end + size(resultsT.(f_name{t}), 1)), :) = resultsT.(f_name{t});
            else
                resultsSplit.within((end+1):(end + size(resultsT.(f_name{t}), 1)), :) = resultsT.(f_name{t});
            end
        end
    end
    
    
    if ~isMainSubTree
        
        if size(resultsSplit.Between, 1) > 1
            mdl = fitglm(resultsSplit.Between(:,1), resultsSplit.Between(:,2));

            yfit = predict(mdl, resultsSplit.Between(:,1));
            plot(resultsSplit.Between(:,1), yfit, 'Color', [119, 123, 126] ./ 255, 'LineStyle', '--');

            txtSubT = strcat(txtSubT, sprintf(' {\\color[rgb]{%f,%f,%f}R^2=%0.4f}\n\r',119 ./ 255,...
                123 ./ 255, 126 ./ 255, mdl.Rsquared.Ordinary));

            r2(indexT, 1) = {mdl.Rsquared.Ordinary};
            r2(indexT, 2) = {'Between Sec Depth All'};

            indexT = indexT + 1;
        end
        
        if size(resultsSplit.within, 1) > 1
            md2 = fitglm(resultsSplit.within(:,1), resultsSplit.within(:,2));

            yfit = predict(md2, resultsSplit.within(:,1));
            plot(resultsSplit.within(:,1), yfit, 'Color', [47,79,79] ./ 255, 'LineStyle', '--');

            txtSubT = strcat(txtSubT, sprintf(' {\\color[rgb]{%f,%f,%f}R^2=%0.4f}\n\r',...
                47./ 255, 79./ 255, 79 ./ 255, md2.Rsquared.Ordinary));

            r2(indexT, 1) = {md2.Rsquared.Ordinary};
            r2(indexT, 2) = {'WithIn Sec Depth All'}; 
            indexT = indexT + 1;
        end
    end
    
    if size(resultsForGLM, 1) > 1
        mdAll = fitglm(cell2mat(resultsForGLM(:,2)), cell2mat(resultsForGLM(:,3)));
        yfitAll = predict(mdAll, cell2mat(resultsForGLM(:,2)));
        plot(cell2mat(resultsForGLM(:,2)), yfitAll, 'Color', [0, 0, 0] ./ 255, 'LineStyle', '--');

        txtSubT = strcat(txtSubT, sprintf(' {\\color[rgb]{%f,%f,%f}R^2=%0.4f}\n\r',0,...
            0, 0, mdAll.Rsquared.Ordinary));

        r2(indexT, 1) = {mdAll.Rsquared.Ordinary};
        r2(indexT, 2) = {'All'}; 
        
        title(txtSubT);
        writetable(cell2table(r2),fullfile(outputpath, ['Rsquared_' num2str(t) '_Depth', num2str(length(classesM)), '_', f_name{t}, '.csv']));
    end
      
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn})
    
    legend(leg, legColor);
    legend('Location','bestoutside')
    fig.Position = [fig.Position(1)-250, fig.Position(2)-100, fig.Position(3) + 250, fig.Position(4) + 100];
    roiwithstar = '*';
    
    for in = 1:length(index_apical)
       roiwithstar = [roiwithstar ' ' selectedRoi{index_apical(in)}]; 
    end
    
    annotation('textbox', [0, 0.2, 0, 0], 'string', roiwithstar, 'FitBoxToText','on')
%     ylim([0,1]);
       
    fileName1 = [outputpath, '\ActivityDistVSDendriticDistForROI_' roiActivityDistanceFunction ,'_eventsSize', roiActivityPeakSize, '_numofTreeDepth', num2str(length(classesM))];
    mysave(fig, fileName1);
    
    plotAndCalcStatisticTest(classesColorName, classesColor, resultsT, leg, legColor, outputpath, classesM, clusterType, roiActivityPeakSize, distType);
 
% ------------------------------------------------------------------------------
    
    %     coutI = 1;
    pictureNames = {fileName1, fileName2};
    
    snapnow;
    fclose('all');
    close all;
    clear resultsT;
end