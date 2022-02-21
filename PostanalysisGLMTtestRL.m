function PostanalysisGLMTtestRL(tableResults, firstIndexPass)
    timeSeg = 1;
    distanceAngle = [];
    degnificantTypes = [];
    for i = 1:size(tableResults)
        cur = load([tableResults.RunLocation{i}, '\glmResultsSummary.mat']);
        gp = load([tableResults.RunLocation{i}, '\glmResHist.mat'], 'generalProperty');
        uniqC = unique(gp.generalProperty.RoiSplit_I1);
        
        [~,maxL] = sort(mean(cur.cont{timeSeg}(gp.generalProperty.RoiSplit_I1(cur.inds{timeSeg}) == uniqC(1), :)), 'descend');
        [~,maxR] = sort(mean(cur.cont{timeSeg}(gp.generalProperty.RoiSplit_I1(cur.inds{timeSeg}) == uniqC(2), :)), 'descend');
        
        [maxValuesAll, maxAll] = sort(mean(cur.cont{timeSeg}(:, :)), 'descend');
        cumsumTestAll = cumsum(maxValuesAll);
        
        finddiffforLR = setdiff(maxL(1:firstIndexPass), maxR(1:firstIndexPass));
        degnificantTypes(end+1) = ~isempty(finddiffforLR);
        
        corrContributionLeft = mean(cur.cont{timeSeg}(gp.generalProperty.RoiSplit_I1(cur.inds{timeSeg}) == uniqC(1), :));
        corrContributionR = mean(cur.cont{timeSeg}(gp.generalProperty.RoiSplit_I1(cur.inds{timeSeg}) == uniqC(2), :));
        angleListPerComponent = 0:(180/(length(corrContributionLeft)-1)):180;
        angleListPerComponent = angleListPerComponent .* (pi ./ 180);
        currROISelectivityIndexLeft = nan(length(corrContributionLeft), 2);
        currROISelectivityIndexR = nan(length(corrContributionR), 2);

        for i_cont = 1:length(corrContributionLeft)
            currROISelectivityIndexLeft(i_cont, 1) = corrContributionLeft(i_cont)*cos(angleListPerComponent(i_cont));
            currROISelectivityIndexLeft(i_cont, 2) = corrContributionLeft(i_cont)*sin(angleListPerComponent(i_cont));
            currROISelectivityIndexR(i_cont, 1) = corrContributionR(i_cont)*cos(angleListPerComponent(i_cont));
            currROISelectivityIndexR(i_cont, 2) = corrContributionR(i_cont)*sin(angleListPerComponent(i_cont));
        end

        leftSelectivity = mean(currROISelectivityIndexLeft);
        RSelectivity = mean(currROISelectivityIndexR);
        
        angleI = atan2(leftSelectivity(2),leftSelectivity(1)).*(180./pi);
        angleK = atan2(RSelectivity(2),RSelectivity(1)).*(180./pi);

        distanceAngle(end+1) = abs(angleK - angleI);
%         currentIndexPass = [];
%         for j = 1:length(cur.typesU)
%             [p,h] = ranksum(cur.cont{timeSeg}(gp.generalProperty.RoiSplit_I1(cur.inds{timeSeg}) == uniqC(1), j), cur.cont{timeSeg}(gp.generalProperty.RoiSplit_I1(cur.inds{timeSeg}) == uniqC(2), j));
%             
%             if h == 1
%                 degnificantTypes(end+1) = cur.typesU(j);
%                 currentIndexPass(end+1) = j;
%             end
%         end
%         
%         tableResults.SigTypes(i) = {cur.typesU(currentIndexPass)};
    end

end

% p_perpredictor_running = [0.15,0.1531,0.0295,0.028,0.0365,0.0381,0.0277,0.0297,0.0167,0.0232,0.2114,0.2098,0.0462];
% predictorsRunning = {'onset', 'offset', 'speed2', 'speed3', 'accel2', 'accel3', 'speed', 'accel', 'walk', 'rest', 'posaccel', 'negaccel', 'xlocation'};
% 
% totalEquales = 0;
% for i = 1:1000
%     selected1_F = randsample(length(predictorsRunning), 1, true, p_perpredictor_running);
%     selected2_F = randsample(length(predictorsRunning), 1, true, p_perpredictor_running);
%     selected1_S = randsample(length(predictorsRunning), 1, true, p_perpredictor_running);
%     selected2_S = randsample(length(predictorsRunning), 1, true, p_perpredictor_running);
%     if selected1_S == selected2_S & selected1_F == selected2_F
%     totalEquales = totalEquales + 1;
%     end
% end
% totalEquales ./ 1000
