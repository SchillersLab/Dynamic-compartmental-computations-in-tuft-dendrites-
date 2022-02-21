function [corGlm, jsGlm] = calcGlmCorrelationMatrix(selectedROI, cont, inds, roiNamesList, timeSeg)
    corGlm(:, :, 1) = ones(length(selectedROI), length(selectedROI));
    corGlm(:, :, 2) = zeros(length(selectedROI), length(selectedROI));
    
    jsGlm(:, :, 1) = zeros(length(selectedROI), length(selectedROI));
    jsGlm(:, :, 2) = zeros(length(selectedROI), length(selectedROI));
    
    for i_roi = 1:length(selectedROI)
        for k_roi = (i_roi + 1):length(selectedROI)
            roiNumI=sscanf(selectedROI{i_roi}, 'roi%05d');
            locationI_glm = find(roiNumI == roiNamesList);
            
            roiNumK=sscanf(selectedROI{k_roi}, 'roi%05d');
            locationK_glm = find(roiNumK == roiNamesList);
            
            contLocationI = find(inds{timeSeg} == locationI_glm);
            contLocationK = find(inds{timeSeg} == locationK_glm);
            
            if isempty(contLocationI) || isempty(contLocationK) 
                corGlm(i_roi, k_roi, 1) = nan;
                corGlm(i_roi, k_roi, 2) = nan;
                corGlm(k_roi, i_roi, 1) = nan;
                corGlm(k_roi, i_roi, 2) = nan;
                
                jsGlm(i_roi, k_roi, 1) = nan;
                jsGlm(i_roi, k_roi, 2) = nan;
                jsGlm(k_roi, i_roi, 1) = nan;
                jsGlm(k_roi, i_roi, 2) = nan;
            else
                [corrR, pR] = corr([cont{timeSeg}(contLocationI, :)', cont{timeSeg}(contLocationK, :)'], 'type', 'Pearson');
                corGlm(i_roi, k_roi, 1) = corrR(1,2);
                corGlm(i_roi, k_roi, 2) = pR(1,2);
                corGlm(k_roi, i_roi, 1) = corrR(1,2);
                corGlm(k_roi, i_roi, 2) = pR(1,2);
                
                M = (cont{timeSeg}(contLocationI, :) + cont{timeSeg}(contLocationK, :)) ./ 2;
                KL1 = nansum(cont{timeSeg}(contLocationI, :) .* log(cont{timeSeg}(contLocationI, :) ./ M));
                KL2 = nansum(cont{timeSeg}(contLocationK, :) .* log(cont{timeSeg}(contLocationK, :) ./ M));
                
                JS = 0.5 * KL1 + 0.5 * KL2;
                jsGlm(i_roi, k_roi, 1) = JS;
                jsGlm(i_roi, k_roi, 2) = 0;
                jsGlm(k_roi, i_roi, 1) = JS;
                jsGlm(k_roi, i_roi, 2) = 0;
            end           
            
        end
    end
end
               