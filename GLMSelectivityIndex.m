function [glmMatrixVector,glmMatrixAngle] = GLMSelectivityIndex(selectedROI, cont, inds, roiNamesList, timeSeg)
    roiSelectivityIndex = nan(length(selectedROI), 2);
    
    contSum = nanmean(cont{timeSeg}, 1);
    [~, max3Index] = sort(contSum,'descend', 'MissingPlacement','last');
    max3Index = max3Index(1:3);   
    
    for i_roi = 1:length(selectedROI)        
        roiNumI=sscanf(selectedROI{i_roi}, 'roi%05d');
        locationI_glm = find(roiNumI == roiNamesList);            
        contLocationI = find(inds{timeSeg} == locationI_glm);
            
        if ~isempty(contLocationI)
            corrContribution = cont{timeSeg}(contLocationI, max3Index);

            angleListPerComponent = 0:(180/(length(corrContribution)-1)):180;
            angleListPerComponent = angleListPerComponent .* (pi ./ 180);
            currROISelectivityIndex = nan(length(corrContribution), 2);
            
            for i_cont = 1:length(corrContribution)
                currROISelectivityIndex(i_cont, 1) = corrContribution(i_cont)*cos(angleListPerComponent(i_cont));
                currROISelectivityIndex(i_cont, 2) = corrContribution(i_cont)*sin(angleListPerComponent(i_cont));
            end
               
            roiSelectivityIndex(i_roi, :) = mean(currROISelectivityIndex);
        end                       
    end
    
    glmMatrixAngle = nan(length(selectedROI), length(selectedROI));
    glmMatrixVector = nan(length(selectedROI), length(selectedROI));
    
    for i_roi = 1:length(selectedROI)
        glmMatrixAngle(i_roi, i_roi) = 0;
        glmMatrixVector(i_roi, i_roi) = 0;
        for k_roi = (i_roi + 1):length(selectedROI)
            distanceVectors = norm(roiSelectivityIndex(i_roi, :) - roiSelectivityIndex(k_roi, :));
            
            angleI = atan2(roiSelectivityIndex(i_roi, 2),roiSelectivityIndex(i_roi, 1)).*(180./pi);
            angleK = atan2(roiSelectivityIndex(k_roi, 2),roiSelectivityIndex(k_roi, 1)).*(180./pi);
                     
            distanceAngle = abs(angleK - angleI);
            
            glmMatrixAngle(i_roi, k_roi) = distanceAngle;
            glmMatrixAngle(k_roi, i_roi) = distanceAngle;
            
            glmMatrixVector(i_roi, k_roi) = distanceVectors;
            glmMatrixVector(k_roi, i_roi) = distanceVectors;
        end
    end
end
               