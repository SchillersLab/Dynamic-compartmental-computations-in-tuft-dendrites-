function [roiActivityDistanceMatrixByH, curEventPearson] = calcROIDistanceInActivity_WindowEventPearson_PerEvent(roiActivity, roiActivityNames, selectedROI, all_event_table, typeCom, clusterNum, clusterCount)
    for e_i = 1:size(all_event_table, 1)
        locationToCompareByH = all_event_table.start(e_i):min(all_event_table.event_end(e_i), all_event_table.pks(e_i) + 20);             
        
         curEventPearson(e_i, :, :) = ones(length(selectedROI), length(selectedROI));
         
         for index = 1:length(selectedROI)
            activitylocation = strcmpi(roiActivityNames, selectedROI{index});

            for secIndex = (index+1):length(selectedROI)
                activitySeclocation = strcmp(roiActivityNames, selectedROI{secIndex});
     
                currentROIActivityByH = roiActivity(locationToCompareByH, activitylocation);
                secROIActivityByH = roiActivity(locationToCompareByH, activitySeclocation);

                corrEventsPeaksROIByH = corr([currentROIActivityByH, secROIActivityByH], 'type', 'Pearson');
                
                if isnan(corrEventsPeaksROIByH(1, 2))
                    corrEventsPeaksROIByH(1, 2) = 0;
                end
                
                curEventPearson(e_i, index, secIndex) = corrEventsPeaksROIByH(1, 2);
                curEventPearson(e_i, secIndex, index) = corrEventsPeaksROIByH(1, 2);
            end
         end        
    end

    roiActivityDistanceMatrixByH = squeeze(mean(curEventPearson, 1, 'omitnan'));   
end