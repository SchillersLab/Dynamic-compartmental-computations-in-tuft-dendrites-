function mainRunnerSaveActivityFileForHS()
    outputpath = "C:\Users\Jackie\Dropbox (Technion Dropbox)\Yara\Layer 5_Analysis\Yara's Data For Shay\SM04\08_18_19_tuft_Final_Version\Analysis\N1\Structural_VS_Functional\cluster4pks10\";
    outputpath = char(outputpath);
   
     re = load([outputpath, '\roiActivityRawData.mat']);
     roiActivity = re.roiActivity;
     roiActivityNames = re.roiActivityNames;
     
     lc = load([outputpath, '\roiActivity_comb.mat'], 'allEventsTable', 'roiActivity_comb');
     
     classes = unique(lc.allEventsTable.clusterByH);
     classes(end + 1) = 0;
     
     size_to_send = 3000;
     tmp = roiActivity(1:size_to_send, :);
     save([outputpath, '\roiActivityRaw_3000F.mat'], 'tmp', 'roiActivityNames')     
     
     tmp = [];
     
     for i = 1:length(classes)
        nameC = ['cluster_', num2str(classes(i))];
        resultsData.(nameC) = ones(size(roiActivity)) .* -100;
        resultsData_comb.(nameC) = ones(size(lc.roiActivity_comb)) .* -100;        
     end
     
     for index = 1:size(lc.allEventsTable, 1)
         nameC = ['cluster_', num2str(lc.allEventsTable.clusterByH(index))];
        
         resultsData.(nameC)(lc.allEventsTable.start(index):(min(lc.allEventsTable.event_end(index), lc.allEventsTable.pks(index) + 20)), :) = roiActivity(lc.allEventsTable.start(index):(min(lc.allEventsTable.event_end(index), lc.allEventsTable.pks(index) + 20)), :);
         resultsData.('cluster_0')(lc.allEventsTable.start(index):(min(lc.allEventsTable.event_end(index), lc.allEventsTable.pks(index) + 20)), :) = roiActivity(lc.allEventsTable.start(index):(min(lc.allEventsTable.event_end(index), lc.allEventsTable.pks(index) + 20)), :);
        
         resultsData_comb.(nameC)(lc.allEventsTable.start(index):(min(lc.allEventsTable.event_end(index), lc.allEventsTable.pks(index) + 20)), :) = lc.roiActivity_comb(lc.allEventsTable.start(index):(min(lc.allEventsTable.event_end(index), lc.allEventsTable.pks(index) + 20)), :);
         resultsData_comb.('cluster_0')(lc.allEventsTable.start(index):(min(lc.allEventsTable.event_end(index), lc.allEventsTable.pks(index) + 20)), :) = lc.roiActivity_comb(lc.allEventsTable.start(index):(min(lc.allEventsTable.event_end(index), lc.allEventsTable.pks(index) + 20)), :);
     
     end
     
     
     
     for i = 1:length(classes)
        nameC = ['cluster_', num2str(classes(i))];
        tmp = resultsData.(nameC);
        tmp(tmp(:, 1) == -100, :) = [];
        
        save([outputpath, '\roiActivityRaw_ByEvents_' num2str(classes(i)) '.mat'], 'tmp', 'roiActivityNames')
        
        tmp = [];
        tmp = resultsData_comb.(nameC);
        tmp(tmp(:, 1) == -100, :) = [];
        
        save([outputpath, '\roiActivityComb_ByEvents_' num2str(classes(i)) '.mat'], 'tmp', 'roiActivityNames')     
     end
         
end