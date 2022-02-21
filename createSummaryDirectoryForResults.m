function createSummaryDirectoryForResults(clusterCount, outputpath, clusterType)
    
    mkdir([outputpath '\SummaryResults\']);
    
    for i_cluster = 0:clusterCount
         if i_cluster == 0
            roiActivityPeakSize = 'All';
        else
            roiActivityPeakSize = ['cluster', num2str(i_cluster)];
         end
        
        copyfile([outputpath '\AnalysisResultsPresentation.pptm'], [outputpath '\SummaryResults\AnalysisResultsPresentation.pptm']);
        
        copyfile([outputpath '\eventsCaSummary.csv'], [outputpath '\SummaryResults\eventsCaSummary.csv']);
        
        filesList_pptm = dir([outputpath '\' roiActivityPeakSize '\' clusterType '\TtestResultsPresentation_*.pptm']);
        filesList_csv = dir([outputpath '\' roiActivityPeakSize '\' clusterType '\DendriticDistVSActivityDistStatistics_*.csv']);
        
        for i = 1:length(filesList_pptm)
            copyfile([filesList_pptm(i).folder, '\', filesList_pptm(i).name], [outputpath '\SummaryResults\' filesList_pptm(i).name]);
        end        
        
        for i = 1:length(filesList_csv)
            copyfile([filesList_csv(i).folder, '\', filesList_csv(i).name], [outputpath '\SummaryResults\' filesList_csv(i).name]);
        end
    end
end