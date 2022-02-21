function createPPForHS
    outputpath = "C:\Users\Jackie\Dropbox (Technion Dropbox)\Yara\Layer 5_Analysis\Shay\SM04\08_18_19_tuft_Final_Version\Analysis\N2\Structural_VS_Functional\2-11-20\Run1\no_behave\HSActivity\HS\";
    outputpath = char(outputpath);
    
    clusterList = {'c0', 'c1', 'c2', 'c3'};
%     clusterList = {'c0', 'c1', 'c2', 'c3', 'c4'};
    

    import mlreportgen.ppt.*
    ppt = Presentation([outputpath '\AnalysisResultsPresentation'], 'AnalysisP.potm');
    open(ppt);
    currentResultsSlideSt = add(ppt, 'Analysis_St');
    currentResultsSlide = add(ppt, 'AnalysisP');
    currentResultsSlide_s = add(ppt, 'AnalysisP');
    
    replace(currentResultsSlideSt.Children(1), Picture([outputpath '\all\GraphWithROI.tif']));       
    replace(currentResultsSlideSt.Children(2), Picture([outputpath '\all\DistMatrixROIStructure.tif']));       
    
    add(currentResultsSlide.Children(end), Paragraph('Cluster By Pks'));
    add(currentResultsSlide_s.Children(end), Paragraph('Cluster By Pks'));
    
    index_presentaion = 1;
    
    for i = 1:length(clusterList)
        replace(currentResultsSlide.Children(end - 1), Picture([outputpath '\' clusterList{i} '\GraphWithROI_2.tif']));
        replace(currentResultsSlide_s.Children(end - 1), Picture([outputpath '\' clusterList{i} '\GraphWithROI_4.tif']));

        replace(currentResultsSlide.Children(index_presentaion), Picture([outputpath '\' clusterList{i} '\DistMatrixActivity_HS.tif']));       
        replace(currentResultsSlide.Children(index_presentaion + 1), Picture([outputpath '\' clusterList{i} '\ActivityDistVSDendriticDistForROI_HS_eventsSize' clusterList{i} '_numofTreeDepth2.tif']));
        replace(currentResultsSlide_s.Children(index_presentaion), Picture([outputpath '\' clusterList{i} '\DistMatrixActivity_HS.tif']));       
        replace(currentResultsSlide_s.Children(index_presentaion + 1), Picture([outputpath '\' clusterList{i} '\ActivityDistVSDendriticDistForROI_HS_eventsSize' clusterList{i} '_numofTreeDepth4.tif']));

        index_presentaion = index_presentaion + 2;
    end
 
    close(ppt);
    
end