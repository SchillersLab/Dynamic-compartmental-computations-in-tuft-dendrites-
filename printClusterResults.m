function printClusterResults(clusterVector, clusterCount, meanCombActivity, all_locationFull_pks, all_locationFull_start, all_locationFull_end, all_locationFull_H, outputpath, clusterType)
    f = figure;
    hold on;
    
    sb1 = subplot(8, 1, 1:6);
    hold on;
    title('Mean Activity only events');
    
    plot(meanCombActivity)
    
    plot(all_locationFull_pks, meanCombActivity(all_locationFull_pks), '*r');
    plot(all_locationFull_start, meanCombActivity(all_locationFull_start), '*b');
    plot(all_locationFull_end, meanCombActivity(all_locationFull_end), '*g');
 
    legend('Activity', 'Peaks', 'StartEvent', 'EndEvents');
    
    for i = 1:clusterCount
        plot(all_locationFull_pks(clusterVector == i), all_locationFull_H(clusterVector == i), 'o');
    end
    
    
    xlim([1, size(meanCombActivity, 1)]);
%     ylim([-1, 5]);
    
    sb2 = subplot(8, 1, 8:8);
    imagesc(meanCombActivity');
    colormap(jet);
    caxis(sb1.YLim);
       
    linkaxes([sb1, sb2], 'x');
    
    mysave(f, [outputpath, '\activity_averagePksHistByMLSpike_cluster_', clusterType]);
end