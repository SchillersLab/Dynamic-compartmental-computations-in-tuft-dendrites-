function [SpikeTrainClusterSecByData] = getClusterForActivity(data, clusterCount)
    [SpikeTrainClusterByData] = kmeans(data, clusterCount, 'Replicates',5, 'MaxIter', 500);

    clusterMaxValue = [];
    for i = 1:clusterCount
        clusterMaxValue(i) = max(data(SpikeTrainClusterByData == i));
    end

    [~, cluster_sort_index] = sort(clusterMaxValue);
    SpikeTrainClusterSecByData = zeros(1, length(SpikeTrainClusterByData));
    for i = 1:clusterCount
        SpikeTrainClusterSecByData(SpikeTrainClusterByData == cluster_sort_index(i)) = i;
    end
%
end