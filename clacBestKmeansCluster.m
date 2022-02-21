function clacBestKmeansCluster(EventsListLocation, Cutoff)
   load(EventsListLocation, 'allEventsTable');
   
   sumDR = zeros(1, 10);
   for k = 1:10
       [~,~,sumd,~] = kmeans(allEventsTable.H, k, 'Replicates',5, 'MaxIter', 500);
       sumDR(k) = sum(sumd);
   end
   
   Var=sumDR(1:end-1)-sumDR(2:end); %calculate %variance explained
   PC=cumsum(Var)/(sumDR(1)-sumDR(end));
   
   [r,~]=find(PC'>Cutoff); %find the best index
   K=r(1,1);
   
   figure;
   plot(PC,'b*--');
   title(['Best K = ', num2str(K)]);
 
end