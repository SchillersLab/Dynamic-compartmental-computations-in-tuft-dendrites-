function plot2Dembedding(examinedInds, outputPath, eventsStr, labelsLUT, labels,embedding, clrs, Method, legendLoc, labelsFontSz)


classes = unique(labels);
figure;
for ci = 1:length(classes)
    scatter(embedding(labels==classes(ci),1)  ,embedding(labels==classes(ci),2),...
        'o', 'MarkerFaceColor', clrs(ci, :), 'MarkerEdgeColor', clrs(ci, :));
    hold all;
end
xlabel('\psi_1', 'FontSize', labelsFontSz), ylabel('\psi_2', 'FontSize', labelsFontSz);
l=legend(labelsLUT, 'Location',legendLoc);

set(gca, 'Box','off');
set(l, 'FontSize',labelsFontSz);
a=get(gcf,'Children');
for ai=1:length(a)
    if strcmp(a(ai).Type, 'Axes')
        setAxisFontSz(a(ai), labelsFontSz);
    end
end
mysave(gcf, fullfile(outputPath, [ Method '2D' eventsStr]));


strs = cellstr(num2str(examinedInds(:)));
text(embedding(:,1),embedding(:,2),strs,'VerticalAlignment','bottom',...
    'HorizontalAlignment','right');
mysave(gcf, fullfile(outputPath, [ Method '2D' eventsStr 'withNumbers']));
end
