function DeleteNanValuesFromHeatMap(figLocation, figName)
    f = openfig(fullfile(figLocation, figName));
    index = ~isnan(f.Children(2).Children(1).CData(1, :));
    f.Children(2).Children(1).CData = f.Children(2).Children(1).CData(index, index);
    f.Children(2).Children(1).AlphaData = f.Children(2).Children(1).AlphaData(index, index);
    mysave(f, fullfile(figLocation, sprintf('Fix_%s', figName)));
end