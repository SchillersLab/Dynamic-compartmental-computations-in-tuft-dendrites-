function PlotData(dataFileLocation, toneTime)
    data_loaded = load(dataFileLocation);
    
    % Create figure
    figure1 = figure;

    % Create axes
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    
    xlim(axes1,[1 size(data_loaded.sLFP,2)]);
    ylim(axes1,[1 size(data_loaded.sLFP,1)]);
    
    % Create image
    imagesc(data_loaded.sLFP,'Parent',axes1);
    colormap jet;
   
    % Create ylabel
    ylabel({'Cell'});

    % Create xlabel
    xlabel({'Time'});


    % Set the remaining axes properties
    set(axes1,'Layer','top');
    % Create colorbar
    colorbar('peer',axes1);

    line([toneTime toneTime], get(gca, 'YLim'), 'Color','k','LineWidth',2, 'LineStyle', ':');
    
end