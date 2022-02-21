function calcSWCTuftEdgeLength(swcFolder, outputpath)
    files = dir([swcFolder, '\*.swc']);
    
    totalLengthNoSoma = zeros(1, length(files));
    totalLength = zeros(1, length(files));
    sumEdge = zeros(1, length(files));
    for i = 1:length(files)
         [gRoi, rootNodeID, selectedROITable] = loadSwcFile(fullfile(files(i).folder,files(i).name), outputpath,false);
         ROI_id_nexus = gRoi.Nodes.ID((gRoi.Nodes.Depth(:,1) == 1));
         ROI_id_soma = gRoi.Nodes.ID((gRoi.Nodes.Depth(:,1) == 0));
         
         idxOutSoma = findedge(gRoi,ROI_id_soma(1),ROI_id_soma(end));
         idxOutNexus1 = findedge(gRoi,ROI_id_soma(end),ROI_id_nexus(1));
         idxOutNexus2 = findedge(gRoi,ROI_id_soma(end),ROI_id_nexus(end));
         
         totalLengthNoSoma(i) = sum(gRoi.Edges.Weight);
         totalLength(i) = sum(gRoi.Edges.Weight);
         sumEdge(i) = sum(gRoi.Edges.Weight) - gRoi.Edges.Weight(idxOutNexus1) - gRoi.Edges.Weight(idxOutNexus2);   
         
         if idxOutSoma ~= 0
            sumEdge(i) = sumEdge(i) - gRoi.Edges.Weight(idxOutSoma);
            totalLengthNoSoma(i) = totalLengthNoSoma(i) - gRoi.Edges.Weight(idxOutSoma);
         end
         
         tabelR(i, 1:4) = {sumEdge(i), totalLengthNoSoma(i), totalLength(i), files(i).name};
         close all;
         
    end
    
    meanE = mean(sumEdge);
    stdE = std(sumEdge);
    meanTotal = mean(totalLength);
    stdETotal = std(totalLength);
    meanNoSoma = mean(totalLengthNoSoma);
    stdENoSoma = std(totalLengthNoSoma);
    save([outputpath, '\Mean_STD_EdgeTuft_All_forallneurons'], 'meanE', 'stdE', 'tabelR', 'totalLengthNoSoma', 'meanNoSoma',...
        'stdENoSoma','totalLength', 'meanTotal', 'stdETotal');
end