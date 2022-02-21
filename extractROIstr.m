function roiNum = extractROIstr(fullstr)

roiNum=sscanf(fullstr, 'ROI:%d Z:1');
if isempty(roiNum)
    roiNum=sscanf(fullstr, 'ROI:Z:1:%d');
end
if isempty(roiNum)
    roiNum=sscanf(fullstr, 'ROI:Z:2:%d');
end
if isempty(roiNum)
    roiNum=sscanf(fullstr, 'ROI:Z:3:%3d:NZ-1');
end
if isempty(roiNum)
    roiNum=sscanf(fullstr, 'ROI:Z:4:%3d:NZ-1');
end
if isempty(roiNum)
    roiNum=sscanf(fullstr, 'ROI %2d:NZ-1');
end
if isempty(roiNum)
    roiNum=sscanf(fullstr, 'AROI:Z1:%d');
end
if isempty(roiNum)
    roiNum=sscanf(fullstr, 'ROI:Z:%d');
end
if isempty(roiNum)
    roiNum=sscanf(fullstr, 'roi%05d');
end
if isempty(roiNum)
    roiNum=sscanf(fullstr, 'apic%05d');
end
if isempty(roiNum)
    roiNum=0;
    warning('Unrecognized neuron name');
end

