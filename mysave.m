function mysave(f,filename, ext)
if isempty(f)
    f=gcf;
end
pth=fileparts(filename);
mkNewFolder(pth);
if ~exist('ext','var') || isempty(ext)
saveas(f,[filename '.fig'],'fig');
saveas(f,[filename '.tif'],'tif');
% print( f, '-depsc', [filename '.eps']);
else
    switch ext
        case 'tif'
           saveas(f,[filename '.tif'],'tif');
        case 'fig'
            saveas(f,[filename '.fig'],'fig');
        case 'eps'
           print( f, '-depsc', [filename '.eps']);
    end
end
        
% print(f,[filename '.pdf'],'-dpdf');
