FolderName = './results_XPCS/power_series/';   % Your destination folder
mkdir(FolderName);

save([FolderName 'results_allscans.mat'],'-v7.3');

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, [FolderName FigName '.fig']);
end
