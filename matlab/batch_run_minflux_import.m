[files, path] = uigetfile("MultiSelect","on", {'*.mat'; '*.*'},...
                          'File Selector');
outdir = uigetdir(path, "Specify output folder");
for i = 1: length(files)
   minflux_import(fullfile(path, files{i}), outdir);
end
