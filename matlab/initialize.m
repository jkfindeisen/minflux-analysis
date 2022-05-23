function initialize()
% Simply adds the folder this file resides in and all sub folders to the
% Matlab path. This is needed so that the scripts can access all functions
% that they want to call.

% add the folder this script is in and all subfolders to the matlab path
folder = fileparts(mfilename('fullpath'));
addpath(genpath(folder));

end