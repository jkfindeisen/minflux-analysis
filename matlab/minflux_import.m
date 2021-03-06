function ROI = minflux_import(file, output_folder, ROI)
% Reads an exported Minflux data file from AI and
%   - exports the valid localizations into a csv file
%   - renders the data set in 2/3D as histogram/Gaussian
%   - exports as tiff image
%   - estimates a few interesting characteristics
%
% Input parameters
%   file    ".mat" file with exported Minflux localizations from an AI
%           Minflux microscope measurement
%   output_folder   where to put the imported data (if not given or empty
%                   takes the folder in which the input file resides
%   ROI     a 2x3 matrix with the first row xyz coordinates of the upper,
%           left corner of the ROI and the second row the lower, right
%           corner
%
% Notes
%
% - If you have multiple channels and you need renderings of the same ROI,
%   just call minflux_import without a ROI on one of them, store the
%   returned ROI and use it (as input parameter) for all other channels.

close all;

if ~exist('drift_correction_time_windows_2D.m', 'file')
    initialize();
end

%% input parameter checking
global defaultPath
if nargin < 1
    if isempty(defaultPath) || all(defaultPath == 0)
        [file, path] = uigetfile();
    else
        [file, path] = uigetfile(defaultPath);
    end
    if file == 0
        return
    end
    defaultPath = path;
    file = fullfile(path, file);
end


fprintf('work on file %s\n', file);
[file_path , base_name] = fileparts(file);
if nargin < 2
   output_folder =uigetdir(file_path, "Specify output folder");
end

if ~endsWith(output_folder, filesep)
    output_folder = [output_folder, filesep]; % add slash if needed
end
if nargin < 3 || isempty(ROI)
    ROI = [];
end

%% parameters to adjust the script (customize here)
correct_z_position = true;      % AI Minflux micrsocopes export z positions wrongly (beginning of 2022), if true, corrects for that
correct_z_position_factor = 0.7;

save_statistics_to_file = true; % if true, stores a figure snapshot of the statistics figures
export_localizations_to_csv = true; % if true, exports localizations to csv file format
export_renderings = true; % if true, exports renderings to various formats
gaussian_render_fwhm = 5e-9; % FWHM of Gaussian peak render (m)


%% load data and if needed, correct z position
minflux = load(file);

if correct_z_position
    minflux.itr.loc(:, :, 3) = minflux.itr.loc(:, :, 3) * correct_z_position_factor;
end

%% calculate statistics and display statistics overview sheet
% calculate statistics
fprintf('calculate Minflux statistics (including intrinsic drift correction)\n');
stat = minflux_calculate_statistics(minflux);

% display statistics
figs = minflux_display_statistics(stat);

if save_statistics_to_file
    exportgraphics(figs(2), [output_folder, base_name, '.statistics.png']);
end

%% export localizations to csv file

if export_localizations_to_csv
    fprintf('export localizations to csv files\n');
    
    % not combined
    fprintf(' localizations\n');
    t = array2table([stat.pos, stat.t, stat.id], 'VariableNames', {'x (m)', 'y (m)', 'z (m)', 't (s)', 'trace id'});
    writetable(t, [output_folder, base_name, '.localizations.csv']);
    
    % not combined, drift corrected
    fprintf(' events\n');
    t = array2table([stat.dpos, stat.t, stat.id], 'VariableNames', {'x (m)', 'y (m)', 'z (m)', 't (s)', 'trace id'});
    writetable(t, [output_folder, base_name, '.localizations-driftcorrected.csv']);
    
    % combined, drift corrected
    fprintf(' drift corrected events\n');
    c = stat.combined;
    t = array2table([c.pos, c.t, c.n, c.std_xyz], 'VariableNames', {'x (m)', 'y (m)', 'z (m)', 't (s)', 'n loc/event', 'sigma x (m)', 'sigma y (m)', 'sigma z (m)'});
    writetable(t, [output_folder, base_name, '.localizations-driftcorrected-combined.csv']);
end

%% export rendered image (tiff, obf, ome-tiff)
sxyz = 2e-9; % pixel size for renderings (m)

if export_renderings
    fprintf('export rendered files to tiff/obf/ome-tiff (can take a while)\n');
    
    % get reasonable ranges if ROI not given
    if isempty(ROI)
        a = 0.01;
        ROI = quantile(stat.pos, [a, 1-a]);
    end
    
    % loop over localizations, drift-corrected localizations, combined
    h1 = render(stat.pos, ROI, stat.is3D, sxyz, gaussian_render_fwhm);
    h2 = render(stat.dpos, ROI, stat.is3D, sxyz, gaussian_render_fwhm);
    h3 = render(stat.combined.pos, ROI, stat.is3D, sxyz, gaussian_render_fwhm);
    
    % and write out
    path = [output_folder, base_name];
    fprintf(' localizations\n');
    export(h1, '.localizations', stat.is3D, path, sxyz);
    fprintf(' events\n');
    export(h2, '.localizations-driftcorrected', stat.is3D, path, sxyz);
    fprintf(' drift corrected events\n');
    export(h3, '.localizations-driftcorrected-combined', stat.is3D, path, sxyz);
    
end

end

function export(h, name, is3D, path, sxyz)
asOme = 0;
% export a set of rendering (histogram, gaussian renderings) to normal tiff
% (if 2D), to obf (if specmy available), to ometiff (if ome available)

% export to tiff (in 3D case a multipage tiff)
resolution = sxyz / 1e-9 * [1,1];
if ~is3D
    h1 = h{1}; % histogram rendering, do not scale
    imwrite(h1, [path, name, '.histogram.tiff'], 'Resolution', resolution);
    
    h2 = h{2}; % gaussian rendering
    imwrite(uint16(h2 / max(h2(:)) * 2^16), [path, name, '.gaussian-rendering.tiff'], 'Resolution', resolution);
else
    h1 = h{1}; % histogram rendering, do not scale
    file = [path, name, '.histogram.tiff'];
    imwrite(h1(:, :, 1), file, 'Resolution', resolution);
    for i = 2 : size(h1, 3)
        imwrite(h1(:, :, i), file, 'WriteMode', 'append', 'Resolution', resolution);
    end
    
    h2 = h{2}; % gaussian rendering
    h2 = uint16(h2 / max(h2(:)) * 2^16);
    file = [path, name, '.gaussian-rendering.tiff'];
    imwrite(h2(:, :, 1), file, 'Resolution', resolution);
    for i = 2 : size(h2, 3)
        imwrite(h2(:, :, i), file, 'WriteMode', 'append', 'Resolution', resolution);
    end
end

% export to obf (AI native format)
try
    h1 = h{1};
    s1 = specmx.Stack(h1);
    s1.Name = [name, '.histogram'];
    s1.PixelLengths = sxyz * ones([1, numel(size(h1))]);
    h2 = h{2};
    s2 = specmx.Stack(h2);
    s2.Name = [name, '.gaussian-rendering'];
    s2.PixelLengths = sxyz * ones([1, numel(size(h2))]);
    specmx.File.write_file([path, name, '.obf'], {s1, s2});
catch ex
    % could not write obf
    fprintf('%s - %s\n', ex.identifier, ex.message);
end

if asOme
    % save as ometiff
    try
        h1 = h{1};
        h2 = h{2};
        bfCheckJavaPath();
        pixelSize = ome.units.quantity.Length(java.lang.Double(sxyz / 1e-6), ome.units.UNITS.MICROMETER);
        meta = createMinimalOMEXMLMetadata(h1); % default dimension order XYZCT
        meta.setPixelsPhysicalSizeX(pixelSize, 0);
        meta.setPixelsPhysicalSizeY(pixelSize, 0);
        meta.setPixelsPhysicalSizeZ(pixelSize, 0);
        nameOriginal = name;
        name = [nameOriginal, '.histogram'];
        bfsave(h1, [path, name, '.ome.tiff'], 'metadata', meta);
        meta = createMinimalOMEXMLMetadata(h2); % default dimension order XYZCT
        meta.setPixelsPhysicalSizeX(pixelSize, 0);
        meta.setPixelsPhysicalSizeY(pixelSize, 0);
        meta.setPixelsPhysicalSizeZ(pixelSize, 0);
        name = [nameOriginal, '.gaussian-rendering'];
        bfsave(h2, [path, name, '.ome.tiff'], 'metadata', meta);
    catch ex
        % could not write ome tiff
        fprintf('%s - %s\n', ex.identifier, ex.message);
    end
end
end

function h = render(pos, R, is3D, sxyz, gaussian_render_fwhm)
% render as histogram and as Gaussian peaks
% h{1} is uint8

assert(gaussian_render_fwhm >= 2*sxyz);

opt1 = struct('type', 'histogram');
opt2 = struct('type', 'fixed_gaussian', 'fwhm', gaussian_render_fwhm);

if is3D
    % 3D rendering needed
    h1 = render_xyz(pos(:, 1), pos(:, 2), pos(:, 3), sxyz, sxyz, sxyz, R(:, 1), R(:, 2), R(:, 3), opt1);
    h2 = render_xyz(pos(:, 1), pos(:, 2), pos(:, 3), sxyz, sxyz, sxyz, R(:, 1), R(:, 2), R(:, 3), opt2);
else
    % 2D rendering
    h1 = render_xy(pos(:, 1), pos(:, 2), sxyz, sxyz, R(:, 1), R(:, 2), opt1);
    h2 = render_xy(pos(:, 1), pos(:, 2), sxyz, sxyz, R(:, 1), R(:, 2), opt2);
end

h = {uint8(h1), h2}; % histogram is uint8

end