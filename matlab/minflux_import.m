function minflux_import(file, output_folder)
% Reads a Minflux data file from AI and
%   - exports the valid localizations into a csv file
%   - renders the data set in 2/3D as histogram/Gaussian
%   - exports as tiff image
%   - estimates a few interesting characteristics
%
% Optionally, for PAINT measurments, split colors by time.

close all;

if ~exist('drift_correction_time_windows_2D.m', 'file')
    initialize();
end

%% default parameters
assert(nargin >= 1, 'Please provide a file path to a *.mat file of exported localizations.');
[file_path , base_name] = fileparts(file);
if nargin < 2
    output_folder = file_path; % use input file folder as output folder if no output folder is given explicitely
end
if ~endsWith(output_folder, filesep)
    output_folder = [output_folder, filesep]; % add slash if needed
end

correct_z_position = true;      % AI Minflux micrsocopes export z positions wrongly (beginning of 2022), if true, corrects for that
correct_z_position_factor = 0.7;

save_statistics_to_file = true; % if true, stores a figure snapshot of the statistics figures
export_localizations_to_csv = true; % if true, exports localizations to csv file format
export_renderings = true; % if true, exports renderings to various formats


%% load data and if needed, correct z position
minflux = load(file);

if correct_z_position
    minflux.itr.loc(:, :, 3) = minflux.itr.loc(:, :, 3) * correct_z_position_factor;
end

%% calculate statistics and display statistics overview sheet
% calculate statistics
stat = minflux_calculate_statistics(minflux);

% display statistics
figs = minflux_display_statistics(stat);

if save_statistics_to_file
    exportgraphics(figs(2), [output_folder, base_name, '.statistics.png']);
end

%% export localizations to csv file

if export_localizations_to_csv
    
    % not combined
    t = array2table([stat.pos, stat.t, stat.id], 'VariableNames', {'x (m)', 'y (m)', 'z (m)', 't (s)', 'trace id'});
    writetable(t, [output_folder, base_name, '.localizations.csv']);
    
    % not combined, drift corrected
    t = array2table([stat.dpos, stat.t, stat.id], 'VariableNames', {'x (m)', 'y (m)', 'z (m)', 't (s)', 'trace id'});
    writetable(t, [output_folder, base_name, '.localizations-driftcorrected.csv']);
    
    % combined, drift corrected
    c = stat.combined;
    t = array2table([c.pos, c.t, c.n, c.std_xyz], 'VariableNames', {'x (m)', 'y (m)', 'z (m)', 't (s)', 'n loc/event', 'sigma x (m)', 'sigma y (m)', 'sigma z (m)'});
    writetable(t, [output_folder, base_name, '.localizations-driftcorrected-combined.csv']);
end

%% export rendered image (tiff, obf, ome-tiff)

if export_renderings
    
    % get reasonable ranges
    a = 0.01;
    R = quantile(stat.pos, [a, 1-a]);
    
    % loop over localizations, drift-corrected localizations, combined
    h1 = render(stat.pos, R, stat.is3D);
    h2 = render(stat.dpos, R, stat.is3D);
    h3 = render(stat.combined.pos, R, stat.is3D);
    
    % and write out
    path = [output_folder, base_name];
    export(h1, '.localizations', stat.is3D, path);
    export(h2, '.localizations-driftcorrected', stat.is3D, path);
    export(h3, '.localizations-driftcorrected-combined', stat.is3D, path);
    
end

end

function export(h, name, is3D, path)
% export a set of rendering (histogram, gaussian renderings) to normal tiff
% (if 2D), to obf (if specmy available), to ometiff (if ome available)

% export to tiff (in 3D case a multipage tiff)
if ~is3D
    h1 = h{1}; % histogram rendering, do not scale
    imwrite(h1, [path, name, '.histogram.tiff']);
    
    h2 = h{2}; % gaussian rendering
    imwrite(uint16(h2 / max(h2(:)) * 2^16), [path, name, '.gaussian-rendering.tiff']);
else
    h1 = h{1}; % histogram rendering, do not scale
    file = [path, name, '.histogram.tiff'];
    imwrite(h1(:, :, 1), file);
    for i = 2 : size(h1, 3)
        imwrite(h1(:, :, i), file, 'WriteMode', 'append');
    end
    
    h2 = h{2}; % gaussian rendering
    h2 = uint16(h2 / max(h2(:)) * 2^16);
    file = [path, name, '.gaussian-rendering.tiff'];
    imwrite(h2(:, :, 1), file);
    for i = 2 : size(h2, 3)
        imwrite(h2(:, :, i), file, 'WriteMode', 'append');
    end
end

sxyz = 2e-9;

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
    % something went wrong
    warning(ex);
end

% save as ometiff
try
    h1 = h{1};
    h2 = h{2};
    bfCheckJavaPath();
    pixelSize = ome.units.quantity.Length(java.lang.Double(sxyz / 1e-6), ome.units.UNITS.MICROMETER);
    meta = createMinimalOMEXMLMetadata(h11); % default dimension order XYZCT
    meta.setPixelsPhysicalSizeX(pixelSize, 0);
    meta.setPixelsPhysicalSizeY(pixelSize, 0);
    meta.setPixelsPhysicalSizeZ(pixelSize, 0);
    name = [name, '.histogram'];
    bfsave(h1, [path, name, '.histogram.ome.tiff'], 'metadata', meta);
    name = [name, '.gaussian-rendering'];
    bfsave(h2, [path, name, '.histogram.ome.tiff'], 'metadata', meta);
catch
end

end

function h = render(pos, R, is3D)
% render as histogram and as Gaussian peaks
% h{1} is uint8

sxyz = 2e-9; % set above
opt1 = struct('type', 'histogram');
opt2 = struct('type', 'fixed_gaussian', 'fwhm', 2*sxyz);

if is3D
    % 3D rendering needed
    h1 = render_xyz(pos(:, 1), pos(:, 2), pos(:, 3), sxyz, sxyz, sxyz, R(:, 1), R(:, 2), R(:, 3), opt1);
    h2 = render_xyz(pos(:, 1), pos(:, 2), pos(:, 3), sxyz, sxyz, sxyz, R(:, 1), R(:, 2), R(:, 3), opt2);
else
    % 2D rendering
    h1 = render_xy(pos(:, 1), pos(:, 2), sxyz, sxyz, R(:, 1), R(:, 2), opt1);
    h2 = render_xy(pos(:, 1), pos(:, 2), sxyz, sxyz, R(:, 1), R(:, 2), opt2);
end

h = {uint8(h1), h2};

end