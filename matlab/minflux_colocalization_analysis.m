function [N_loc_cluster_reduced, pos_std_cluster_reduced] = minflux_colocalization_analysis(files)
global defaultPath

%% Paramerers for the analysis
% minimal localization dbscan
minpts = 3;
% distance for dbscan
epsilon = 10e-9; % 10 nm!

% Maximal size of std of all localization included in a cluster to be further analyzed for colocalization
coloc_maximal_std_cluster = 3e-8; % 30 nm!
% minimal number of localizations per cluster  to be further analyzed for colocalization
coloc_minimal_N_loc_cluster = 10;



if nargin < 1
    if isempty(defaultPath) || all(defaultPath == 0)
        [files, path] = uigetfile("MultiSelect","on");

    else
        [files, path] = uigetfile(defaultPath, "MultiSelect","on");
    end
    if iscell(files) 
        defaultPath = path;
        files = fullfile(path, files); 
    elseif files == 0
        return
    end
end

%% Read files and clean up of point clouds with a dbscal per track
for  i = 1: numel(files)
    [mat{i}, valid{i}, id{i}, pos{i}] = read_mat_file(files{i}, epsilon,  minpts);
end

%% From here we consider several proteins jointly
%% Join cluster for several point clouds
pos_all = [];
for i = 1:numel(files)
    pos_all = [pos_all; pos{i}];
    id_cluster = minflux_cluster_identify(pos_all, epsilon, minpts);
    uid_cluster = unique(id_cluster(id_cluster > 0));
    N_cluster = numel(uid_cluster);
end

%% statistics on number of localizations in each cluster
stat_cluster = zeros(N_cluster, 7);
for i = 1 : N_cluster
    ix = id_cluster == i;
    stat_cluster(i, :) = [sum(ix), mean(pos_all(ix, :), 1), std(pos_all(ix, :), [], 1)];
end

% Number of localizations per cluster
N_loc_cluster = stat_cluster(:, 1);
% standard deviation of localizations per cluster
pos_std_cluster = sqrt(sum(stat_cluster(:, 5:7).^2, 2));


%% Perform a filtering of cluster based on size and number of localizations
% the rationale is to remove big aggregates that may not be single vesicles
ix = pos_std_cluster < coloc_maximal_std_cluster & N_loc_cluster > coloc_minimal_N_loc_cluster;
fprintf('reduce to %d clusters\n', sum(ix));
id_cluster_reduced = uid_cluster(ix);
N_cluster_reduced = numel(id_cluster_reduced);
pos_std_cluster_reduced = pos_std_cluster(ix);
pos_cluster_red

%% Count colocalizations in each cluster
%works only for 2 proteins for the moment
% split dataset again and count how much in each cluster
n = zeros(N_cluster_reduced, 2);
for i = 1 : N_cluster_reduced
    Nprot1 = size(pos{1},1)
    n1 = sum(id_cluster(1:Nprot1) == id_cluster_reduced(i));
    n2 = sum(id_cluster(Nprot1+1:end) == id_cluster_reduced(i));
    N_coloc_cluster_reduced(i, :) = [n1, n2];
end

id_coloc_objects  = (N_coloc_cluster_reduced(:,1)>0).*(N_coloc_cluster_reduced(:,2)>0);
pos_std_cluster_reduced(id_coloc_objects>0);

%% Plotting
% Plot of projected data in 2D 
close('all')
figure(1);
all_axes = [];
for  i = 1: numel(files)
    ax{i} = subplot(1,numel(files),i);
    all_axes = [all_axes, ax{i}];
    plot(pos{i}(:, 1), pos{i}(:, 2), '.', 'MarkerSize', 4);
    daspect([1,1,1]);
    title(['2D proj. loc. channel ' + num2str(i)] );
end
linkaxes(all_axes);

% Plot 3D localizations
figure(2);
all_axes = [];
for  i = 1: numel(files)
    ax{i} = subplot(1,numel(files),i);
    all_axes = [all_axes, ax{i}];
    plot3(pos{i}(:, 1), pos{i}(:, 2), pos{i}(:, 3), '.', 'MarkerSize', 4);
    daspect([1,1,1]);
    grid on;
    box on;
    view(-44,15);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(['3D loc. channel ' num2str(i)]);
end
linkaxes(all_axes);

figure(3); 
plot(N_loc_cluster, pos_std_cluster, '.');
xlabel('Number of localization per cluster')
ylabel('Std of loc. ~ size of cluster (m)')
title('Size of clusters')

% Plot clusters each cluster has a color
figure(4);
hold on;
for i = 1 : N_cluster
    idx = id_cluster == uid_cluster(i);
    plot3(pos_all(idx, 1), pos_all(idx, 2), pos_all(idx, 3), '.', 'MarkerSize', 4);
end
daspect([1,1,1]);
title(['Joint clusters ' num2str(N_cluster)]);
box on;
grid on;
view(20, 30);

% Plot cleaned reduced clusters
figure(4);
hold on;
for i = 1 : N_cluster_reduced
    idx = id_cluster == uid_cluster(i);
    plot3(pos(idx, 1), pos_all(idx, 2), pos_all(idx, 3), '.', 'MarkerSize', 4);
end
daspect([1,1,1]);
title(['Joint clusters ' num2str(N_cluster)]);
box on;
grid on;
view(20, 30);

figure(4);
hold on;
for i = 1 : numel(files)
    plot3(pos_all(idx, 1), pos_all(idx, 2), pos_all(idx, 3), '.', 'MarkerSize', 4);
end
daspect([1,1,1]);
title(['Joint clusters ' num2str(N_cluster)]);
box on;
grid on;
view(20, 30);


% Plot clusters colored by single protein and joint clusters
figure;
plot(N_coloc_cluster_reduced(:, 1), N_coloc_cluster_reduced(:, 2), '.');
xlabel('channel 1');
ylabel('channel 2');
daspect([1,1,1]);
xlim([0, 1000]);
ylim([0, 1000]);
title(sprintf('ch1 %d, ch2 %d, both %d', ...
   sum(N_coloc_cluster_reduced(:, 1) > 0 & N_coloc_cluster_reduced(:, 1) == 0), ...
    sum(N_coloc_cluster_reduced(:, 1) & N_coloc_cluster_reduced(:, 2) > 0), ...
    sum(N_coloc_cluster_reduced(:, 1)  > 0 &N_coloc_cluster_reduced(:, 2) > 0)));



% hold on;
% for i = 1 : Nc
%     idx = id == uid(i);
%     plot3(pos(idx, 1), pos(idx, 2), pos(idx, 3), '.', 'MarkerSize', 4);
% end
% daspect([1,1,1]);
% title(num2str(Nc));
% box on;
% grid on;
% view(20, 30);
% 
% % statistics on number of localizations in each cluster
% s = zeros(Nc, 7);
% for i = 1 : Nc
%     ix = id == i;
%     s(i, :) = [sum(ix), mean(pos(ix, :), 1), std(pos(ix, :), [], 1)];
% end
% n = s(:, 1);
% sa = sqrt(sum(s(:, 5:7).^2, 2));
% figure; plot(n, sa, '.');
% 
% ix = sa < 3e-8 & n > 10;
% 
% fprintf('reduce to %d clusters\n', sum(ix));
% ids = uid(ix);
% Nc = numel(ids);
% 
% % split dataset again and count how much in each cluster
% n = zeros(Nc, 2);
% for i = 1 : Nc
%     n1 = sum(id(1:N) == ids(i));
%     n2 = sum(id(N+1:end) == ids(i));
%     n(i, :) = [n1, n2];
% end
% 
% figure;
% plot(n(:, 1), n(:, 2), '.');
% xlabel('channel 1');
% ylabel('channel 2');
% daspect([1,1,1]);
% xlim([0, 1000]);
% ylim([0, 1000]);
% title(sprintf('ch1 %d, ch2 %d, both %d', sum(n(:, 1) > 0 & n(:, 2) == 0), sum(n(:, 1) == 0 & n(:, 2) > 0), sum(n(:, 1) > 0 & n(:, 2) > 0)));

end

function [mat, valid, id, pos] = read_mat_file(filename,  epsilon, minpts)
    % read the file 
    mat = load(filename);
    % non combined localizations
    valid = mat.vld(:);
    pos = squeeze(mat.itr.loc(valid, end, :));
    id = mat.tid(valid);
    id = id(:); % to have everything Nx1

    % filter according to dbscan
    pos = minflux_filter(pos, id,  epsilon, minpts);

end

    function plot_point_clouds()
% % co-registration of point clouds
% figure();
% ax1 = subplot(1,2,1);
% plot(pos1(:, 1), pos1(:, 2), '.', 'MarkerSize', 4);
% daspect([1,1,1]);
% ax2 = subplot(1,2,2);
% plot(pos2(:, 1), pos2(:, 2), '.', 'MarkerSize', 4);
% daspect([1,1,1]);
% linkaxes([ax1, ax2]);
% 
% fig = figure();
% fig.Position = [100, 100, 1000, 1000];
% ax1 = subplot(1,2,1);
% plot3(pos1(:, 1), pos1(:, 2), pos1(:, 3), '.', 'MarkerSize', 4);
% daspect([1,1,1]);
% grid on;
% box on;
% view(-44,15);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title('channel 1');
% ax2 = subplot(1,2,2);
% plot3(pos2(:, 1), pos2(:, 2), pos2(:, 3), '.', 'MarkerSize', 4);
% daspect([1,1,1]);
% grid on;
% box on;
% view(-44,15);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title('channel 2');
% 
% fig = figure();
% fig.Position = [100, 100, 1000, 1000];
% hold on;
% plot3(pos1(:, 1), pos1(:, 2), pos1(:, 3), 'b.', 'MarkerSize', 4);
% plot3(pos2(:, 1), pos2(:, 2), pos2(:, 3), 'r.', 'MarkerSize', 4);
% daspect([1,1,1]);
% grid on;
% box on;
% view(-44,15);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title('channel 2');
% 
end