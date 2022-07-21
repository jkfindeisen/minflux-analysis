function minflux_colocalization_analysis()
% Minflux colocalization analysis

% file 1
file = 'C:\Users\ssamban1\Desktop\Minflux\ATG9-Syp\10-05-22\ROI3\220510-170636_3Dminflux_Syp_Paint_2nM_P1_16%640_PH06.mat';
m1 = load(file);

% non combined localizations
valid = m1.vld(:);
pos1 = squeeze(m1.itr.loc(valid, end, :));
id1 = m1.tid(valid);
id1 = id1(:); % to have everything Nx1
pos1 = minflux_filter(pos1, id1);

% file 2
file = 'C:\Users\ssamban1\Desktop\Minflux\ATG9-Syp\10-05-22\ROI3\220510-175507_3Dminflux_ATG_Paint_2nM_P2_16%640_PH06.mat';
m2 = load(file);

% non combined localizations
valid = m2.vld(:);
pos2 = squeeze(m2.itr.loc(valid, end, :));
id2 = m2.tid(valid);
id2 = id2(:); % to have everything Nx1
pos2 = minflux_filter(pos2, id2);

% co-registration of point clouds
figure();
ax1 = subplot(1,2,1);
plot(pos1(:, 1), pos1(:, 2), '.', 'MarkerSize', 4);
daspect([1,1,1]);
ax2 = subplot(1,2,2);
plot(pos2(:, 1), pos2(:, 2), '.', 'MarkerSize', 4);
daspect([1,1,1]);
linkaxes([ax1, ax2]);

fig = figure();
fig.Position = [100, 100, 1000, 1000];
ax1 = subplot(1,2,1);
plot3(pos1(:, 1), pos1(:, 2), pos1(:, 3), '.', 'MarkerSize', 4);
daspect([1,1,1]);
grid on;
box on;
view(-44,15);
xlabel('x');
ylabel('y');
zlabel('z');
title('channel 1');
ax2 = subplot(1,2,2);
plot3(pos2(:, 1), pos2(:, 2), pos2(:, 3), '.', 'MarkerSize', 4);
daspect([1,1,1]);
grid on;
box on;
view(-44,15);
xlabel('x');
ylabel('y');
zlabel('z');
title('channel 2');

fig = figure();
fig.Position = [100, 100, 1000, 1000];
hold on;
plot3(pos1(:, 1), pos1(:, 2), pos1(:, 3), 'b.', 'MarkerSize', 4);
plot3(pos2(:, 1), pos2(:, 2), pos2(:, 3), 'r.', 'MarkerSize', 4);
daspect([1,1,1]);
grid on;
box on;
view(-44,15);
xlabel('x');
ylabel('y');
zlabel('z');
title('channel 2');

% cluster identification on both point clouds
N = size(pos1, 1);
pos = [pos1; pos2];
id = minflux_cluster_identify(pos);

uid = unique(id(id > 0));
Nc = numel(uid);

fig = figure();
hold on;
for i = 1 : Nc
    idx = id == uid(i);
    plot3(pos(idx, 1), pos(idx, 2), pos(idx, 3), '.', 'MarkerSize', 4);
end
daspect([1,1,1]);
title(num2str(Nc));
box on;
grid on;
view(20, 30);

% statistics on number of localizations in each cluster
s = zeros(Nc, 7);
for i = 1 : Nc
    ix = id == i;
    s(i, :) = [sum(ix), mean(pos(ix, :), 1), std(pos(ix, :), [], 1)];
end
n = s(:, 1);
sa = sqrt(sum(s(:, 5:7).^2, 2));
figure; plot(n, sa, '.');

ix = sa < 3e-8 & n > 10;

fprintf('reduce to %d clusters\n', sum(ix));
ids = uid(ix);
Nc = numel(ids);

% split dataset again and count how much in each cluster
n = zeros(Nc, 2);
for i = 1 : Nc
    n1 = sum(id(1:N) == ids(i));
    n2 = sum(id(N+1:end) == ids(i));
    n(i, :) = [n1, n2];
end

figure;
plot(n(:, 1), n(:, 2), '.');
xlabel('channel 1');
ylabel('channel 2');
daspect([1,1,1]);
xlim([0, 1000]);
ylim([0, 1000]);
title(sprintf('ch1 %d, ch2 %d, both %d', sum(n(:, 1) > 0 & n(:, 2) == 0), sum(n(:, 1) == 0 & n(:, 2) > 0), sum(n(:, 1) > 0 & n(:, 2) > 0)));

end

