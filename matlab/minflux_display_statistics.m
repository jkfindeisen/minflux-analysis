function figs = minflux_display_statistics(stat)
% Creates a figure with some overview characteristics of Minflux data

% TODO efo, efc

assert(nargin == 1);
figs = [];

%% parameters

% plot information
lw = 2; % linewidth of median line

% Colors for plot
cf = [0.3010 0.7450 0.9330]; % histogram face color
ce = [0 0.4470 0.7410]; % histogram edge color
cm = [0.6350 0.0780 0.1840]; % plot median line color

%% 2D xy histogram rendering of drift corrected positions

% x, y coordinates (in nm)
x = stat.dpos(:, 1) / 1e-9;
y = stat.dpos(:, 2) / 1e-9;

% reasonable boundaries (rounded to next 100nm)
a = 0.02;
R = 100;
Rx = quantile(x, [a, 1-a]);
Rx = [floor(Rx(1)/R)*R, ceil(Rx(2)/R)*R];
Ry = quantile(y, [a, 1-a]);
Ry = [floor(Ry(1)/R)*R, ceil(Ry(2)/R)*R];

% 2D render of x,y (histogram, 1nm pixel size)
sxy = 1;
[h, xh, yh] = render_xy(x, y, sxy, sxy, Rx, Ry);

% display
fig = figure(345);
clf('reset');
fig.WindowState = 'maximized';
im = imagesc(xh, yh, h.');
im.Parent.YDir = 'normal';
axis image;
colormap(hot);
caxis([0, max(h(:))* 0.5]);
hold on;
xlabel('x (nm)');
xlim(Rx);
ylabel('y (nm)');
ylim(Ry);
title(sprintf('%d events recorded in %.1f min (%.1f / s)', numel(unique(stat.id)), stat.t(end)/60, numel(unique(stat.id)) / stat.t(end)));

figs = [figs; fig];

%% interesting histograms
fig = figure(346);
fig.Position = [100, 100, 1600, 1000];
clf('reset');

% background frequency - plot

subplot(3, 5, 1);
hold on;
m = median(stat.fbg);

g = 0:2e3:6e4;
histogram(stat.fbg, g, 'FaceColor', cf, 'EdgeColor',ce);
xlim(g([1,end]));

plot(m*[1,1],ylim(),'LineWidth',3,'Color',cm);

decorate('background signal (Hz)', 'occurence', sprintf('median: %.0f kHz', m / 1e3));


%% CFR center-frequency-ratio - plot
subplot(3, 5, 2);
hold on
m = median(stat.cfr);

g = -0.5:0.05:1.5;
histogram(stat.cfr,g,'FaceColor',	cf, 'EdgeColor',ce);
xlim(g([1,end]));

plot(m*[1,1],ylim(),'LineWidth',lw,'Color',cm);

decorate('center-frequency-ratio', 'occurence', sprintf('median: %.2f', m));

%% time between events (t_scan)
subplot(3, 5, 3);
hold on;

m = median(stat.t_scan);

g = 0:0.2:10;
histogram(stat.t_scan,g,'FaceColor',	cf, 'EdgeColor', ce);
xlim(g([1,end]));

plot(m*[1, 1],ylim(),'LineWidth',lw,'Color',cm);
decorate('time elapsed between events t_{scan} (s)', 'occurence', sprintf('median t_{scan}= %.2f s', m));

%% time within events (t_loc)
subplot(3, 5, 4);
hold on;

m = median(stat.t_loc);

g = 0:0.1:6;
histogram(stat.t_loc,g,'FaceColor', cf, 'EdgeColor', ce);
xlim(g([1,end]));

plot(m*[1, 1],ylim(),'LineWidth',lw,'Color',cm);
decorate('time elapsed within events t_{loc} (s)', 'occurence', sprintf('median t_{loc}= %.2f s', m));

%% number of iteration repeats for each trace
subplot(3,5,5);
hold on;

m = mean(stat.combined.n);

g = 0:1:50;
h = histogram(stat.combined.n, g, 'FaceColor', cf, 'EdgeColor', ce);
h.Parent.YScale = 'log';
xlim(g([1,end]));
plot(m*[1, 1],ylim(),'LineWidth',lw,'Color',cm);
decorate('number iteration repeats / trace', 'occurence', sprintf('mean n=%.2f', m));


%% plot sigma x, y, z, r
T = 5;
sigmas = stat.combined.std_xyz(stat.combined.n >= T, :); % all events with at least T localizations
sigmas = [sigmas, sqrt(sum(sigmas(:, 1:2).^2, 2)/2)]; % adds sigma_r
sigmas = sigmas / 1e-9; % now in nm
med = median(sigmas, 1);

g = 0:0.5:20; % limit of x-range for histogram of sigma_r
labels = {'x', 'y', 'z', 'r'};
for i = 1 : 4
    subplot(3, 5, 5+i);
    histogram(sigmas(:, i),g,'FaceColor',	cf, 'EdgeColor',ce);
    hold on;
    plot(med(i)*[1,1],ylim(),'LineWidth',lw,'Color',cm);
    xlim(g([1,end]));
    decorate(sprintf('\\sigma_%s (nm)', labels{i}), 'occurence', sprintf('median \\sigma_%s = %.1f nm ', labels{i}, med(i)));
end

%% std error of combined localizations
s = stat.combined.ste_xyz / 1e-9; % in nm
sr = sqrt((s(:,1).^2+s(:,2).^2)/2);
sz = s(:, 3);

g = 0:0.5:10; % limit of x-range for histogram of sigma_r

subplot(3,5,11);
histogram(sr,g,'FaceColor', cf, 'EdgeColor',ce);
hold on;
plot(median(sr)*[1,1],ylim(),'LineWidth',lw,'Color',cm);
xlim(g([1,end]));
decorate('\epsilon_r (nm)', 'occurence', sprintf('median \\epsilon_r = %.1f nm', median(sr)));

subplot(3,5,12);
histogram(sz,g,'FaceColor', cf, 'EdgeColor',ce);
hold on;
plot(median(sz)*[1,1],ylim(),'LineWidth',lw,'Color',cm);
xlim(g([1,end]));
decorate('\epsilon_z (nm)', 'occurence', sprintf('median \\epsilon_z = %.1f nm', median(sz)));


%% FRC
frc = stat.frc;
qi = frc.qi{1};
ci = mean(cat(2, frc.ci{:}), 2);
resolution = mean(frc.resolution);
subplot(3,5,13);
plot(qi, ci);
decorate('k', 'FRC', sprintf('est. res. %.1f nm', resolution/1e-9));

frc = stat.frc_combined;
qi = frc.qi{1};
ci = mean(cat(2, frc.ci{:}), 2);
resolution = mean(frc.resolution);
subplot(3,5,14);
plot(qi, ci);
decorate('k', 'FRC of combined', sprintf('est. res. %.1f nm', resolution/1e-9));

figs = [figs; fig];

%% drift
d = stat.drift;
subplot(3,5,15);
hold on;
for i = 2 : size(d, 2)
    plot(d(:, 1), d(:, i)/1e-9, 'DisplayName', labels{i-1});
end
decorate('time (s)', 'Est. drift (nm)');
legend('x', 'y', 'z', 'Location', 'best');

end

function decorate(labelx, labely, plot_title)
% some often used  functionality together

if nargin > 2
    title(plot_title);
end

xlabel(labelx);
ylabel(labely);
grid on;
box on;
pbaspect([1 1 1]);
end