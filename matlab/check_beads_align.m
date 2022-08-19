function main
fname = '220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b'
indata = load(strcat('C:\Users\apoliti\Desktop\mfluxtest\analysis\', fname , '\' , fname ,'.mat'))
fn = fieldnames(indata);
scatter_plot(1, indata, 'lnc', {'blue', 'red'}', 0.15)
scatter_plot(2, indata, 'loc', {'blue', 'red'}', 0.15)
scatter_plot(3, indata, 'ltr', {'blue', 'red'}', 0.15)
scatter_plot(4, indata, 'lre', {'blue', 'red'}', 0.15)

function scatter_plot(fig, indata, pos, col, al)
    fn = fieldnames(indata);
    figure(fig)
    clf;
    hold;
    for i=1:numel(fn)
        scatter3(indata.(fn{i}).(pos)(:,1), indata.(fn{i}).(pos)(:,2), indata.(fn{i}).(pos)(:,3), ...
            4, 'filled', 'MarkerFaceColor',col{i});
    end
    alpha(al)
    daspect([1,1,1]);
    grid on;
    box on;
    view(-44,15);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(pos);
end

refBeads = load(strcat('C:\Users\apoliti\Desktop\mfluxtest\analysis\', fname , '\' , fname ,'_ref.mat'))
nrBeads = 2
figure(5)
ax1 = subplot(1,3,1);
plot3(refBeads.lnc(:, 1), refBeads.lnc(:, 2), refBeads.lnc(:, 3), '.', 'MarkerSize', 10);
daspect([1,1,1]);
grid on;
box on;
view(-44,15);
xlabel('x');
ylabel('y');
zlabel('z');
title(['3D lnc channel']);

ax2 =  subplot(1,3,2);
plot3(refBeads.ltr(:, 1), refBeads.ltr(:, 2), refBeads.ltr(:, 3), '.', 'MarkerSize', 10);
daspect([1,1,1]);
grid on;
box on;
view(-44,15);
xlabel('x');
ylabel('y');
zlabel('z');
title(['3D ltr channel']);

ax3 =  subplot(1,3,3);
plot3(refBeads.lre(:, 1), refBeads.lre(:, 2), refBeads.lre(:, 3), '.', 'MarkerSize', 10);
daspect([1,1,1]);
grid on;
box on;
view(-44,15);
xlabel('x');
ylabel('y');
zlabel('z');
title(['3D lre channel']);
%linkaxes([ax1, ax2, ax3])


end


% %%
% refBeads = load('C:\Users\apoliti\Desktop\mfluxtest\220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b\refBeads.mat')
