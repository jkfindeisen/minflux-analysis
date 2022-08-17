function main
indata = load('C:\Users\apoliti\Desktop\mfluxtest\220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_04b\220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_04b.mat')
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

refBeads = load('C:\Users\apoliti\Desktop\mfluxtest\220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b\refBeads.mat')

end

% %
% figure(1)
% al = 0.15;
% clf
% s1 = scatter3(indata.(fn{1}).lnc(:, 1), indata.(fn{1}).lnc(:, 2), indata.(fn{1}).lnc(:, 3), 4, 'filled', 'MarkerFaceColor','blue');
% hold
% s2 = scatter3(indata.(fn{2}).lnc(:, 1), indata.(fn{2}).lnc(:, 2), indata.(fn{2}).lnc(:, 3), 4, 'filled', 'MarkerFaceColor','red');
% alpha(s1, al)
% alpha(s2, al)
% 
% 
% %%
% figure(2)
% clf
% plot3(indata.(fn{1}).loc(:, 1), indata.(fn{1}).loc(:, 2), indata.(fn{1}).loc(:, 3), '.', 'MarkerSize', 4, 'MarkerFaceColor', 'blue')
% hold
% plot3(indata.(fn{2}).loc(:, 1), indata.(fn{2}).loc(:, 2), indata.(fn{2}).loc(:, 3), '.', 'MarkerSize', 4, ...
%     'MarkerFaceColor', 'red')
% 
% daspect([1,1,1]);
% grid on;
% box on;
% view(-44,15);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title(['3D loc, default registration']);
% 
% 
% figure(3)
% clf
% plot3(indata.(fn{1}).ltr(:, 1), indata.(fn{1}).ltr(:, 2), indata.(fn{1}).ltr(:, 3), '.', 'MarkerSize', 4, 'MarkerFaceColor', 'blue')
% hold
% plot3(indata.(fn{2}).ltr(:, 1), indata.(fn{2}).ltr(:, 2), indata.(fn{2}).ltr(:, 3), '.', 'MarkerSize', 4, ...
%     'MarkerFaceColor', 'red')
% 
% daspect([1,1,1]);
% grid on;
% box on;
% view(-44,15);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title(['3D ltr, post aquisition translation']);
% 
% 
% 
% 
% 
% 
% %%
% refBeads = load('C:\Users\apoliti\Desktop\mfluxtest\220309_VGlut_paint_2nM_3DMINFLUX_16p_PH0_6_05b\refBeads.mat')
% figure(4)
% ax1 = subplot(1,3,1);
% plot3(refBeads.lnc(:, 1), refBeads.lnc(:, 2), refBeads.lnc(:, 3), '.', 'MarkerSize', 4);
% daspect([1,1,1]);
% grid on;
% box on;
% view(-44,15);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title(['3D lnc channel ' num2str(i)]);
% 
% ax2 =  subplot(1,3,2);
% plot3(refBeads.ltr(:, 1), refBeads.ltr(:, 2), refBeads.ltr(:, 3), '.', 'MarkerSize', 4);
% daspect([1,1,1]);
% grid on;
% box on;
% view(-44,15);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title(['3D lnc channel ' num2str(i)]);
% linkaxes([ax1, ax2])
% 
% ax2 =  subplot(1,3,3);
% plot3(refBeads.lre(:, 1), refBeads.lre(:, 2), refBeads.lre(:, 3), '.', 'MarkerSize', 4);
% daspect([1,1,1]);
% grid on;
% box on;
% view(-44,15);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title(['3D lnc channel ' num2str(i)]);
% linkaxes([ax1, ax2])