% Setup
DB = Meta_Setup(DB, 10);
MC_Setup = Meta_Activation_FWE('setup', DB);

% Enter:
% 1
% SelfOthe
% [1 4]
% [1 -1]
% Other-Self

%% Plot
MC_Setup.connames{1} = 'SelfOthe';
Meta_plot_points_on_slices(DB, MC_Setup)
enlarge_axes(gcf, 1.4)
enlarge_axes(gcf, 1.2)
scn_export_papersetup(700); saveas(gcf, 'otherblue_selfred_points_slices_1.fig');
enlarge_axes(gcf, 1.2)
ph = findobj(gcf, 'Type', 'line');
set(ph, 'MarkerFaceColor', 'none');

%% Detail of one slice, up close

%f2 = create_figure;  % just do this once

% click on a slice, then:
ax1 = gca;
figure(f2); clf
ax2 = copyobj(ax1, f2);
axes(ax2)
colormap gray
camzoom(8)
ph = findobj(ax2, 'Type', 'line', 'color', 'b')
set(ph, 'MarkerFaceColor', 'b');
ph = findobj(ax2, 'Type', 'line', 'color', 'r')
set(ph, 'MarkerFaceColor', 'r');
ph = findobj(ax2, 'Type', 'line', 'color', 'g')
set(ph, 'MarkerFaceColor', 'g');
axis off


