%% Saving figures in better resolution
% first save figures as matlab figures, then just switch out the file and save name

open hubbleprime.fig

ax = gca;                                   % keep axis' as manually set
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.ZTickMode = 'manual';
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
ax.ZLimMode = 'manual';

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 6];

print('hubbleprime','-djpeg','-r300')            % save as .jpeg file, 300 DPI (dots-per-inch) resolution

