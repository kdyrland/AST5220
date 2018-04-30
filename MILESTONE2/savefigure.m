%% Saving figures in better resolution
% first save figures as matlab figures and switch out file- and save name

open visibilityfct.fig

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

print('visibilityfct','-djpeg','-r300')
% save as .jpeg file, 300 DPI (dots-per-inch) resolution

