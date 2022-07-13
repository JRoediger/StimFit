function [faxchild,txt,fig_cont,fig_main] = openvizzfig(monrev,optdist,options,sidename)


figure('Name',options.figname,'NumberTitle','off','Units','normalized','Position',[0.1148,0.1847,0.8020,0.6097])
load contfig_segmented

xvals = [0,options.ampl(end)];
yvals = [1,options.nconts];

% Child 8
faxchild{10} = subplot(2,9,18);
imshow(fig_main)

% Child 7
faxchild{9} = subplot(2,9,13);
imshow(fig_main)

% Child 6
faxchild{6} = subplot(2,9,15:17);
xlabel('iteration')
ylabel('optimal outcome')
xticks([])
yticks([])
title(['Multipolar Settings (', sidename{2},')'])
faxchild{8} = animatedline;

% Child 5
faxchild{5} = subplot(2,9,10:12);
xlabel('iteration')
ylabel('optimal outcome')
xticks([])
yticks([])
title(['Multipolar Settings (', sidename{1},')'])
faxchild{7} = animatedline;


% Child 4
faxchild{4} = subplot(2,9,9);
imshow(fig_main)

% Child 3
faxchild{3} = subplot(2,9,4);
imshow(fig_main)

% Child 2
faxchild{2} = subplot(2,9,6:8);
xlabel('amplitude [mA]')
ylabel('contact')
xlim(xvals)
ylim(yvals)
imagesc(nan(options.nconts,length(options.ampl)))
set(gca,'YDir','normal')
title(['Monopolar Settings (', sidename{2},')'])

% Child 1
faxchild{1} = subplot(2,9,1:3);
xlabel('amplitude [mA]')
ylabel('contact')
xlim(xvals)
ylim(yvals)
imagesc(nan(options.nconts,length(options.ampl)))
set(gca,'YDir','normal')
title(['Monopolar Settings (', sidename{1},')'])

colormap(flipud(pink))
txt = annotation('textbox', [0.13, 0.05, 0.8, 0], 'string', 'Loading figure...','FontSize',10,'EdgeColor','none');
fax = gcf;

set(fax, 'MenuBar', 'none', 'ToolBar', 'none');
