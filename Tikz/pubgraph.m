%publish graph function (makes graphs pretty and publishable)
function pubgraph(fh,FS,LW,CL)
%fh is the handles to the figure containing the graph
%FS is the Font Size 
%LW is the Line Width
%CL is Color

%% ------------------------------------------------------------------------
% cool ideas: http://www.mikesoltys.com/2011/09/16/tool-of-the-week-prettyer-matlab-plots/
%--------------------------------------------------------------------------
%%
figure(fh)                                         	%pull the figure forwards
axs = findall(fh, 'Type', 'axes');                 	%get the axes on the figure
set(axs, 'FontSize', FS, 'LineWidth', LW,'Color',CL)%make everything on the axis correct
set(findall(fh, 'Type', 'text'), 'FontSize', FS); 	%make all other text correct ,'FontName','Calibri'
set(findall(fh, 'Type', 'line'), 'LineWidth',LW);	%make all other lines correct
set(fh,'Color',CL)                                  %set the figure background color

%for each graph on the figure
for i=1:length(axs)
%get the curent axis
    ax=axs(i); 
%     axis(ax,'tight')
%make sure y axis isn't too tight
    v = ylim;                                     	%get current y axis limits
    sc = 0.02*diff(v);                           	%add in 2% each side
    ylim([v(1)-sc , v(2)+sc])                     	%apply it to the graphs
    % v=xlim; sc=0.02*diff(v); xlim([v(1)-sc v(2)+sc])	%make sure x axis isn't too tight
    
%format X and Y axis tickmarks so all numbers are displayed with the same precision
    YTick = str2num(get(ax,'YTickLabel'));          %get the current Y axis labels
    XTick = str2num(get(ax,'XTickLabel'));          %get the current X axis labels
    %check what precision the y axis is
    if ~any(rem(YTick,1)),
        yfmt='%-5.0f'; 
    elseif ~any(rem(YTick*10,1)),
        yfmt='%-5.1f'; 
    elseif ~any(rem(YTick*100,1)),
        yfmt='%-5.2f';
    elseif ~any(rem(YTick*1000,1)),
        yfmt='%-5.2f';
    else
        yfmt='%-5.4f';
    end
    %check what precision the x axis is
    if ~any(rem(XTick,1)),
        xfmt='%-5.0f'; 
    elseif ~any(rem(XTick*10,1)),
        xfmt='%-5.1f'; 
    elseif ~any(rem(XTick*100,1)),
        xfmt='%-5.2f';
    elseif ~any(rem(XTick*1000,1)),
        xfmt='%-5.3f';
    else
        xfmt='%-5.4f';
    end
    %set the axis precision
    set(ax,'YTick',get(ax,'YTick'),'YTickLabel',num2str(YTick,yfmt));%set all the y labels to the same precision
    set(ax,'XTick',get(ax,'XTick'),'XTickLabel',num2str(XTick,xfmt));%set all the x labels to the same precision
end