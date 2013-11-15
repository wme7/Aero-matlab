function fh = setaxes(varargin)
%SETAXES   creates handles of nested functions for fixing problems with MATLAB figures.
%  
%   SETAXES with no output argument make automatic adjustments to the  
%   current axes position propery if parts of the axes is clipped from 
%   sides due to small figure size.
%
%   SETAXES(AX) takes the axes with handle AX for the current axes.
%
%   This first step works for common cases. For complete control over
%   axes position and other features, use the advanced syntax below.
%
%   FH = SETAXES(..) makes no automatic adjustments and returns a structure  
%   FH of 9 function handles for 2-D view and 3 function handles for 3-D view.
%   A brief description of the purpose of each function is:
%
%     Function                    Purpose of the nested function handle
%   -----------------    ---------------------------------------------------
%        'FH.xoffset'        To move and stretch the axes in horizontal direction
%        'FH.yoffset'        To move and stretch the axes in vertical direction
%   'FH.legendmoveresize'    To shrink and relocate the legend
%
%
%   Function                 Purpose of the nested function handle
%     (2D-view only)
%   -----------------   ---------------------------------------------------
%     'FH.xtick2text'     To replace XTickLabels by text objects
%     'FH.ytick2text'     To replace YTickLabels by text objects
%   'FH.xlabelcorner'     To place XLabel at the appropriate corner 
%   'FH.ylabelcorner'     To place YLabel at the appropriate corner   
%     'FH.axesarrows'     To add arrows to the axes lines
%     'FH.d1annotation'   To facilitate 1D annotation objects in data space units 
%     'FH.d2annotation'   To facilitate 2D annotation objects in data space units 
%     
%    
%
%
%    FH.xoffset() adjusts the axes position in the horizontal direction to
%    make up for any clipped (left or right) part of the axes based on some
%    default settings.
%
%
%
%    FH.yoffset() adjusts the axes position in the vertical direction to
%    make up for any clipped (lower or upper) part of the axes based on 
%    some default settings.
%
%    FH.axesmoveresize([L B R T]) adds the array [L B -R -T] to the   
%    'OuterPosition' property  of the axes. 
%
%
%    FH.legendmoveresize(SHX,SHY) shrinks the legend width by a scalar factor 
%    SHX and its height by a scalar factor SHY. SHX and SHY must be positive
%    reals not exceeding 1.
%    The fontsize and linewidth of the legend and its children are resized 
%    by a factor MIN(SHX,SHY) for consistant look.
%
%       FH.legendmoveresize(SHX) is the same as FH.legendshrink(SHX,SHX) and    
%    FH.legendshrink() is the same as FH.legendshrink(.75)
%    
%    FH.legendmoveresize(SHX,SHX,X,Y) additionally relocates the lower-left corner
%    of the legend at (X,Y) in the data space.
%    
%
%
%    FH.xtick2text(OFFSET) replaces the tick labels on the X-Axis by text
%    objects where OFFSET is the factor of the axis height used as the   
%    vertical margin from the X-Axis. 
%
%       FH.xtick2text() is the same as FH.xtick2text(0.03). 
%
%    FH.xtick2text(..., PROP1, VALUE1, PROP2, VALUE2, ...)
%       sets the values of the specified properties of the text objects.
%    
%    HT = FH.xtick2text(...) returns in HT the handles of the text objects.
% 
% 
%    FH.ytick2text(OFFSET) replaces the tick labels on the Y-Axis by text
%    objects where OFFSET is the factor of the axis width used as the   
%    horizontal margin from the Y-Axis. 
%
%       FH.ytick2text() is the same as FH.ytick2text(0.02). 
%
%    FH.ytick2text(..., PROP1, VALUE1, PROP2, VALUE2, ...)
%       sets the values of the specified properties of the text objects.
%    
%    HT = FH.ytick2text(...) returns in HT the handles of the text objects. 
% 
% 
%    FH.xlabelcorner() places the XLABEL in an appropriate corner. The most   
%    common corner is the right-bottom when the X-Axis is to the bottom 
%    and the ticks are increasing from left to right. The corner changes 
%    accordingly for other combinations of axes properties.
% 
%    FH.xlabelcorner(PROP1, VALUE1, PROP2, VALUE2, ...)
%    sets the values of the specified properties of the XLabel object.
%
%
%    FH.ylabelcorner() places the YLABEL in an appropriate corner. The most   
%    common corner is the top-left when the Y-Axis is to the left 
%    and the ticks are increasing upward. The corner changes accordingly  
%    for other combinations of axes properties. 
% 
%    FH.ylabelcorner(PROP1, VALUE1, PROP2, VALUE2, ...)
%    sets the values of the specified properties of the YLabel object.
% 
% 
%    FH.axesarrows() adds arrow annotation objects to X- and Y-Axis. The default
%    'HeadLength' and 'HeadWidth' are both 5. Line width and Color are 
%    inherited from the the corresponfing axis line.
%
%    FH.axesarrows('xx') adds arrow annotation objects to X-Axis only. 
%    FH.axesarrows('yy') adds arrow annotation objects to Y-Axis only. 
%    
%    FH.axesarrows(..., PROP1, VALUE1, PROP2, VALUE2, ...)
%    sets the values of the specified properties of the arrow objects.
%
%    HA = FH.axesarrows(...) returns in HA the handles of the arrow objects.
%
%    Note that since axesarrows() changes the axes position property, any
%    annotation objects should be added after axesarrows().
%
%
%    [XFIG, YFIG] = FH.d1annotation(X,Y) transforms (axes) data coordinates (X,Y) 
%    into normalized (figure) coordinates (XFIG, YFIG) for placing 
%    1D annotation objects (arrow, doublearrow, textarrow) into data space.
%
%    PFIG = FH.d2annotation(P) transforms the (axes) data position 
%    P =[X, Y, W, H] into normalized (figure) position 
%    P =[XFIG, YFIG, WFIG, HFIG] for placing 2D annotation objects 
%    (rectangle, textbox, ellipse) into data space 
%
%=======================================================================
%                   EXAMPLE trying to cover a few features:
%-----------------------------------------------------------------------
% 
% figure('Units', 'centimeters', ...
%        'Position', [2 2 8 6], ...
%        'PaperUnits', 'centimeters', ...       
%        'PaperSize', [8 6], ...
%        'PaperPositionMode', 'auto')
% 
% ax = axes('FontSize', 8, ...
%            'LineWidth', 0.4, ...
%            'Box', 'off', ...
%            'TickDir', 'out' );
% %    
% hold all
% 
% t = linspace(0,5,51);
% n = 1:49;
% rand('state', sum(100*clock));
% W = [0 cumsum(-log(rand(size(n)))./n)];
% hPlot(1) = stairs(W, [n 50], '-b');
% hPlot(2) = line(t, exp(t), 'LineStyle', '--', 'Color', 'r');
% axis([0 5 0 50])
% 
% xlabel('$T=\mathrm{Exp}(\frac{1}{X})$', 'Interpreter','LaTeX');
% ylabel('$X$', 'Interpreter','LaTeX');
% legend(hPlot, {'sample $n$'; 'mean $\mu(T)$'}, 'Interpreter','LaTeX')
% set(ax, 'XTickLabel', {'$0$';'$\tau_1$';'$\tau_2$'; ...
%                          '$\tau_3$';'$\tau_4$';'$\tau_5$'}, ...
%         'YTickLabel', {'$0$';'$X_1$';'$X_2$';'$X_3$';'$X_4$';'$X_5$'} )
% 
% fh = setaxes;
% fh.xtick2text()
% fh.ytick2text()
% fh.ylabelcorner()
% fh.yoffset();
% fh.axesarrows()
% fh.legendmoveresize(.8,.8,.2,20)
% [xfig,yfig] = fh.d1annotation([1 2], [45 45]);
% har = annotation('doublearrow', xfig , yfig);
% set(har, 'LineWidth', .25, 'Color', 'k', ...
%          'Head1Length', 3.5, 'Head2Length', 3.5, ...
%          'Head1Width', 3.5, 'Head2Width', 3.5, ...
%          'Head1Style', 'vback3', 'Head2Style', 'vback3')
%      
% text(1.5, 48, '$\Delta \tau$', ...
%             'VerticalAlignment', 'bottom', ...
%             'HorizontalAlignment', 'center', 'Interpreter','LaTeX');
% 
% print -depsc setaxes_example.eps
% print -dpdf setaxes_example.pdf
%
%
%=========================================================================

% Mukhtar Ullah
% mukhtar.ullah@informatik.uni-rostock.de
% Sep 18, 2008

%=========================================================================
[ax,args,nargs] = extractoptional([], @isnumeric, varargin{:});

if isempty(ax)
    ax = gca;
elseif ~isscalar(ax) || ~ishandle(ax) || ~strcmp(get(ax,'type'),'axes')
    error('Invalid Axes handle')
end

if nargs>0 && rem(nargs,2)~=0
    error('Incorrect number of input arguments')
else
    args = {args{:}, 'Units', 'normalized'};
end

set(ax, args{:});

axprop = get(ax, {'Position', ...
                  'OuterPosition', ...
                  'LooseInset', ...
                  'XTick', ...
                  'YTick', ...                  
                  'XLim', ...
                  'YLim', ...
                  'Xlabel', ...
                  'Ylabel', ...
                  'Title'} );
              
numvar = numel(axprop);

axprop  = [axprop, get(ax, {'XAxisLocation', ...
                            'YAxisLocation', ...
                            'XDir', ...
                            'YDir', ...
                            'XScale', ...
                            'YScale', ...
                            'XTickLabelMode', ...
                            'YTickLabelMode'} )];
               
[axpos, outpos, looseinset, xtick, ytick, axlim(1:2), axlim(3:4), ...  
                hxlab, hylab, htitle, xxloc, yyloc] = axprop{1:numvar+2};

proptf = num2cell(strcmp(axprop(numvar+1:end), ...
   {'bottom', 'left', 'normal', 'normal', 'log', 'log', 'auto', 'auto'}));
                  
[isxxbottom, isyyleft, isxdirnormal, isydirnormal, ...
          isxxlog, isyylog, autoxticklab, autoyticklab] = proptf{:};

[dxt, dyt, dxtfig, dytfig, xtickroom, ytickroom] = deal(0);                                        
isxlabelcorner = false;
isylabelcorner = false;
                                        
notxxloc = setdiff({'top', 'bottom'}, xxloc);
notyyloc = setdiff({'left', 'right'}, yyloc);
set([hxlab, hylab, htitle], {'Tag'}, {'XLabel'; 'YLabel'; 'Title'})
xlabpos = get(hxlab, 'Position');
ylabpos = get(hylab, 'Position');

[axwidth, axheight, ...
          xlabhalign, ylabhalign, xlabvalign, ylabvalign] = falign();
                
[az,el]=view(ax);
isax2D = (az==0) && (el==90);
%--------------------------------------------------------------------------

if nargout==0
    yoffset()
    xoffset()
else
    names = {'xoffset', 'yoffset', 'axesmoveresize', 'legendmoveresize'};
    
    if isax2D
        names = { names{:}, 'xtick2text', 'ytick2text', ...
                            'xlabelcorner', 'ylabelcorner', ...
                            'axesarrows', 'd1annotation', 'd2annotation'};
    end
    
    handles = cellfun(@(s) eval(['@', s]), names, 'UniformOutput', false); 
    fh = cell2struct(handles, names, 2);
end

%==========================================================================
    function xoffset()
        if ~isylabelcorner && get(hylab, 'Rotation')==0
            set(hylab, 'HorizontalAlignment', ylabhalign)
        end
        tightinset = get(ax, 'TightInset');
        set(ax, 'LooseInset', max(looseinset, tightinset))
        if isax2D
            restoreleftright(hylab)
        end
    end
%--------------------------------------------------------------------------
    function yoffset()
        tightinset = get(ax, 'TightInset');
        set(ax, 'LooseInset', max(looseinset, tightinset))
        if isax2D
            restoreupdown(hxlab)
            restoreupdown(htitle)
        end
    end
%--------------------------------------------------------------------------
    function axesmoveresize(offset)
        if nargin == 0            
            disp('Drag a rectangle in the axes to define OuterPosition:')
            k = waitforbuttonpress;
            if k == 1
                return
            end            
            fig = get(ax, 'Parent');
            oldunits = get(fig, 'Units');
            set(fig, 'Units', 'normalized')
            outpos = rbbox;
            set(fig, 'Units', oldunits)
        else            
            outpos = outpos + [1 1 -1 -1].*offset;
        end
        set(ax, 'OuterPosition', outpos)
        drawnow
        axpos = get(ax, 'Position');
    end
    function legendmoveresize(shx,shy,locx,locy)        
        
        [hleg,hlegobj] = legend(ax);
        
        if isempty(hleg)
            return
        end
        
        nargs = nargin;
        if nargs == 0
            set(hleg, 'Units', 'points','ButtonDownFcn','selectmoveresize')
            disp('Move/Resize Legend and then press any key to continue:')
            pause
            return
        end
        
        if isempty(shx)
            shx = .75;
        end

        if nargs < 2 || isempty(shy)
            shy = shx;
        end

        if any([shx shy] == 0) || any(abs([shx shy]) > 1)
            error('SHX and SHY must be positive reals not exceeding 1')
        end

        shobj = min(shx,shy);

        set(hleg, 'FontSize', shobj*get(hleg, 'FontSize'), ...
                  'LineWidth', shobj*get(hleg, 'LineWidth')    );

        hlegline = findobj(hlegobj, 'Type', 'line');

        leglinelinewidth = ...
            cellfun(@(x) shobj*x, get(hlegline, 'LineWidth'));

        set(hlegline, {'LineWidth'}, num2cell(leglinelinewidth))

        oldunits = get(ax, 'Units');
        set([ax,hleg], 'Units', 'points')
        axpospt = get(ax, 'Position');
        set(ax, 'Units', oldunits)
        legpos = get(hleg, 'Position');

        legrightx = legpos(1) + legpos(3);
        legrighty = legpos(2) + legpos(4);

        axrightx = axpospt(1) + axpospt(3);
        axrighty = axpospt(2) + axpospt(4);

        shx = shx*legpos(3);
        shy = shy*legpos(4);

        isOutsideLeft = legrightx < axpospt(1);
        isOutsideRight = legpos(1) > axrightx;
        isInsideRight = ~isOutsideLeft && ~isOutsideRight && ...
            (legpos(1)-axpospt(1) > axrightx-legrightx);

        isOutsideBottom = legrighty < axpospt(2);
        isOutsideTop = legpos(2) > axrighty;
        isInsideTop = ~isOutsideBottom && ~isOutsideTop && ...
            (legpos(2)-axpospt(2) > axrighty-legrighty);

        if isInsideRight && isOutsideLeft
            legpos(1) = legpos(1) + legpos(3) - shx;
        end


        if isInsideTop && isOutsideBottom
            legpos(2) = legpos(2) + legpos(4) - shy;
        end

        legpos(3:4) = [shx shy];

        if exist('locx', 'var') && ~isempty(locx)
            [locxfig, locyfig] = dsxy2figxy(locx, axlim(3));
            legpos(1) = locxfig*axpospt(3)/axpos(3);
        end

        if exist('locy', 'var') && ~isempty(locy)
            [locxfig, locyfig] = dsxy2figxy(axlim(1),locy);
            legpos(2) = locyfig*axpospt(4)/axpos(4);
        end

        set(hleg, 'Position', legpos)
        
    end  
 %--------------------------------------------------------------------------   
    function hOut = xtick2text(varargin)
        [xtickoffset, pvpairs] = extractoffset(.03, varargin{:});
        hxt = [];
        xticklab = cellstr(get(ax, 'XTickLabel'));
        xticklab = matchticks(xtick, xticklab);
        if ~isempty(xticklab)  
            set(ax, 'XTickLabel', []);
            
            if xor(isxxbottom,isydirnormal)
                yloc = axlim(4) + xtickoffset*axheight;
            else
                yloc = axlim(3) - xtickoffset*axheight;
            end
            
            yloc = yloc(ones(size(xtick)));

            if isyylog
                yloc = 10.^yloc;
                xticklab = append10(xticklab, autoxticklab, pvpairs);
            end

            hxt = text(xtick, yloc, xticklab);
            
            set(hxt(strmatch('$', xticklab)), 'Interpreter', 'latex')
            
            set(hxt, 'HorizontalAlignment', 'center', ...
                     'VerticalALignment', xlabvalign, pvpairs{:});
            %--------------------------------------------------------------
            oldunits = get(hxt(1), 'Units');
            set(hxt, 'Units', 'normalized')
            xtickext = get(hxt, 'Extent');
            xtickroom = max(cellfun(@(x) x(4), xtickext));
            set(hxt, 'Units', oldunits)
            
            xtickext = get(hxt, 'Extent');
            xlabposnow = get(hxlab, 'Position');
            
            if isxxbottom
                xlabposnow(2) = ...
                    min(cellfun(@(x) x(2), xtickext)) + xtickoffset;
            else
                xlabposnow(2) = ...
                    max(cellfun(@(x) x(2)+x(4), xtickext)) + xtickoffset;
            end
            
            set(hxlab, 'Position', xlabposnow)

        end
        if nargout
            hOut = hxt;
        end        
    end
%--------------------------------------------------------------------------
    function hOut = ytick2text(varargin)
        [ytickoffset, pvpairs] = extractoffset(.02, varargin{:});           
        hyt = [];
        yticklab = cellstr(get(ax, 'YTickLabel'));
        yticklab = matchticks(ytick, yticklab);
        if ~isempty(yticklab)
            set(ax, 'YTickLabel', []);
            
            if xor(isyyleft,isxdirnormal)
                xloc = axlim(2) + ytickoffset*axwidth;
            else
                xloc = axlim(1) - ytickoffset*axwidth;
            end
            
            xloc = xloc(ones(size(ytick)));
            
            if isxxlog
                xloc = 10.^xloc;
                yticklab = append10(yticklab, autoyticklab, pvpairs);
            end

            hyt = text(xloc, ytick, yticklab);
            
            set(hyt(strmatch('$', yticklab)), 'Interpreter', 'latex')            

            set(hyt, 'VerticalAlignment', 'middle', ...
                     'HorizontalALignment', ylabhalign, pvpairs{:})
            %--------------------------------------------------------------
            oldunits = get(hyt(1), 'Units');
            set(hyt, 'Units', 'normalized')
            ytickext = get(hyt, 'Extent');
            ytickroom = max(cellfun(@(x) x(3), ytickext));
            set(hyt, 'Units', oldunits)
            
            ytickext = get(hyt, 'Extent');
            ylabposnow = get(hylab, 'Position');
            
            if isyyleft
                ylabposnow(1) = min(cellfun(@(x) x(1), ytickext));
            else
                ylabposnow(1) = max(cellfun(@(x) x(1)+x(3), ytickext));
            end
            
            set(hylab, 'Position', ylabposnow)

        end
        if nargout
            hOut = hyt;
        end        
    end
%--------------------------------------------------------------------------
    function xlabelcorner(varargin)        
        if isxxlog
            xlabpos(1:2) = 10.^xlabpos(1:2);
        end        

        set(hxlab, ...
            'Position', xlabpos, ...
            'VerticalALignment', 'middle', ...
            'HorizontalALignment', xlabhalign, varargin{:});

        isxlabelcorner = true;       
        axpos = get(ax, 'Position');
        restoreleftright(hxlab)        
    end
%--------------------------------------------------------------------------
    function ylabelcorner(varargin)
        if isyylog
            ylabpos(1:2) = 10.^ylabpos(1:2);
        end

        set(hylab, ...
            'Position', ylabpos, ...
            'Rotation', 0, ...
            'VerticalALignment', ylabvalign, ...
            'HorizontalALignment', 'center', varargin{:});
        
        isylabelcorner = true;        
        axpos = get(ax, 'Position');
        restoreupdown(hylab)        
    end
%--------------------------------------------------------------------------
    function hOut = axesarrows(varargin)
        args = varargin;
        xxarrow = true;
        yyarrow = true;
        
        if nargin>0
            switch args{1}
                case 'xx'
                    yyarrow = false;
                    args(1) = [];
                case 'yy'
                    xxarrow = false;
                    args(1) = [];
                otherwise
                    error('Unknown option')
            end
        end
        
       [left, bottom, right, top] = deal(0); 
                
        if xxarrow && ~isxlabelcorner
            if isyyleft
                right = dxtfig;
            else
                left = dxtfig;
            end
        end

        if yyarrow && ~isylabelcorner
            if isxxbottom
                top = dytfig;
            else
                bottom = dytfig;
            end
        end
        
        axesmoveresize([left, bottom, right, top])
        
        if isyyleft
            xxarx = axpos(1)+axpos(3)+[0 dxtfig];
            yyarx = axpos([1 1]);
        else
            xxarx = axpos(1)+[axpos(3) -dxtfig];
            yyarx = axpos([1 1]) + axpos([3 3]);
        end

        if isxxbottom
            yyary = axpos(2)+axpos(4)+[0 dytfig];
            xxary = axpos([2 2]);
        else
            yyary = axpos(2)+[axpos(4) -dytfig];
            xxary = axpos([2 2]) + axpos([4 4]);
        end
        
        if xxarrow
            h(1) = annotation('arrow', xxarx, xxary, ...
                                      'Color', get(ax, 'XColor'));
        end

        if yyarrow
            h(2) = annotation('arrow', yyarx, yyary, ...
                                      'Color', get(ax, 'YColor'));
        end
        
        set(nonzeros(h), 'HeadLength', 5, 'HeadWidth', 5, ...
                         'LineWidth', get(ax, 'LineWidth'), args{:})

        if nargout==1
            hOut = h;
        end
    end
%--------------------------------------------------------------------------
    function [xfig,yfig] = d1annotation(x,y)
        if nargin == 0
            disp('Click two points on the graph to place 1D annotation objects:')
            [x, y] = ginput(2);           
        end
        [xfig, yfig] = dsxy2figxy(x, y);
    end
%--------------------------------------------------------------------------
    function [bbox] = d2annotation(abox)
        if nargin == 0
            disp(['Drag a rectangle in the axes ', ...
                  'to define the bounding box of 2D annotation objects:'])
            k = waitforbuttonpress;
            point1 = get(ax,'CurrentPoint');    % button down detected
            rbbox;
            point2 = get(ax,'CurrentPoint');    % button up detected
            point1 = point1(1,1:2);              % extract x and y
            point2 = point2(1,1:2);
            abox(1:2) = min(point1,point2);
            abox(3:4) = abs(point1-point2);
        end
        [bbox] = dsxy2figxy(abox);
    end
%--------------------------------------------------------------------------
    function varargout = dsxy2figxy(varargin)
        error(nargchk(1,2,nargin,'struct'));
        args = varargin;

        if isscalar(args)          % Must be a 4-element POS vector
            pos = args{1};
            if ~any(size(pos)==4)
                error('Position vector must have 4 rows or 4 columns')
            end
        else
            [x,y] = deal(args{:});  % Two tuples (start & end points)
        end

        %--------------------------------------
        if exist('x','var')

            if isxxlog, x = log10(x); end

            if isyylog, y = log10(y); end

            if isxdirnormal
                dx = x - axlim(1);
            else
                dx = axlim(2) - x;
            end

            if isydirnormal
                dy = y - axlim(3);
            else
                dy = axlim(4) - y;
            end

            varargout{1} = dx*axpos(3)/axwidth + axpos(1);
            varargout{2} = dy*axpos(4)/axheight + axpos(2);

        else

            if isxxlog, pos([1 3]) = log10(pos([1 3]));  end

            if isyylog, pos([2 4]) = log10(pos([2 4]));  end

            if isxdirnormal
                dxpos = pos(1) - axlim(1);
            else
                dxpos = axlim(2) - pos(1);
            end

            if isydirnormal
                dypos = pos(2) - axlim(3);
            else
                dypos = axlim(4) - pos(2);
            end

            pos(1) = dxpos/axwidth*axpos(3) + axpos(1);
            pos(2) = dypos/axheight*axpos(4) + axpos(2);
            pos(3) = pos(3)*axpos(3)/axwidth;
            pos(4) = pos(4)*axpos(4)/axheight;
            varargout{1} = pos;

        end

    end
%==========================================================================
    function restoreupdown(lab)
        withlab = ~isempty(get(lab, 'String')) && ...
                            strcmp(get(lab, 'Visible'), 'on');
        if withlab
            oldunits = get(lab, 'Units');
            set(lab, 'Units', 'normalized')
            labext = get(lab, 'Extent');

            if strcmpi(get(lab, 'Tag'), 'YLabel')
                labext = labext.*axpos([3 4 4 3]);
            else
                labext = labext.*axpos([3 4 3 4]);
            end

            extroom = labext(4);

            if strcmpi(get(lab, 'Tag'), 'XLabel')
                set(lab, 'Visible', 'off')
                keyword = xxloc;
                extroom = 1.2*max(extroom, xtickroom*axpos(4));
            else
                keyword = notxxloc;
            end

            tightinset = get(ax, 'TightInset');
            set(lab, 'Visible', 'on', 'Units', oldunits)
            if strcmp(keyword, 'bottom')
                looseinset(2) = max(looseinset(2), tightinset(2) + extroom);
            else
                %                 looseinset(4) = max(looseinset(4), labext(2)+labext(4)-axpos(4));
                looseinset(4) = max(looseinset(4), tightinset(4) + extroom);
            end

            set(ax, 'LooseInset', looseinset)
            axpos = get(ax, 'Position');
        end
    end
%--------------------------------------------------------------------------
    function restoreleftright(lab)
        withlab = ~isempty(get(lab, 'String')) && ...
                                strcmp(get(lab, 'Visible'), 'on');
        if withlab
            oldunits = get(lab, 'Units');
            set(lab, 'Units', 'normalized')
            labext = get(lab, 'Extent');

            if strcmpi(get(lab, 'Tag'), 'YLabel')
                set(lab, 'Visible', 'off')
                labext = labext.*axpos([3 4 4 3]);
                if get(lab, 'Rotation')==90
                    extroom = labext(4);
                else
                    extroom = labext(3);
                end
                extroom = max(extroom, ytickroom*axpos(3));
                keyword = yyloc;
            else
                labext = labext.*axpos([3 4 3 4]);
                extroom = labext(3);
                keyword = notyyloc;
            end
            
            extroom = 1.2*extroom;

            tightinset = get(ax, 'TightInset');
            set(lab, 'Visible', 'on', 'Units', oldunits)

            if strcmp(keyword, 'left')
                looseinset(1) = max(looseinset(1), tightinset(1) + extroom);
            else
                %                 looseinset(3) = max(looseinset(3), labext(1)+labext(3)-axpos(3));
                looseinset(3) = max(looseinset(3), tightinset(3) + extroom);
            end

            set(ax, 'LooseInset', looseinset)
            axpos = get(ax, 'Position');
        end
    end
%--------------------------------------------------------------------------
    function [ axwidth, axheight, xlabhalign, ylabhalign, ...
                                  xlabvalign, ylabvalign ] = falign()

        if isxxlog
            axlim(1:2) = log10(axlim(1:2));
            xlabpos(1) = log10(xlabpos(1));
        end

        if isyylog
            axlim(3:4) = log10(axlim(3:4));
            ylabpos(2) = log10(ylabpos(2));
        end

        axwidth = diff(axlim(1:2));
        axheight = diff(axlim(3:4));

        dxtfig = .02;
        dxt = max(0.05*axwidth, dxtfig*axwidth/axpos(3));
        dxtfig = dxt*axpos(3)/axwidth;

        dytfig = .02;
        dyt = max(0.05*axheight, dytfig*axheight/axpos(4));
        dytfig = dyt*axpos(4)/axheight;

        if isyyleft
            xlabhalign = 'left';
            ylabhalign = 'right';
            if isxdirnormal
                xlabpos(1) = axlim(2) + 1.5*dxt;
                ylabpos(1) = axlim(1);
            else
                xlabpos(1) = axlim(1) + 1.5*dxt;
                ylabpos(1) = axlim(2);
            end
        else
            xlabhalign = 'right';
            ylabhalign = 'left';
            if isxdirnormal
                xlabpos(1) = axlim(1) - 1.5*dxt;
                ylabpos(1) = axlim(2);
            else
                xlabpos(1) = axlim(2) - 1.5*dxt;
                ylabpos(1) = axlim(1);
            end
        end

        if isxxbottom
            xlabvalign = 'top';
            ylabvalign = 'bottom';
            if isydirnormal
                ylabpos(2) = axlim(4)  + 1.5*dyt;
                xlabpos(2) = axlim(3);
            else
                ylabpos(2) = axlim(3)  + 1.5*dyt;
                xlabpos(2) = axlim(4);
            end
        else
            xlabvalign = 'bottom';
            ylabvalign = 'top';
            if isydirnormal
                ylabpos(2) = axlim(3)  - 1.5*dyt;
                xlabpos(2) = axlim(4);
            else
                ylabpos(2) = axlim(4)  - 1.5*dyt;
                xlabpos(2) = axlim(3);
            end
        end

    end
%--------------------------------------------------------------------------
    function [tickoffset, pvpairs] = extractoffset(tickoffset, varargin)
        pvpairs = varargin;
        npv = nargin - 1;
        if npv>0 && isnumeric(pvpairs{1})
            tickoffset = pvpairs{1};
            pvpairs(1) = [];
            npv = npv - 1;
        end

        if npv>0 && (rem(npv,2) ~= 0)
            error('Incorrect number of input arguments')
        end

        textpropnames = {'FontAngle', 'FontName', ...
            'FontUnits', 'FontSize', 'FontWeight'};

        textpv = [textpropnames; get(ax, textpropnames)];
        pvpairs = [textpv(:)' pvpairs];
    end
%--------------------------------------------------------------------------
    function ticklab = matchticks(tick, ticklab)
        ntick = numel(tick);
        nticklab = numel(ticklab);
        q = floor(ntick/nticklab);
        r = rem(ntick,nticklab);
        ticklab = ticklab([repmat(1:nticklab, 1, q) 1:r]);
    end
%--------------------------------------------------------------------------
    function ticklab = append10(ticklab, autoticklab, pvpairs)
        if autoticklab
            interpreter = get(0, 'DefaultTextInterpreter');
            if ~isempty(pvpairs)
                tparams = lower(pvpairs(1:2:end-1));
                tvals = pvpairs(2:2:end);
                k = strmatch('Interpreter', tparams, 'exact');
                if ~isempty(k)
                    interpreter = tvals{k+1};
                end
            end
            if strcmp(interpreter, 'tex')
                f10 = @(str) ['10^{', str, '}'];
            else
                f10 = @(str) ['$10^{', str, '}$'];
            end
            ticklab = cellfun(f10, ticklab, 'UniformOutput', false);
        end
    end

%--------------------------------------------------------------------------
    function [y, args, nargs]= extractoptional(y, validator, varargin)

        args = varargin;
        nargs = length(args);

        if nargs > 0 && validator(args{1})
            y = args{1};
            args = args(2:end);
            nargs = nargs-1;
        end
    end

end
