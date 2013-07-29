% Set Manuel's color preferences for programming with Matlab and Fortran

%% Setup the color pallet
                                   
% divide by 256 for Matlab         
sol.base03  = [ 59  59  59] / 256; 
sol.base02  = [  7  54  66] / 256; 
sol.base01  = [ 88 110 117] / 256; 
sol.base00  = [101 123 131] / 256; 
sol.base0   = [128 255 128] / 256; 
sol.base1   = [147 161 161] / 256; 
sol.base2   = [238 232 213] / 256; 
sol.base3   = [253 246 227] / 256; 
sol.yellow  = [181 137   0] / 256; 
sol.orange  = [255,128,  0] / 256; 
sol.red     = [220  50  47] / 256; 
sol.magenta = [211  54 130] / 256; 
sol.violet  = [108 113 196] / 256; 
sol.blue    = [  0 128 255] / 256; 
sol.cyan    = [ 42 161 152] / 256; 
sol.green   = [133 153   0] / 256;
sol.gray    = [160,160,160] / 256;

%% Setup Color Scheme for Light/Dark Options
    sycb = 0;           % Set 'Use system colors' checkbox to False
    txc  = sol.base0;   % Set 'Text' color
    bgc  = sol.base03;  % Set 'Background' color
    kwd  = sol.orange;  % Set 'Keywords' color
    cmt  = sol.gray;    % Set 'Comments' color - spec calls for base1, not sure about that...
    str  = sol.blue;    % Set 'Strings' color
    ustr = sol.orange;  % Set 'Unterminated strings' color
    scmd = sol.yellow;  % Set 'System commands' color
    errs = sol.red;     % Set 'Errors' color
    hyp  = sol.blue;    % Set 'Hyperlinks' color
    warn = sol.orange;  % Set 'Warnings' color
    afhb = 0;           % Set 'Autofix highlight' checkbox to False
    afh  = sol.base00;  % Set 'Autofix highlight' color
%     ahib = 0;           % Set 'Automatically highlight' checkbox to False
%     ahi  = sol.violet;  % Set 'Automatically highlight' color
    vwss = sol.green;   % Set 'Variables with shared scope' color
    clhb = 0;           % Set 'Highlight cells' checkbox to False
    hclb = 1;           % Set 'Highlight current line' checkbox to True
    hcl  = sol.base02;  % Set 'Highlight current line' color
    slnb = 1;           % Set 'Show line numbers' checkbox to True
    rtlb = 1;           % Set 'Show line' checkbox in Right-hand text limit to True

%% Reset colors
%     sycb = 1;           % Set 'Use system colors' checkbox to True
% Need to get colors setup to return to default

%% Desktop tool colors
com.mathworks.services.Prefs.setBooleanPref('ColorsUseSystem',sycb); clear('sycb')

com.mathworks.services.Prefs.setColorPref('ColorsText',java.awt.Color(txc(1), txc(2), txc(3))); 
com.mathworks.services.ColorPrefs.notifyColorListeners('ColorsText'); clear('txc')

com.mathworks.services.Prefs.setColorPref('ColorsBackground',java.awt.Color(bgc(1), bgc(2), bgc(3))); 
com.mathworks.services.ColorPrefs.notifyColorListeners('ColorsBackground'); clear('bgc')

%% Setup MATLAB syntax highlighting colors
com.mathworks.services.Prefs.setColorPref('Colors_M_Keywords',java.awt.Color(kwd(1), kwd(2), kwd(3)));
com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Keywords'); clear('kwd')

com.mathworks.services.Prefs.setColorPref('Colors_M_Comments',java.awt.Color(cmt(1), cmt(2), cmt(3)));
com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Comments'); clear('cmt')

com.mathworks.services.Prefs.setColorPref('Colors_M_Strings',java.awt.Color(str(1), str(2), str(3)));
com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Strings'); clear('str')

com.mathworks.services.Prefs.setColorPref('Colors_M_UnterminatedStrings',java.awt.Color(ustr(1), ustr(2), ustr(3)));
com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_UnterminatedStrings'); clear('ustr')

com.mathworks.services.Prefs.setColorPref('Colors_M_SystemCommands',java.awt.Color(scmd(1), scmd(2), scmd(3)));
com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_SystemCommands'); clear('scmd')

com.mathworks.services.Prefs.setColorPref('Colors_M_Errors',java.awt.Color(errs(1), errs(2), errs(3)));
com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Errors');clear('errs')

%% Other colors
com.mathworks.services.Prefs.setColorPref('Colors_HTML_HTMLLinks',java.awt.Color(hyp(1), hyp(2), hyp(3)));
com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_HTML_HTMLLinks'); clear('hyp')

%% Code analyzer colors
com.mathworks.services.Prefs.setColorPref('Colors_M_Warnings',java.awt.Color(warn(1), warn(2), warn(3)));
com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Warnings'); clear('warn')

com.mathworks.services.Prefs.setBooleanPref('ColorsUseMLintAutoFixBackground',afhb);
com.mathworks.services.Prefs.setColorPref('ColorsMLintAutoFixBackground',java.awt.Color(afh(1), afh(2), afh(3)));
com.mathworks.services.ColorPrefs.notifyColorListeners('ColorsMLintAutoFixBackground'); clear('afhb','afh')

%% Variable and function colors
% com.mathworks.services.Prefs.setBooleanPref('Editor.VariableHighlighting.Automatic',ahib);
% com.mathworks.services.Prefs.setColorPref('Editor.VariableHighlighting.Color',java.awt.Color(ahi(1), ahi(2), ahi(3)));
% com.mathworks.services.ColorPrefs.notifyColorListeners('Editor.VariableHighlighting.Color'); clear('ahib','ahi')

% *** need to add checkbox option for variables with shared scope here
com.mathworks.services.Prefs.setColorPref('Editor.NonlocalVariableHighlighting.TextColor',java.awt.Color(vwss(1), vwss(2), vwss(3)));
com.mathworks.services.ColorPrefs.notifyColorListeners('Editor.NonlocalVariableHighlighting.TextColor');clear('vwss')

%% Cell display options
com.mathworks.services.Prefs.setBooleanPref('EditorCodepadHighVisible',clhb); clear('clhb')
% Highlight cells color is 'Editorhighlight-lines'

%% Editor/Debugger General display options
com.mathworks.services.Prefs.setBooleanPref('Editorhighlight-caret-row-boolean',hclb);
com.mathworks.services.Prefs.setColorPref('Editorhighlight-caret-row-boolean-color',java.awt.Color(hcl(1), hcl(2), hcl(3)));
com.mathworks.services.ColorPrefs.notifyColorListeners('Editorhighlight-caret-row-boolean-color'); clear('hclb','hcl')
com.mathworks.services.Prefs.setBooleanPref('EditorShowLineNumbers',slnb); clear('slnb')
com.mathworks.services.Prefs.setBooleanPref('EditorRightTextLineVisible',rtlb);
%EditorRightTextLineLimit=I80
%EditorRightTextLimitLineWidth=I1

%% Cleanup
clear('sol'); clear('rtlb')