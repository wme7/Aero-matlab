function install_connector()
%INSTALL_CONNECTOR  Install the MATLAB Connector.
%   INSTALL_CONNECTOR Install the files for the MATLAB Connector into your MATLAB installation.

%   Copyright 2010-2011 The MathWorks, Inc.

% turn off connector if running
try
    com.mathworks.util.ClassLoaderBridge.findClass('com.mathworks.matlabserver.embeddedwebserver.JettyWebServer');
    running = com.mathworks.matlabserver.embeddedwebserver.JettyWebServer.isRunning();
    if running
        connector('off');
        clear('java');
    end
catch e %#ok<NASGU>
end

% Determine location of this file
currentFile = mfilename('fullpath');
[currentDir, currentScriptName] = fileparts(currentFile);

% Define variables used in the script
productName = 'MATLAB Connector';
licenseAgreementFile = fullfile(currentDir, 'MATLAB_Connector_SLA.txt');
toolboxStr = 'toolbox';
connectorStr = 'connector';
connectorRoot = fullfile(matlabroot, toolboxStr, connectorStr);
connectorHelpFile = fullfile(connectorRoot, connectorStr, 'html', 'bsjg133.html');
helpFileStr = 'http://www.mathworks.com/help/matlabmobile';

% make sure user is using at least MATLAB R2012a
if verLessThan('matlab', '${matlab.minimum.version}')
   error('%s does not run on MATLAB versions earlier than R2012a.', productName);
end

% check if the jvm is available
if (~usejava('jvm'))
    error('%s requires Java to run.', currentScriptName)
end

% Display License Agreement
if (exist(licenseAgreementFile, 'file') ~= 2)
    error('The %s license agreement file does not exist on your machine.\nThis is an invalid installation.', productName)
end

type(licenseAgreementFile)
reply = input('Do you accept this license agreement? Y/N [N]: ', 's');
msg = sprintf('You must accept the license agreement if you wish to continue installing. \n...Exiting the %s Installation.', productName);
if(isempty(regexpi(reply,'^(y|yes)$')))
    disp(msg)
    return
end

doinstall
install_connector_path_update(false);
rehash toolboxcache

% Wrap-up
msg = sprintf('\nInstallation of %s is complete.\n', productName);
disp(msg)

msg = sprintf('To start the connector type \"connector on\".');
disp(msg)

msg = sprintf('\nFor more information see the <a href="matlab:web(''%s'')">Setup and User''s Guide</a>.', helpFileStr);
disp(msg)
end

function doinstall
productName = 'MATLAB Connector';
% check if the user can write to matlabroot
tempFileInMLRoot = tempname(matlabroot);
fid = fopen(tempFileInMLRoot, 'w');
if (fid == -1)
    failed = true;
    
    if ispc
        try
            proc = System.Diagnostics.Process;
            % EXE to run
            proc.StartInfo.FileName = fullfile(matlabroot,'bin','matlab');
            % Arguments to the EXE
            proc.StartInfo.Arguments = '/wait /noslpash /nodesktop /r install_connector_helper,exit';
            % Run-as admin
            proc.StartInfo.Verb = 'runas';
            proc.StartInfo.WindowStyle = ...
                System.Diagnostics.ProcessWindowStyle.Hidden;
            msg = sprintf(['\n------------------------------------------------------------------------\n' ...
                'The MATLAB Connector Installer is launching a separate MATLAB session,\n'...
                'and will close that session when the installation is complete.\n'...
                '------------------------------------------------------------------------']);
            disp(msg);
            proc.Start(); % Start
            proc.WaitForExit(); % Wait for the process to end
            failed = proc.ExitCode ~= 0;
        catch %#ok<CTCH>
        end
    end
    if (failed)
        error('\nError: You do not have write permission to the MATLAB Installation folder (%s).\nThe %s Installer cannot proceed.', matlabroot, productName);
    end
    return
else
    fclose(fid);
    delete(tempFileInMLRoot);
    install_connector_helper;
end
end