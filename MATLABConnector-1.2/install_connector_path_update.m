function install_connector_path_update(shouldSave)

productName = 'MATLAB Connector';

% add directories from addon .phl file to pathdef.m, current path
msg = sprintf('Adding %s folders to the MATLAB Path...', productName);
disp(msg)

% stash current path, pathdef before re-creating pathdef
current_path = path;
saved_path = pathdef;

% turn off duplicate path warning before modifying path
w_state = warning('off', 'MATLAB:dispatcher:pathWarning');

% recreate pathdef to get newly added .phl file into pathdef
% This requires that a phl file be created/be present for the connector
restoredefaultpath;
path(saved_path, path);
if(shouldSave)
    if (savepath ~= 0)
        disp('Warning: Unable to save modified path to file.')
        msg = sprintf('To have %s on the path for future MATLAB sessions, you will need to save the path to a different file.', productName);
        disp(msg)
    end
end

% rebuild current path with newly added paths
path(current_path, path);

% restore duplicate path warning
warning(w_state);
end