function remove_subfolder(funcname,subfolder_name)

% remove subfolders of the folder containing funcname.m from the Matlab path, but keep the folder
% containing funcname.m on the path.


loc = which(funcname);
prog_name = dir(loc);
prog_name = prog_name.name;
dir_name = loc(1:(end - numel(prog_name)));
check_this = [dir_name subfolder_name];

% is it on the path?
pathCell = regexp(path, pathsep, 'split');
onPath = any(strcmpi(check_this, pathCell));
if onPath == 1
   rmpath(check_this)
else
   % not on the path, so don't have to do anything
end