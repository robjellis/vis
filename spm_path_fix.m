function spm_path_fix(tar_prog)

% code by RJE to (1) make sure select matlab functions are in the path and
% (2) make sure the target directory is on top (so that any other uses of
% function names (e.g., "nanmean") will not be used incorrectly)
%
% tar_prog will be a character string of the target program (e.g., 'vis')
%
% rje | 2013.03.16, updated 2017 Nov


% look for spm.m
if exist('spm.m','file') ~= 2
    
   % get compatible version of Matlab
    matver = version('-release');

    set1 = {'13SP1','13SP2','14','14SP1','14SP2'}; % SPM5 only
    set2 = {'14SP3','2006a','2006b'}; % SPM5 or SPM8
    set3 = {'2007a','2007b','2008a','2008b','2009a','2009b','2010a','2010b','2011a','2011b','2012a','2012b','2013a','2013b','2014a','2014b','2015a'}; % SPM8 or SPM12
    set4 = {'2015b','2016a','2016b','2017a','2017b',}; % SPM12 only

    if max(strcmp(matver,set1)) == 1
        ok_ver = 'SPM5';
    elseif max(strcmp(matver,set2)) == 1
        ok_ver = 'Either SPM5 or SPM8';
    elseif max(strcmp(matver,set3)) == 1
        ok_ver = 'Either SPM8 or SPM12';
    elseif max(strcmp(matver,set4)) == 1
        ok_ver = 'SPM12';
    else
        ok_ver = ['An SPM version compatible with Matlab ' matver];
    end 
    
   fprintf(['\n Error: ' ok_ver ' must be on the Matlab path. ',...
            '\n Add the *main* SPM directory now. <' tar_prog '> will then relaunch.\n\n'])
   pause(0.5)
   
   spm_parent = uigetdir(pwd,'Select the main SPM directory');
   addpath(spm_parent); % need the main SPM folder added, but NOT all the subdirectories, toolboxes, etc.
   
   % now that SPM is on the path, we can determine its version
   spmver = spm('ver');
   spmver = spmver(4:end); % cut out 'SPM'
   
   if strcmp(spmver,'2') || strcmp(spmver,'5')
      batch_present = 0; 
   elseif strcmp(spmver,'8') || strcmp(spmver,'12b') || strcmp(spmver,'12')
      batch_present = 1; 
   end
   
   if batch_present == 1 % then we need matlabbatch
       spm_target = [spm_parent filesep 'matlabbatch']; % the target directory
       addpath_recurse(spm_target); % will add the subdirectories of this
   end
   
   % get the directory of tar_prog
   loc = which(tar_prog);
   prog_name = dir(loc);
   prog_name = prog_name.name;
   dir_name = loc(1:(end - numel(prog_name)));
   
   % now move to the top - this is important so that we always use RJE functions by default!
   path(dir_name,path)
   
   % also add subdirectories
   addpath_recurse(dir_name)
   
   return % don't need to relaunch the program
   
else % still need to move the target directory to the top of the path!
    
   % get the directory of tar_prog
   loc = which(tar_prog);
   prog_name = dir(loc);
   prog_name = prog_name.name;
   dir_name = loc(1:(end - numel(prog_name)));
   
   % now move to the top
   path(dir_name,path)
   
   % also add subdirectories
   addpath_recurse(dir_name)
   
   % don't need to relaunch the program
end