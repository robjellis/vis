
%
%    This program is a set of image visualization tools for MR images, and is 
%    made available the neuroimaging community as copyright freeware. 
%
%    You can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.  
%    See the GNU General Public License for more details.
%
%    A copy of the GNU General Public License is included in the main directory for vis: vis/gpl3.txt.
%          
%    NOTES:
%
%    * Portions of the code call SPM functions (http://www.fil.ion.ucl.ac.uk/spm), 
%      which is also released under the GNU General Public License.  SPM
%      must be installed and running for vis to work
%
%    * The MATLAB "statistics" toolbox must be installed for vis to function.
%
%    * Documentation is available at http://tools.robjellis.net
%
%    * For more information, please contact Rob Ellis at
%      robjellis@gmail.com
%
%
%    INSTALLATION:
%
%    1. Unzip vis.zip. In Linux, use this command: unzip <path>/vis.zip
%    2. Place the unzipped folder within the "toolbox" folder in your SPM directory.
%       This should result in vis being available upon rebooting of SPM.
%    3. To run the program, type "vis" at the MATLAB prompt. SPM must be active in the path,
%       but does not need to be actively running.
%
%
%    Copyright (C) 2011 by Robert J Ellis
%   
    
version = '2012.11.16';

fprintf(['\n Visualized Statistics toolbox for SPM software \n Version: ' version '\n Copyright (C) 2011 Robert J. Ellis \n\n']);



% get version
  ver = spm('ver');
  
  if str2num(ver(4)) == 2
     ver = 'SPM2';
  elseif str2num(ver(4)) == 5
     ver = 'SPM5'; 
  elseif str2num(ver(4)) == 8
     ver = 'SPM8';
  end
  

%=============================================   
% conc main menu
%=============================================      

vistitle = strcat('Visualization tools (v. ',version,')');

choice = menu(vistitle,...
              'Histogram outline(s)',...
              'Surface plot (axial)',...
              'Scatter plot', ...
              'Q-Q plot',...
              'Bland-Altman plot',...
              '  Extract parameter values  ', ...
              ' Parcellate clusters ', ...
              '[change directory]',...
              '[clear and exit]');

 active_choices = 7; % easier to preserve last two options
 
if choice == 1 || choice == 6 || choice == 7
    %if choice == 6
        % only works for SPM8 right now
       %if strcmp(ver,'SPM8')
           % OK
       %else
       %    fprintf(' ** Note: Only works for SPM8 at this time. ** \n\n');
       %    return
       %end
    %end
    
    [fstats, Xvals, Yvals] = vis_getvol(choice);
    if choice == 1
    fprintf('Finished. \n\n Variables may now be called: \n    "fstats": min, mean, max, std for each file; \n    "Xvals": unique histogram x-values for each file; \n    "Yvals": unique histogram y-values for each file.\n\n');
    elseif choice == 6
    fprintf('Finished. \n\n Variable "fstats" may now be called (min, mean, max, std for each file) \n\n');
    end
elseif choice == 3 || choice == 4 || choice == 5
    
    if choice == 4
          % requires stats toolbox
           if exist('quantile') == 2
               % OK
           else
               fprintf(' ** The MATLAB stats toolbox does not appear to be installed.\n    Q-Q plot is unavailable without the stats toolbox.\n\n')
               return
           end 
        
    end
    [corr_stats] = vis_corrs(choice);
     
elseif choice == 2
    [slice_output] = vis_slice;
   
elseif choice == active_choices + 1
    % change to a directory
    newdir = cfg_getfile(1,'dir','Select new directory ...');
    newdir = char(newdir);
    cd(newdir)
    
elseif choice == active_choices + 2  % clear all variables and return to the prompt
    
    clear version vistitle ver choice corr_stats
    return
end     

return
