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
%      which is also released under the GNU General Public License.  
%      SPM must be installed and running for vis to work
%
%    * KDE estimation via:
%      http://www.mathworks.com/matlabcentral/fileexchange/14034-kernel-density-estimator
%  
%    * Color generation for KDE and histogram plots via:
%      http://www.mathworks.com/matlabcentral/fileexchange/29702
%
%    * High-quality figure exporting via:
%      http://www.mathworks.com/matlabcentral/fileexchange/23629
%
%    * The MATLAB "statistics" toolbox may be required for some functions.
%
%    * Documentation is available at http://tools.robjellis.net
%
%    * For more information, please contact Rob Ellis at robjellis@gmail.com
%
%
%    INSTALLATION:
%
%    1. Unzip vis.zip. In Linux, use this command: unzip <path>/vis.zip
%    2. Place the unzipped folder within the "toolbox" folder in your SPM directory.
%       This should result in vis being available upon rebooting of SPM.
%    3. To run the program, type "vis" at the MATLAB prompt. 
%
%
%    Copyright (C) 2011 by Robert J Ellis | http://robjellis.net
%   
    
%% path stuff
% must start here to get imcalc at the top of the path
spm_path_fix('vis.m'); % spm_path_fix is by RJE

% get save date of this file
loc = which('vis.m'); % note: must have a unique name and NOT an internal variable name

file_info = dir(loc); % note: there cannot be a variable "dir" active in Matlab or this won't work
save_date = file_info.date;
save_date = save_date(1:11); % trim off the timestamp for clarity


% fprintf('  _                    _       ')
% fprintf('\n (_)                  | |      ')
% fprintf('\n  _ _ __ ___   ___  __| | __ _ ')
% fprintf(['\n | | ''_ ` _ \\ / _ \\/ _` |/ _` |  Exploratory Data Analysis using ' spm('ver') ])
% fprintf(['\n | | | | | | |  __/ (_| | (_| |  This software version: ' save_date ]')
% fprintf('\n |_|_| |_| |_|\\___|\\__,_|\\__,_|  (C) Robert J Ellis (http://tools.robjellis.net) \n')

    % get spm version and revision
    [vv, rr] = spm('ver');
    
fprintf('        _       ')
fprintf('\n       (_)    ')
fprintf('\n __   ___ ___ ')
fprintf(['\n \\ \\ / / / __|   Visualized statistics and EDA toolbox (using ' vv ' rev. ' rr ')'])
fprintf(['\n  \\ V /| \\__ \\   This software version: ' save_date ])
fprintf('\n   \\_/ |_|___/   (C) Robert J Ellis (http://tools.robjellis.net)')
fprintf('\n                ')                        

%=============================================   
% main menu
%=============================================      

choice = input(['\n Choose a menu option below, or hit <ENTER> to quit: ',...
                '\n    [1] Probability density function (KDE)',... % opt 1
                '\n    [2] Histogram outline',... % opt 2
                '\n    [3] Surface plot (axial)',... % opt 3
                '\n    [4] 2D plots (Scatter, Bland-Altman, Q-Q plot)  ', ... % opt 4
                '\n    [5] Six-figure plot (X vs. Y)',... % opt 5
                '\n    [6] Extract parameter values from ROI(s)', ... % opt 6
                '\n    [7] Parcellate activation clusters', ... % opt 7
                '\n    [8] [change directory]',... % opt 8
                '\n     --> ']); 
          
% Old version
% vistitle = strcat('Visualization Tools (' ,save_date,')');
% 
% choice = menu_cent(vistitle,...
%               'Kernel density estimation of PDF',... % opt 1
%               'Histogram outline',... % opt 2
%               'Surface plot (axial)',... % opt 3
%               '  2D plots (Scatter, Bland-Altman, Q-Q plot)  ', ... % opt 4
%               'Six-plot (X vs. Y)',... % opt 5
%               'Extract parameter values', ... % opt 6
%               'Parcellate clusters', ... % opt 7
%               '[change directory]',... % opt 8
%               '[quit]'); % opt 9

if isempty(choice)
   % remove copies of SPM functions from the path to avoid errors
   remove_subfolder('vis.m','spm_fxns')
   
   fprintf('    Goodbye! \n\n');
   return
end


active_choices = 7; % easier to preserve last two options
 
if choice == 1 || choice == 2 || choice == 6 || choice == 7
    % histograms, extract, parcellate
    [fstats] = vis_getvol(choice);
    

elseif choice == 4 || choice == 5
    % 2D plots, Six-plot
    [plot_stats] = vis_corrs(choice);
     
elseif choice == 3
    [slice_output] = vis_slice;
   
elseif choice == active_choices + 1
    % change to a directory
    newdir = cfg_getfile(1,'dir','Select new directory ...');
    newdir = char(newdir);
    cd(newdir)
    
elseif choice < 1 || choice > 8
   fprintf(' Error: Input not recognized. Terminating. \n')
   return
end     

%% easy exporting of figures
% http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig?s_eid=PSM_3901
% https://sites.google.com/site/oliverwoodford/software/export_fig

if choice <= 5
   exp_fig = input(' Export current figure as .png?\n   [1] Yes \n   [2] No <ENTER> \n   --> ');
   
   if isempty(exp_fig)
       exp_fig = 2;
   end
   
   if exp_fig == 1
       % pdf requires ghostscript; see https://sites.google.com/site/oliverwoodford/software/export_fig
       % exp_type = input('    File type:  [1] .pdf; [2] .png: ');
       
      exp_type = 2;
      
      exp_name = input('     Desired file name: ','s');
      
      if isempty(exp_name)
          exp_name = 'vis_output';
      end
      
      if exp_type == 1
          exp_inst = ['export_fig ' exp_name '.pdf -q101']; % lossless compression
      elseif exp_type == 2
          exp_inst = ['export_fig ' exp_name '.png -m2.5'];
      end
      
      % do it
      fprintf('     Working ...')
      eval(exp_inst)
      
      fprintf([' Done. \n Figure "' exp_name '.png" saved to the working directory (''pwd''). \n ']);
   end
   
end

% after vis_getvol is finished, then we display this:
if choice == 1 || choice == 2
    fprintf('\n Variable "fstats" contains all related statistics.\n\n');
end

% clear unneeded variables
clear active_choices choice exp_fig file_info loc save_date

% remove copies of SPM functions from the path to avoid errors
remove_subfolder('vis.m','spm_fxns')
