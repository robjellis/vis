function [plot_stats] = vis_corrs(choice)
%
% Image comparison utility: x-y plot of values in two images on a voxel-by-voxel basis.
%
% Optimized for SPM5/8/12 and requires it to be running in the background.
%
% Assumes both images have the SAME DIMENSIONS.
%
% last update = 2017.11.06

% Copyright (C) 2011 by Robert J Ellis
    

%fprintf('\n\n Pearson correlation between two (sets of) images at every voxel. \n 1. Select [Group1 and Group2] files. \n 2. Select whether to exclude zeros from 1st, 2nd, or neither image. \n 3. Result will plot in a new figure. \n\n');

file_sel = input('\n File selection \n   [1] Single pair of files from one directory \n   [2] Separate selection of Group1 and Group2 files \n   [3] Single group, homotopic correlation \n   --> ');

if file_sel == 1
    files = spm_select(2,'image','Select a pair of files:',[],pwd,'.*');
    files1 = files(1,:);
    files2 = files(2,:);
    numf = 1;
elseif file_sel == 2
    files1=spm_select([1 Inf],'image','Select x-axis file(s):',[],pwd,'.*');
    files2=spm_select([1 Inf],'image','Select y-axis file(s):',[],pwd,'.*');
    
    if size(files1,1) == size(files2,1)
    % OK
        numf = size(files1,1);
    else
        fprintf('\n Warning: x-axis files and y-axis files do not have the same number. Terminating.\n\n')
        return
    end

elseif file_sel == 3
    files1=spm_select([1 Inf],'image','Select (x-axis) file(s):',[],pwd,'.*');
    numf = size(files1,1);
end

rnd = 0;

% select the masking template

use_mask = input(' Use implicit mask? \n   [1] Yes \n   [2] No <ENTER> \n   --> ');
if isempty(use_mask)
   use_mask = 2;
end

if use_mask == 1
    if file_sel == 1 || file_sel == 2
        mask=spm_select(1,'image','Specify inclusive mask:',[],pwd,'.*');
    elseif file_sel == 3
        mask=spm_select(1,'image','Specify a *single hemisphere* inclusive mask:',[],pwd,'.*');
    end
 vm = spm_read_vols(spm_vol(mask));
 vm(isnan(vm)) = 0;  % replaces NaNs with zero
 vm = vm ~= 0; % binarize non-zero values
 [x1 y1 z1] = size(vm);
 
  vm2 = vm(:);
 
end 

minex = input(' Only include values above V? \n   [1] Yes \n   [2] No <ENTER> \n   --> ');
if isempty(minex)
   minex = 2;
end

if minex == 2
    
    excz = input(' Exclude voxels with zeros in \n   [1] First image \n   [2] Second image \n   [3] Either image <ENTER> \n   --> ');
    if isempty(excz)
       excz = 3;
    end

elseif minex == 1
    minval = input(' Enter minimum value V: ');
    excz = 3;
end

% useful variables


% plot figs; optional; perhaps we just want to calculate the stats

if numf > 1
    plotfigs = input('Show plots? \n   [1] Yes <ENTER>\n   [2] No \n   --> ');
else
    plotfigs = 1;
end

if isempty(plotfigs) 
   plotfigs = 1;
end

% record important stats

r_val = zeros(numf,1);
p_val = zeros(numf,1);
z_val = zeros(numf,1);
nc_val = zeros(numf,1); % normalized correlation of Jenkinson et al. (2002)
i1_mn = zeros(numf,1);
i2_mn = zeros(numf,1);
i1_sd = zeros(numf,1);
i2_sd = zeros(numf,1);
BA_mean = zeros(numf,1); % for B-A plot
BA_error = zeros(numf,1); % for B-A plot

% *********************
% file loop

fprintf('\n Working ... ');

% make separate plots for PAIRS of files
if numf == 1
   fignum = 1001;
else
   fignum = 1001:(1001+numf-1); 
end

for i = 1:numf % works on successive PAIRS of files

% read in both volumes
 
 namef1 = strtrim(files1(i,:)); % strtrim gets rid of white space characters
 %namef1b = namef1;
  if namef1(numel(namef1) - 1) == ','   % gets rid of ",1" at end of file name (as needed)
     namef1 = namef1(1:(numel(namef1)-2)); 
  end
     namef1b_file = dir(namef1);
     [namef1b_file] = namef1b_file.('name');
     
 vol1 = spm_read_vols(spm_vol(namef1));

if file_sel == 1 || file_sel == 2  
     namef2 = strtrim(files2(i,:));
     %namef2b = namef2;
     if namef2(numel(namef2) - 1) == ','
         namef2 = namef2(1:(numel(namef2)-2));
     end
         namef2b_file = dir(namef2);
         [namef2b_file] = namef2b_file.('name');
         
     vol2 = spm_read_vols(spm_vol(namef2));
     
elseif file_sel == 3
    % just flip the file1 volume
    vol2 = flipdim(vol1,1); % flip in x-dimension
    namef2b_file = strcat(namef1b_file,'_x-flip');
    
end

 
 % compare sizes of volumes
 
 if sum(size(vol1) - size(vol2)) == 0
     %fprintf('\n \n The scatter plot images have the same dimensions ...\n');
 else
     fprintf('\n\n\n Error: Volumes do NOT have the same dimensions. Terminating. \n\n\n');
     return
 end
 
 % make sure same size with mask
if use_mask == 1

    if sum(size(vm) - size(vol1)) == 0
        %fprintf('\n Mask and images have the same dimensions ... \n');
    else
        fprintf('\n\n\n Error: Mask and images do NOT have the same dimensions. Terminating. \n\n\n');
        return
    end
end



% turn each volume into a single column of data, which will be the x and y
% coordinates of the plot we want to see

v1s = vol1(:);
v2s = vol2(:);


% turn NaNs into zeros

v1s(isnan(v1s)) = 0;
v2s(isnan(v2s)) = 0;

%size(v1s);
%size(v2s);



% retain file values where the mask is not 0 
if use_mask == 1
v1s = v1s(vm2 ~= 0);  
v2s = v2s(vm2 ~= 0);  

end

ss = v1s + v2s;

% retain values from the x and y columns where their sum does not equal 0;
% gets rid of non-brain zero values

v1b = v1s(ss ~=0);
v2b = v2s(ss ~=0);


% exclude zeros option


%if excz == '3'
%    v1 = v1b;           % retain all values
%    v2 = v2b;
    
if excz == 1
    % exclude zeros based on 1st image mask
    v1 = v1b(v1b ~= 0);
    v2 = v2b(v1b ~= 0);

elseif excz == 2
    % exclude zeros based on 2nd image mask
    v1 = v1b(v2b ~= 0);
    v2 = v2b(v2b ~= 0);
    
elseif excz == 3
    % only include non-zero voxels in BOTH images
    v1t = v1b(v1b ~= 0);
    v2t = v2b(v1b ~= 0);
    v1 = v1t(v2t ~= 0);
    v2 = v2t(v2t ~= 0);
    
end

if minex == 1
    
    v1t = v1b(v1b > minval);           % retain values above minval
    v2t = v2b(v1b > minval);

    v1 = v1t(v2t > minval);
    v2 = v2t(v2t > minval);
end

im1vals = v1;
im2vals = v2;

v3 = [v1, v2];

if isempty(v3)
    fprintf('\n\n Error: There are no voxels that meet this criteria when zeros are excluded. Respecify. \n\n');
    vis_corrs
end

min1 = floor(min(v1));
max1 = ceil(max(v1));
min2 = floor(min(v2));
max2 = ceil(max(v2));

min3 = min(min1,min2);
max3 = max(max1,max2);

bb = [min1,max1,min2,max2];
min4 = min(bb);
max4 = max(bb);

%the pearson correlation between v1 and v2; outputs r and p values
% use "corrcoef" rather than "corr" because the former does not require the
% stats toolbox
[rval pval] = corrcoef(v1,v2);  
[rval pval] = corrcoef(1,2);

% the normalized correlation [cost function] of Jenkinson et al. (2002)
nc = sum(v1.*v2) / (sqrt(sum(v1.^2))*sqrt(sum(v2.^2)));


r_val(i) = rval;
p_val(i) = pval;
z_val(i) = atanh(rval);  % Fisher's z transform: z = 0.5*ln((1+r)/(1-r)) = atanh(r)
nc_val(i) = nc;
i1_mn(i) = mean(v1);
i2_mn(i) = mean(v2);
i1_sd(i) = std(v1);
i2_sd(i) = std(v2);


% for Bland Altman
        
% transformation for Bland-Altman plot
% see http://en.wikipedia.org/wiki/Bland%E2%80%93Altman_plot

ba1 = (im1vals + im2vals) / 2;
ba2 =  im1vals - im2vals;

% temporary
ba_mean = mean(ba2);
ba_err = 1.96 * std(ba2);

% save these
BA_mean(i) = ba_mean;
BA_error(i) = ba_err;


% =========================
% figures; don't have to plot, just optional

    if plotfigs == 1

        if choice == 4
           % just the 2D plots
           %figname = '2D plots';
           plot_eda(v1,v2,[4 5 6],fignum(i),namef1b_file,namef2b_file);
        elseif choice == 5
           % six-plot
           %figname = 'Six-plot';
           plot_eda(v1,v2,[1 2 3 4 5 6],fignum(i),namef1b_file,namef2b_file);
        end


    end

end % file loop

stats.r_vals = r_val;
stats.p_vals = p_val;
stats.z_vals = z_val;
stats.nc_vals = nc_val;
stats.i1_mn = i1_mn;
stats.i1_sd = i1_sd;
stats.i2_mn = i2_mn;
stats.i2_sd = i2_sd;
stats.i1_vals = v1; % just the most recent v1
stats.i2_vals = v2; % just the most recent v2

if choice == 5
    stats.BA_mean = BA_mean;
    stats.BA_error = BA_error;
end

%=== print results
plot_stats = stats;

fprintf('\n Analysis complete. \n Table of statistics results saved as "plot_stats".\n\n');


