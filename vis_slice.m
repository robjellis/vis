function [output] = vis_slice
%
% Description:
% To create a surface plot of an axial single brain slice (default is z = 0).
%
% last update = 2012.09.15


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

% clear all variables
clear vol file org slice

fprintf('\n 3-D surface plot of voxel values in a single axial slice.\n');

V = spm_vol(spm_select(1,'image','Select the image to analyze: ',[],pwd,'.*'));

zval = input('\n Enter z slice location in mm, or type 0 (zero) for z-origin: ');


aftrans = V.mat; % the affine transformation matrix
xflip = sign(aftrans(1,1));    % -1 means flipped along x-axis
 
invV = inv(V.mat);
origin = invV(1:3,4); % the origin

vol =  spm_read_vols(V);  % the 3D volume
size_vol = size(vol);        % [x y z] dimensions

% if qfactor = -1, need to flip along the x-axis
if xflip == -1
   vol = flipdim(vol,1); 
end


if zval == 0
   z = vol(:,:,origin(3));
   
else
   zval2 = floor(zval/aftrans(3,3));
   z = vol(:,:,zval2+origin(3));
end

zstr = num2str(zval);

% now plot the surface

figure(20)
surf(z)
xlabel('y-axis'); ylabel('x-axis'); zlabel('Voxel intensity'); 
title(strcat('Surface plot of voxel values at z = ',zstr,' mm'));
colorbar;

% output useful information
output.size = size_vol;
output.origin = origin';
output.xflip = xflip;

fprintf(['\n\n ** Note the orientation in which x-axis and y-axis values increase in MATLAB: \n\n    x-axis: 0 to ' num2str(size_vol(1)) ' (left to right); \n    y-axis: 0 to ' num2str(size_vol(2)) ' (posterior to anterior); \n    z-axis: 0 to ' num2str(size_vol(3)) ' (inferior to superior).\n\n']);
