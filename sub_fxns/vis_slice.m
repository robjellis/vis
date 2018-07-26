function [output] = vis_slice
%
% Description:
% To create a surface plot of an axial single brain slice (default is z = 0).
%
% last update = 2012.09.15

% Copyright (C) 2011 by Robert J Ellis   

% clear all variables
clear vol file org slice

fprintf('\n 3-D surface plot of voxel values in a single axial slice.\n');

V = spm_vol(spm8_select(1,'image','Select the image to analyze: ',[],pwd,'.*'));

zval = input('\n Enter z slice location in mm, <ENTER> for z = 0 (origin): ');

if isempty(zval)
    zval = 0;
end

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

fprintf(['\n\n ** Note the orientation in which x-axis and y-axis values increase in MATLAB: \n\n    x-axis: 0 to ' num2str(size_vol(1)) ' (Left to Right); \n    y-axis: 0 to ' num2str(size_vol(2)) ' (Posterior to Anterior); \n    z-axis: 0 to ' num2str(size_vol(3)) ' (Inferior to Superior).\n\n']);
