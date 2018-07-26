function flip = spm_flip_analyze_images

% This is a copy of spm_flip_analyze_images.m from SPM8/SPM12, for use in
% code by RJ Ellis (http://robjellis.net).
% 
% The value of "flip" is set to either 0 or 1. SPM has a default of "flip = 1",
% meaning that the indices of the voxels are stored as "left is right";
% flipping the images about the y-axis will then restore "left is left".
%
% Unless you are working with legacy image data, it is unlikely that you will
% ever need to edit this function. 
%
% Keeping "flip = 0" is necessary for working with images from the
% face_rfx dataset (http://www.fil.ion.ucl.ac.uk/spm/data/face_rfx/).

flip = 0; % flip = 0 means that Analyze images will be read in and written out without flipping.

