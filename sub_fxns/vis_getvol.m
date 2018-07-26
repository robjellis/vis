function [fstats, Xvals, Yvals] = vis_getvol(choice)

%
% Description:
% A simple program to read in all the (masked) values in a particular image;
% useful for getting percentiles or other stats.
%
% RJ Ellis | version = 2013.11.07

% Copyright (C) 2011 by Robert J Ellis

% clear variables
clear dim1 vmat v1b x1 y1 z1;

% get current directory
curdir = pwd;

% *************************
% 2013.11.07 - remap this option
if choice == 7
    choice = 8;
end
% *************************

if choice == 1 || choice == 2 || choice == 6
    % select the image to mask
    files = spm_select([1 Inf],'image','Select image(s):',[],pwd,'.*');
    numf = size(files,1);
    % select the masking template

    if choice == 1|| choice == 2
        use_mask = input(' Use inclusive mask / ROI mask?: [1] Yes; [2] No <ENTER>: ');
         if isempty(use_mask)
            use_mask = 2; % note: vis will still exclude voxels with 0 or NaN as their value
         end

    elseif choice == 6
        % ROI masking is assumed
        use_mask = 1;

    end
    
elseif choice == 8
    im_type = input(' Type of image: \n   [1] Statistic image (can apply a threshold); \n   [2] Cluster image of activations (integer labelled): ');
    use_mask = 1;
end


if use_mask == 1
    ver = spm('ver');
    if strcmp(ver,'SPM8') || strcmp(ver,'SPM12')
        % OK
    else
       fprintf('\n ** Ensure that masks have same dimensions as target image.\n');
    end
    
    if choice == 1 || choice == 2
        masks=spm_select(1,'image','Select a single inclusive mask / ROI mask:',[],pwd,'.*'); % only 1 mask at this point (RJE can change later)
    
    elseif choice == 6
        masks=spm_select([1 Inf],'image','Select ROI(s) images:',[],pwd,'.*'); % can have multiple ROIs (columns)
    
    elseif choice == 8 % we actually treat the clusters of interest as the mask(s)
        if im_type == 1
            masks=spm_select(1,'image','Select a *single* statistic image:',[],pwd,'.*'); 
        elseif im_type == 2
            masks=spm_select(1,'image','Select a *single* cluster image:',[],pwd,'.*'); 
        end    
            % need to see if it is binarized
            mvol = spm_read_vols(spm_vol(masks(1,:)));
            mvol(isnan(mvol)) = 0;
            
            % how many unique values in this image?
            mvals = numel(unique(mvol));
            
            if mvals < 200
                % assume cluster image (including a simple binarized image)
                im_type = 2;
                numm = numel(unique(mvol)) - 1; % just in case there are missing integers; also, we ignore 0
                mvals = unique(mvol);
                mvals = mvals(mvals~=0); % ignore zero
                
            else
                % if there are many unique values, assume it contains an
                % actual brain activaton map
                im_type = 1;
                im_thr = input(' Enter minimum height threshold (T- or Z-value) for statistic image: ');
                % will threshold later
                   
                numm = 1; % treat all suprathreshold voxels as belonging to a SINGLE cluster for now
                mvals = 1; % to make the code work
            end
            
    end
    
    if choice == 1 || choice == 2 || choice == 6
        numm = size(masks,1); % number of image masks
        nregs = 1; % to make code work
    elseif choice == 8
        % number of masks defined above
    end
    maskdir = fileparts(masks(1,:));
elseif use_mask == 2
    numm = 1; % we will just mask the image with itself
    maskdir = fileparts(files(1,:));
    
end

% now get the atlas
if choice == 1 || choice == 2 || choice == 6
    nregs = 1; % to make code work
elseif choice == 8
    files = spm_select(1,'image','Select the desired atlas:',[],pwd,'.*');
    numf = size(files,1);
    % how many regions?
    avol = spm_read_vols(spm_vol(files(1,:)));
    avol(isnan(avol)) = 0;
    avals = unique(avol); % number of unique atlas values (include zero, will see why)
    nregs = numel(avals);
    
end
if choice == 1
    % we don't have to worry about bins; kde.m uses a default Range and
    % number of points in the mesh (2e12)
    
    calc_perc = 0; % just to make code work
    nperc = 1; % to make code work
    minV = 2; % to make code work
    
elseif choice == 2
    commonX = input(' Set x-axis values automatically (for each image)?: [1] Yes <ENTER>; [2] No: ');
    
    if isempty(commonX)
        commonX = 1;
    end

    if commonX == 2
       xvals = input(' Enter the parameters [min:step:max] for x-vals, in [ ]: ');
       nbins = numel(xvals);

    elseif commonX == 1
            choosex = input(' Use 50 bins per image? [1] Yes <ENTER>; [2] No: ');   
            if isempty(choosex)
               choosex = 1;
            end

            if choosex == 1
               nbins = 50;
            elseif choosex == 2
               nbins = input(' How many bins?: ');
            end

    end

    calc_perc = 0; % just to make code work
    nperc = 1; % to make code work
    minV = 2; % to make code work
    
elseif choice == 6
    nbins = 1;
    calc_perc = 1;
    percit = input('\n Evaluate percentiles? \n   [1] Yes \n   [2] No <ENTER> \n   --> ');
    
    if isempty(percit)
       percit = 2;
    end
    
    if percit == 1
        percit2 = input(' Method \n   [1] Fixed for all images \n   [2] Unique for each file \n   --> ');
        if percit2 == 1
            perc = input(' Enter a *single* percentile value (e.g., 50): ');
            perc = zeros(numf,1) + perc;
            nperc = 1;
        elseif percit2 == 2
            perc = input([' Paste in a [' num2str(numf) ' x 1] or [' num2str(numf) ' x 1] vector of percentile values: ']);
            % get into a column
            perc = perc(:);
            nperc = 1; % for simplicity, we assume just a single percentile per file
            
            if numf ~= size(perc,1) % should be equal to the number of files
                % bad
                fprintf(' Warning: the number of files and percentiles does not match.\n\n');
                return
            end
        end
        
        %if exist('prctile')
           %fprintf(' (Stats toolbox is installed; will use MATLAB prctile function.)\n')
           %prc_ch = 1;
        %else
           %fprintf(' (Stats toolbox is not installed; will use SPM prctile function.)\n') 
           prc_ch = 2; % use this by default now
        %end
        
    elseif percit == 2
        calc_perc = 0; % just to make code work
        nperc = 1; % just to make code work
    end
    
   
   minV = input('\n Exclude values < V \n   [1] Yes \n   [2] No <ENTER> \n   --> '); 
   
   if isempty(minV)
       minV = 2;
   end
   
   if minV == 1
      minT = input('\n   Enter the threshold V: ');
   end
   
elseif choice == 8 % just to make code work
   calc_perc = 0; % to make code work
   nperc = 1; % to make code work
   minV = 2; % to make code work
   nbins = 1;
end

% the list below is for choice == [2 6 7]
percs_data = nan(numf,nperc,numm);
percs_kde = nan(numf,nperc,numm);
numvox = zeros(numf,numm);
mins = zeros(numf,numm);
means = zeros(numf,numm);
means_pos = zeros(numf,numm);
maxs = zeros(numf,numm);
stds = zeros(numf,numm);
skews = zeros(numf,numm);
sums = zeros(numf,numm);
fzeros = zeros(numf,numm);
dims = zeros(numf,3);
origin = zeros(numf,3);

if choice == 1
   kde_mesh = 2^12;
   Xvals = zeros(kde_mesh,numf);  % will be unique set of xvals and yvals per image
   Yvals = zeros(kde_mesh,numf);
elseif choice == 2
   Xvals = zeros(nbins,numf);  % will be unique set of xvals and yvals per image
   Yvals = zeros(nbins,numf);
end



% the following matrix is just for choice == 8
cls_vox = zeros(nregs+1,numm+1); % will count the number of voxels for each atlas region (rows), separate for each cluster (column)
cls_vox(1,1) = NaN; % first column will be atlas index; first row witll be cluster image value

fprintf('\n Working ... ');

% a few more things for cluster parcellation first ...

% first, let's make sure that the images and masks have the same
% dimensions. We just check the first mask, and hope for the best.

if use_mask == 1
    chkmask = strtrim(masks(1,:));
    cmv = spm_read_vols(spm_vol(chkmask));

    chkfile = strtrim(files(1,:));
    cfv = spm_read_vols(spm_vol(chkfile));

    if sum(size(cmv) - size(cfv)) == 0
      % we are ok
      use_res = 0;
    else
      % need to reslice
      use_res = 1;

      fprintf(' \n Note: Mask(s) have a different dimension than the target image. Reslicing ...\n');

      % masks will have the prefix "rs_"
      % assume nearest neighbor for simplicity!
      % some code borrowed from Tor Wager's reslice_imgs.m

        flags = struct('interp', 0, ... % rje: neearest neighbor
            'mask', 0, ...              % do not mask
            'mean', 0, ...              % do not write mean image
            'hold', -1, ...             % TW: i don't think this is used anymore
            'which', 1, ...             % reslice 2nd-nth only
            'wrap', [1 1 0]', ...       % the default (rje: SPM says [1 1 0] for fMRI)
            'prefix','rs_' ...           
            );

        P = str2mat(chkfile,masks);
        spm_reslice(P, flags) % copy of spm_reslice for SPM12

    end

end

% ===================================
if choice == 1 || choice == 2 || choice == 6
    for m = 1:numm
        % -------------------------------
        % mask loop

        if use_mask == 1

            if use_res == 1
                % use the resliced mask
                mask = strtrim(masks(m,:)); % cut out text whitespace in filename

                xx = find(mask == filesep);
                dd = mask(1:max(xx));
                ff = mask(max(xx)+1:end);
                mask = [dd 'rs_' ff];
            elseif use_res == 0
                % use the original mask
                mask = strtrim(masks(m,:));
            end


            vm = spm_read_vols(spm_vol(mask)); % one mask at a time
            [x1 y1 z1] = size(vm);
            vm2 = vm(:);
            vm2(isnan(vm2)) = 0;
            
            % 2013.11.07 - not sure why this is here!
            if choice == 8
                % threshold the image?
                if im_type == 1
                    vm2 = vm2 >= im_thr; % now binarized at this threshold
                elseif im_type == 2
                    % OK
                end
            end

        % x-axis orientation
        namemask = spm_vol(mask);
        mmat=namemask.('mat');
        mflip = sign(mmat(1,1));    % -1 means flipped along x-axis

        end 


        % ******************
        % file loop

        % set colors
        cc = distinguishable_colors(numf); % http://www.mathworks.com/matlabcentral/fileexchange/29702

        for f = 1:numf
            if f == 1
                maxY = NaN; % leave here
                minX = NaN;
                maxX = NaN;
            end
        file = strtrim(files(f,:));  % strtrim removes white spaces in cases where file names have different lengths 

        % read in the image
         namefile = spm_vol(file);
         namef2 = namefile.('fname');
         namef = dir(namef2);
         namef = namef.('name');
         vf = spm_read_vols(spm_vol(file));
         

        vmat=namefile.('mat');
        fflip = sign(vmat(1,1));    % -1 means flipped along x-axis

        origin2 = inv(vmat);
        origin2 = origin2(1:3,4);
        origin2 = origin2';
        origin(f,:) = origin2; % rje confirmed this is correct
        vox_dims(f,:) = abs([vmat(1,1) vmat(2,2) vmat(3,3)]);
        dims(f,:) = size(vf);
        [x2 y2 z2] = size(vf);

        % now x-flip the file if needed
        if use_mask == 1
            if mflip == fflip
               % no need to flip
            else
               vf = flipdim(vf,1); 
            end
        end

         vf2 = vf(:);
         vf2(isnan(vf2)) = 0;


        % === get the file values
        if use_mask == 1
            fvals = vf2(vm2 ~= 0);  % retain file values where the *mask* is not 0
        elseif use_mask == 2
            fvals = vf2(vf2 ~= 0);  % retain all non-zero fvals
        end

        % if there is only one MASK, save the file values

        if f == 1 
            if numm == 1
                fvals_all = nan(numel(fvals),numf);
            elseif numm > 1
                fvals_all = nan;
            end
        end

        % now, make sure we get rid of zeros in fvals (e.g., if the file were
        % originally a binarized image
        fvals = fvals(fvals ~= 0);

        % also, we may be dealing with small residual values
        if minV == 1
            fvals = fvals(fvals >= minT);
        elseif minV == 2
            % retain all current fvals
        end

        % note: fvals may be empty after this
        if isempty(fvals) 
            fvals = NaN; % so that mean and min etc will return NaN
            numvox(f,m) = 0;
        else
            numvox(f,m) = numel(fvals);

        end

        % percentiles
        if calc_perc == 1
            if prc_ch == 1
                % use matlab percentile
                percs_data(f,1:nperc,m) = prctile(fvals,perc(f,:));
            elseif prc_ch == 2
                % use NIST percentile
                percs_data(f,1:nperc,m) = prctile_nist(fvals,perc(f,:))';
            end
            
            % KDE method by RJE
            kde_res = prctile_kde(fvals,'b',perc(f,:),0); % this function uses 0 to 100 as input value
            percs_kde(f,1:nperc,m) = kde_res.kde_percentiles;
        elseif calc_perc == 0
            % do nothing; e.g., can't do percentiles if the matlab "stats" toolbox is not
            % installed
        end


        mins(f,m) = min(fvals); 
        means(f,m) = mean(fvals); 
        means_pos(f,m) = mean(fvals(fvals>0));  % mean of positive voxels
        maxs(f,m) = max(fvals);
        stds(f,m) = std(fvals); 
        stds(f,m) = std(fvals);    

        %if exist('skewness') == 2  
            skews(f,m) = skewness_vec(fvals); % skewness_vec is by RJE 
        %else
            %skews(f,m) = NaN; 
        %end

        sums(f,m) = sum(fvals);
        fzeros(f,m) = sum(fvals == 0);
        if numm == 1
            fvals_all(1:numel(fvals),f) = fvals;
        end
        
        if choice == 1
            % KDE plot: [bandwidth,density,xmesh,cdf]=kde(data,n,MIN,MAX)
            [bandwidth,density,xmesh]=kde(fvals,kde_mesh); % use defaults for KDE for now
            
            figure(101) 
            plot(xmesh,density,'LineWidth',1,'color',cc(f,:),'DisplayName',namef)
            Xvals(:,f) = xmesh(:);
            Yvals(:,f) = density(:);
            
            
            % store the extreme value across files
            maxY = max(maxY,max(density(:)));
            minX = min(minX,min(xmesh(:)));
            maxX = max(maxX,max(xmesh(:)));
            
            hold on
            
            % if the last file ...
            if f == numf
                % adjust the global axis
                axis([minX maxX 0 maxY*1.05])
            end
            
        elseif choice == 2
            % === now do the histogram
            figure(102)
            if commonX == 2
               % xvals defined already
            elseif commonX == 1
                xvals = linspace(min(fvals), max(fvals),nbins);  % note: this will ADJUST with every file
            end
            yvals = histc(fvals,xvals);
            plot(xvals,yvals,'LineWidth',1,'color',cc(f,:),'DisplayName',namef);
            hold on  

            Xvals(:,f) = xvals(:);
            Yvals(:,f) = yvals(:);
            

        end

        end % file loop


    end

    if choice == 1 || choice == 2
        if choice == 1
            lt = title('KDE of probability density');
            ly = ylabel('Probability density');
            
        elseif choice == 2
            lt = title('Histogram outline');
            ly = ylabel('Count'); 
        end

        lx = xlabel('Parameter value'); whitebg('w'); box on;
        set(lx,'Interpreter','none'); set(ly,'Interpreter','none'); set(lt,'Interpreter','none'); 
        legend1 = legend('show');
        set(legend1,'Interpreter','none','Location','NorthEastOutside');
        set(gcf, 'Color', 'w');
        
        hold off

    end
end

% =========================================
if choice == 8
        % read in the cluster / thresholded statistic image
        if use_res == 1
            % use the resliced image
            mask = strtrim(masks(1,:)); % cut out text whitespace in filename
            
            xx = find(mask == filesep);
            dd = mask(1:max(xx));
            ff = mask(max(xx)+1:end);
            mask = [dd 'rs_' ff];
        elseif use_res == 0
            % use the original mask
            mask = strtrim(masks(1,:));
        end
        
        % x-axis orientation
        namemask = spm_vol(mask);
        mmat=namemask.('mat');
        mflip = sign(mmat(1,1));    % -1 means flipped along x-axis
        
        vm = spm_read_vols(spm_vol(mask)); % one mask at a time
        [x1 y1 z1] = size(vm);
        vm2 = vm(:);
        vm2(isnan(vm2)) = 0;       

    cls_vox(2:end,1) = avals;
    cls_vox(1,2:end) = mvals;

    for m = 1:numm % for each cluster ...
        mtemp = vm2 == mvals(m); % the value of the current cluster

        %for f = 1:numf % only one atlas ...
            if m == 1 % only need to do this once
                file = strtrim(files(1,:));  % strtrim removes white spaces in cases where file names have different lengths 

                % read in the image
                 namefile = spm_vol(file);
                 namef2 = namefile.('fname');
                 namef = dir(namef2);
                 namef = namef.('name');
                 vf = spm_read_vols(spm_vol(file));

                vmat=namefile.('mat');
                fflip = sign(vmat(1,1));    % -1 means flipped along x-axis

                origin2 = inv(vmat);
                origin2 = origin2(1:3,4);
                origin2 = origin2';
                origin = origin2; % rje confirmed this is correct
                dims = size(vf);
                [x2 y2 z2] = size(vf);

                % now x-flip the file if needed
                if mflip == fflip
                    % OK
                else
                   vf = flipdim(vf,1); 
                end

                 vf2 = vf(:);
                 vf2(isnan(vf2)) = 0;
            end % read in the atlas

                % === apply the cluster mask to the atlas

                fvals = vf2(mtemp == 1);  % retain file values where the *mask* is not 0
                %fvals = fvals(fvals ~= 0); we need to retain 0 values to
                %count voxels outside of atlas regions


                % note: fvals may be empty after this
                if isempty(fvals) 
                    fvals = NaN; % so that mean and min etc will return NaN
                else
                    % now let's tally the regions
                     for a = 1:nregs
                         cls_vox(a+1,m+1) = sum(fvals == avals(a)); % num of voxels with this value
                     end    

                end

        %end
    end % mask loop
    
% now eliminate atlas regions that are empty

cls_voxT = cls_vox(2:end,2:end);
cls_sum = sum(cls_voxT,1); % column sum
at_sum = sum(cls_voxT,2); % row sum

% now get percentages
cls_voxP = bsxfun(@rdivide,cls_voxT,cls_sum) * 100;

cls_perc = cls_vox;
cls_perc(2:end,2:end) = cls_voxP;

at_keep = find(at_sum > 0) + 1; % + 1 restores the proper index
at_keep(numel(at_keep)+1) = 1; % also need the fist row
at_keep = sort(at_keep);

cls_vox = cls_vox(at_keep,:);
cls_perc = cls_perc(at_keep,:);


fvals_all = NaN;

fprintf(['\n Finished.\n',...
         '\n Cluster parcellation voxel counts are stored in variable "fstats.cls_vox".',... 
         '\n    * The first column indicates the numbered atlas region.',... 
         '\n    * Activation clusters are parcellated in separate columns, with row 1 indicating the original cluster number. ',...
         '\n Cluster parcellation percentages are stored in variable "fstats.cls_perc"',...
         '\n    * The table structure corresponds to that of fstats.cls_vox.',...
         '\n The cluster sums may be called via "fstats.cls_sum".\n']);

end % choice == 8

if choice == 1 || choice == 2 || choice == 6
    cls_sum = NaN;
    fprintf(['\n Finished.',... 
            '\n The variable "fstats" contains parameter values within each ROI [columns] from each image [rows].\n']);

end

numvox = round(numvox); % to get into integers

fstats.matrix_dims = dims;
fstats.matrix_origins = origin;
fstats.vox_dims = vox_dims;
fstats.ROI_voxels = numvox;
fstats.mins = mins;
fstats.means = means;
fstats.means_pos = means_pos;
fstats.maxs = maxs;
fstats.stds = stds;
fstats.skews = skews;
%fstats.sums = sums;
fstats.zeros = fzeros;
fstats.percs_data = percs_data;
fstats.percs_kde = percs_kde;
fstats.vals = fvals_all;

if choice == 1 || choice == 2
    fstats.plot_Xvals = Xvals;
    fstats.plot_Yvals = Yvals;
elseif choice == 8
    fstats.cls_vox = cls_vox;
    fstats.cls_sum = cls_sum;
    fstats.cls_perc = cls_perc;
end

% finally, delete the resliced masks
cd(maskdir)
delete('rs_*')

