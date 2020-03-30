function quality_report_funct(varargin)
% QUALITY_REPORT_FUNCT produces the user-specified functional MRI quality
% metric report
% Note: Functionals must be in the form of 4D NIfTI's (.nii or .nii.gz)
% 
% INPUTS (all are specified as name-value pairs):
% 
% REQUIRED:
% 
% 'parent_folder',      full path to the main directory in which data is
%                       stored
% 
% 'funct_prefix',       a string denoting the file prefix shared by all
%                       functionals (the parent_folder will be searched for
%                       files beginning with this string)
% 
% OPTIONAL:
% 
% 'output_dir',         a string denoting either a full path or a partial
%                       path where quality report outputs are to be
%                       written; if 'output_dir' is a partial path, the
%                       specified directory will be created within the same
%                       folder in which a given functional file is found
%   
% 'funct_ext',          a string denoting the file extension of functionals
%                       e.g., '.nii','.nii.gz'
%                       Default: '.nii'
% 
% 'exclude_str',        a string which if found within a file name will
%                       cause that file to be excluded; CASE SENSITIVE
%                       note: this can also be a cell array containing
%                       multiple strings to exclude
% 
% 'metric',             vector of 1/0 integers or booleans specifying which
%                       quality metrics are desired; order-specific options
%                       are [tSNR, SFNR, tSBNR, SNR-Funct, Mean Power];
%                       thus, if one desires tSNR and SFNR, the input
%                       should be specified as, <'metric',[1, 1, 0, 0, 0]>
%                       Default: [1, 0, 0, 0, 0]
% 
% 'qc_orientation',     a vector of 1/0 integers or booleans specifying 
%                       which orientation(s) user would like to output 
%                       slice mosaic .png files of quality metrics; 
%                       order-specific orientation options are: 
%                       [axial, sagittal, coronal]. Thus, if one desires 
%                       sagittal slices only, the input should be 
%                       <'qc_orientation',[0,1,0]>
%                       Default: [1,0,0]
% 
% 'startSlice',         a vector of integers specifying the first slice in
%                       each dimension that contains desired structures
%                       (e.g., brain); e.g., [36,1,45] would denote to
%                       begin with slice 36 for axial, 1 for sagittal, and
%                       45 for coronal
%                       default: an automatic estimate for this
% 
% 'endSlice',           a vector of integers specifying the last slice in
%                       each dimension that contains desired structures
%                       (e.g., brain); e.g., [36,24,45] would denote to
%                       end with slice 36 for axial, 24 for sagittal, and
%                       45 for coronal
%                       default: an automatic estimate for this
% 
% 'write_qc_volume',    boolean: write a 3D NIfTI for requested quality
%                       metrics? (Default: False)
% 
% 'slice_by_time',      a vector of 1/0 integers or booleans specifying 
%                       whether to plot signal for a given slice across 
%                       time; order-specific orientation options are: 
%                       [axial, sagittal, coronal]. Thus, if one desires to
%                       visualize a sagittal slice and a coronal slice
%                       across time, the input should be:
%                       <'slice_by_time',[0,1,1]>
%                       Default: [1,0,0]
% 
% 'which_slice_by_time', a vector of positive integers specifying which 
%                        slice is desired for 'slice_by_time' plots; Each 
%                        element of the vector corresponds to the slice 
%                        requested for each respective dimension as in
%                        'slice_by_time' input. Thus, if one desires
%                        sagittal slice 23 and coronal slice 47, the input
%                        should read: [NaN, 23, 47]. If an orientation was
%                        0/false for 'slice_by_time', any input for that
%                        orientation here will be ignored (thus, NaN input
%                        is not necessary). 
%                        Default: approximate middle slice
% 
% 'startTime',          integer: the first TR to use for the slice_by_time 
%                       mosaic(s) (default: first TR)
% 
% 'endTime',            integer: the last TR to use for the slice_by_time 
%                       mosaic(s) (default: last TR)
% 
% 'masking',            integer specifying type of masking: 1 = intensity
%                       thresholding (default), 2 = image region (a central
%                       cubic region is taken to contain signal, and outer
%                       edges of image are considered as background noise
% 
% 'printRes',           integer specifying resolution of .png outputs in
%                       dpi
%                       Default: 150
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Elliot Layden, 2018, The University of Chicago
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify Function Path and Add Helper Scripts:
script_fullpath = mfilename('fullpath');
[script_path,~,~] = fileparts(script_fullpath);
addpath(genpath(script_path))

% Retrieve name-value pair inputs:
inputs = varargin;
parsed = struct('parent_folder',[],'funct_prefix',[],'output_dir',[],...
    'funct_ext','.nii','exclude_str',[],'metric',[1, 0, 0, 0, 0],...
    'qc_orientation',[1,0,0],'write_qc_volume',0,'slice_by_time',[1,0,0],...
    'which_slice_by_time',[0,0,0],'masking',[],'printRes',150,...
    'startSlice',[],'endSlice',[],'startTime',[],'endTime',[]); 
poss_input = {'parent_folder','funct_prefix','output_dir','funct_ext',...
    'exclude_str','metric','qc_orientation','write_qc_volume','slice_by_time',...
    'which_slice_by_time','masking','printRes','startSlice','endSlice',...
    'startTime','endTime'}; 
input_ind = zeros(1,length(poss_input));
for i = 1:length(poss_input)
    j = find(strcmp(poss_input{i},inputs));
    if ~isempty(j)
        input_ind(i) = j;
        input1 = inputs{input_ind(i)+1};
        parsed.(poss_input{i}) = input1;
    end
end

if isempty(parsed.which_slice_by_time)
    parsed.which_slice_by_time = zeros(1,3);
end

% Search through all subdirectories of parent_folder, noting functs:
allSubDir = strsplit(genpath(parsed.parent_folder),';');
for i = 1:length(allSubDir)
    if isdir(allSubDir{i}) %#ok
        listing = dir(fullfile(allSubDir{i},[parsed.funct_prefix,'*',parsed.funct_ext]));
        if ~isempty(listing)
            for j = 1:length(listing)
                dontExclude = checkExclude(listing(j).name);
                if dontExclude
                    try
                        [~,fname] = fileparts(listing(j).name);
                        try
                            funct = load_nii(fullfile(allSubDir{i},fullfile(listing(j).name)));
                            writeType = 1;
                        catch
                            funct = load_untouch_nii(fullfile(allSubDir{i},fullfile(listing(j).name)));
                            writeType = 2;
                            warning(['Non-orthogonal shearing detected in affine matrix for ',fullfile(allSubDir{i},fullfile(listing(j).name))]) 
                        end
                        % For now, just do simple tSNR:
                        functDim = funct.hdr.dime.dim(2:5);
                        tSNRstat = nanmean(funct.img,4) ./ std(funct.img,0,4,'omitnan');
                        
                        % Get 3D funct:
                        funct3D = funct; funct3D.img = mean(funct.img,4);
                        funct3D.hdr.dime.dim(1) = 3;
                        funct3D.hdr.dime.dim(5) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% QC Metric Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        tSNR = funct3D; tSNR.img = tSNRstat;
                        statName = 'tSNR';
                        tSNRVec = tSNRstat(tSNRstat(:)>0);
                        sort_tSNR = sort(tSNRVec);
                        climMin = sort_tSNR(round(.2*numel(tSNRVec)));
                        climMax = sort_tSNR(round(.999*numel(tSNRVec)));
                        clim = [climMin, climMax];
                        
                        % Axial
                        if parsed.qc_orientation(1)
                            orientationName = 'Axial';
                            if isempty(parsed.startSlice) || parsed.startSlice<1 || parsed.startSlice>functDim(3)
                                botLimit = round(functDim(3)*.25);
                            else
                                botLimit = parsed.startSlice(1);
                            end
                            if isempty(parsed.endSlice) || parsed.endSlice<1 || parsed.endSlice>functDim(3)
                                topLimit = round(functDim(3)*.72);
                            else
                                topLimit = parsed.endSlice(1);
                            end
                            sliceRange = topLimit-botLimit+1;
                            % Determine slices and dimension for mosaic:
                            if (sliceRange>=15) 
                                slices = round(linspace(botLimit,topLimit,15));
                                axesDim = [3,5];
                            elseif (sliceRange>=10) 
                                slices = round(linspace(botLimit,topLimit,9));
                                axesDim = [2,5];
                            elseif (sliceRange>=5) 
                                slices = round(linspace(botLimit,topLimit,5));
                                axesDim = [1,5];
                            else
                                slices = botLimit:topLimit;
                                axesDim = [1,slices+1];
                            end
                            if isempty(parsed.output_dir) % did user specify a path for output?
                                [handles] = neuroimage_editor_mosaic('background',funct3D,'overlay',tSNR,'slices',slices,'dimension',3,'slice_label_pos',2,'colorbar',1,...
                                    'colormap','jet','title',[statName,'_',orientationName],'slice_locator',0,'print_res',parsed.printRes,'overlay_clim',clim,...
                                    'figure_pos',[.1,.038,.83,.91],'axes_dim',axesDim,'print',fullfile(allSubDir{i},[fname,'_',statName,'_',orientationName,'.png']));
                            elseif exist(parsed.output_dir,'dir')==7 % did user specify full path?
                                [handles] = neuroimage_editor_mosaic('background',funct3D,'overlay',tSNR,'slices',slices,'dimension',3,'slice_label_pos',2,'colorbar',1,...
                                    'colormap','jet','title',[statName,'_',orientationName],'slice_locator',0,'print_res',parsed.printRes,'overlay_clim',clim,...
                                    'figure_pos',[.1,.038,.83,.91],'axes_dim',axesDim,'print',fullfile(parsed.output_dir,[fname,'_',statName,'_',orientationName,'.png']));
                            else % did user specify a relative path to be created?
                                fpath = fullfile(allSubDir{i},parsed.output_dir);
                                mkdir(fpath);
                                [handles] = neuroimage_editor_mosaic('background',funct3D,'overlay',tSNR,'slices',slices,'dimension',3,'slice_label_pos',2,'colorbar',1,...
                                    'colormap','jet','title',[statName,'_',orientationName],'slice_locator',0,'print_res',parsed.printRes,'overlay_clim',clim,...
                                    'figure_pos',[.1,.038,.83,.91],'axes_dim',axesDim,'print',fullfile(fpath,[fname,'_',statName,'_',orientationName,'.png']));
                            end
                            close(handles.figure); 
                        end
                        
                        % Sagittal
                        if parsed.qc_orientation(2)
%                            imagesc(flipud(squeeze(funct.img(40,:,:,5))'))
                            orientationName = 'Sagittal';
                            if isempty(parsed.startSlice) || length(parsed.startSlice)<2
                                botLimit = round(functDim(2)*.22);
                            else
                                botLimit = parsed.startSlice(1);
                            end
                            if isempty(parsed.endSlice) || length(parsed.endSlice)<2
                                topLimit = round(functDim(2)*.78);
                            else
                                topLimit = parsed.endSlice(1);
                            end
                            sliceRange = topLimit-botLimit+1;
                            % Determine slices and dimension for mosaic:
                            if (sliceRange>=15) 
                                slices = round(linspace(botLimit,topLimit,15));
                                axesDim = [3,5];
                            elseif (sliceRange>=10) 
                                slices = round(linspace(botLimit,topLimit,9));
                                axesDim = [2,5];
                            elseif (sliceRange>=5) 
                                slices = round(linspace(botLimit,topLimit,5));
                                axesDim = [1,5];
                            else
                                slices = botLimit:topLimit;
                                axesDim = [1,slices+1];
                            end
                            if isempty(parsed.output_dir) % did user specify a path for output?
                                [handles] = neuroimage_editor_mosaic('background',funct3D,'overlay',tSNR,'slices',slices,'colorbar',1,...
                                    'figure_pos',[.1,.038,.83,.91],'dimension',2,'slice_label_pos',2,'overlay_clim',clim,... % 'rotate',-90,
                                    'colormap','jet','title',[statName,'_',orientationName],'slice_locator',0,'print_res',parsed.printRes,...
                                    'axes_dim',axesDim,'print',fullfile(allSubDir{i},[fname,'_',statName,'_',orientationName,'.png']));
                            elseif exist(parsed.output_dir,'dir')==7 % did user specify full path?
                                [handles] = neuroimage_editor_mosaic('background',funct3D,'overlay',tSNR,'slices',slices,'colorbar',1,...
                                    'figure_pos',[.1,.038,.83,.91],'dimension',2,'slice_label_pos',2,'overlay_clim',clim,... % 'rotate',-90,
                                    'colormap','jet','title',[statName,'_',orientationName],'slice_locator',0,'print_res',parsed.printRes,...
                                    'axes_dim',axesDim,'print',fullfile(parsed.output_dir,[fname,'_',statName,'_',orientationName,'.png']));
                            else % did user specify a relative path to be created?
                                fpath = fullfile(allSubDir{i},parsed.output_dir); mkdir(fpath);
                                [handles] = neuroimage_editor_mosaic('background',funct3D,'overlay',tSNR,'slices',slices,'colorbar',1,...
                                    'figure_pos',[.1,.038,.83,.91],'dimension',2,'slice_label_pos',2,'overlay_clim',clim,... % 'rotate',-90,
                                    'colormap','jet','title',[statName,'_',orientationName],'slice_locator',0,'print_res',parsed.printRes,...
                                    'axes_dim',axesDim,'print',fullfile(fpath,[fname,'_',statName,'_',orientationName,'.png']));
                            end
                            close(handles.figure); 
                        end
                        
                        % Coronal
                        if parsed.qc_orientation(3)
                            % imagesc(flipud(squeeze(funct.img(:,30,:,5))'))
                            orientationName = 'Coronal';
                            if isempty(parsed.startSlice) || length(parsed.startSlice)<3
                                botLimit = round(functDim(1)*.15);
                            else
                                botLimit = parsed.startSlice(1);
                            end
                            if isempty(parsed.endSlice) || length(parsed.endSlice)<3
                                topLimit = round(functDim(1)*.75);
                            else
                                topLimit = parsed.endSlice(1);
                            end
                            sliceRange = topLimit-botLimit+1;
                            % Determine slices and dimension for mosaic:
                            if (sliceRange>=15) 
                                slices = round(linspace(botLimit,topLimit,15));
                                axesDim = [3,5];
                            elseif (sliceRange>=10) 
                                slices = round(linspace(botLimit,topLimit,9));
                                axesDim = [2,5];
                            elseif (sliceRange>=5) 
                                slices = round(linspace(botLimit,topLimit,5));
                                axesDim = [1,5];
                            else
                                slices = botLimit:topLimit;
                                axesDim = [1,slices+1];
                            end
                            if isempty(parsed.output_dir) % did user specify a path for output?
                                [handles] = neuroimage_editor_mosaic('background',funct3D,'overlay',tSNR,'slices',slices,'colorbar',1,...
                                    'figure_pos',[.1,.038,.83,.91],'dimension',1,'slice_label_pos',2,'overlay_clim',clim,...
                                    'colormap','jet','title',[statName,'_',orientationName],'slice_locator',0,'print_res',parsed.printRes,...
                                    'axes_dim',axesDim,'print',fullfile(allSubDir{i},[fname,'_',statName,'_',orientationName,'.png']));
                            elseif exist(parsed.output_dir,'dir')==7 % did user specify full path?
                                [handles] = neuroimage_editor_mosaic('background',funct3D,'overlay',tSNR,'slices',slices,'colorbar',1,...
                                    'figure_pos',[.1,.038,.83,.91],'dimension',1,'slice_label_pos',2,'overlay_clim',clim,...
                                    'colormap','jet','title',[statName,'_',orientationName],'slice_locator',0,'print_res',parsed.printRes,...
                                    'axes_dim',axesDim,'print',fullfile(parsed.output_dir,[fname,'_',statName,'_',orientationName,'.png']));
                            else % did user specify a relative path to be created?
                                fpath = fullfile(allSubDir{i},parsed.output_dir); mkdir(fpath);
                                [handles] = neuroimage_editor_mosaic('background',funct3D,'overlay',tSNR,'slices',slices,'colorbar',1,...
                                    'figure_pos',[.1,.038,.83,.91],'dimension',1,'slice_label_pos',2,'overlay_clim',clim,...
                                    'colormap','jet','title',[statName,'_',orientationName],'slice_locator',0,'print_res',parsed.printRes,...
                                    'axes_dim',axesDim,'print',fullfile(fpath,[fname,'_',statName,'_',orientationName,'.png']));
                            end
                            close(handles.figure); 
                        end
                        
%%%%%%%%%%%%%%%%%%%%%%%%%% Slice by Time Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%
                        if any(parsed.slice_by_time)
                            % Get Time Range:
                            if ~isempty(parsed.startTime) && parsed.startTime>0 && parsed.startTime<=functDim(4)
                                startTime = parsed.startTime;
                            else
                                startTime = 1;
                            end
                            if ~isempty(parsed.endTime) && parsed.endTime>0 && parsed.endTime<=functDim(4) && parsed.endTime>startTime
                                endTime = parsed.endTime;
                            else
                                endTime = functDim(4);
                            end
                            timeRange = endTime-startTime+1;
                            % Determine how many times for mosaic:
                            if (timeRange>=15) 
                                times = round(linspace(startTime,endTime,15));
                                axesDim = [3,5];
                            elseif (timeRange>=10) 
                                times = round(linspace(startTime,endTime,9));
                                axesDim = [2,5];
                            elseif (timeRange>=5) 
                                times = round(linspace(startTime,endTime,5));
                                axesDim = [1,5];
                            else
                                times = startTime:endTime;
                                axesDim = [1,times+1];
                            end
                            timesLabels = cell(1,length(times));
                            for ixx = 1:length(times)
                                timesLabels{ixx} = num2str(times(ixx));
                            end
                        end
                        
                        if parsed.slice_by_time(1)
                            orientationName = 'Axial';
                            % Get Slice:
                            if ~isempty(parsed.which_slice_by_time) && parsed.which_slice_by_time(1)>0 && parsed.which_slice_by_time(1)<=functDim(3)
                                whichSlice = parsed.which_slice_by_time(1);
                            else
                                whichSlice = round(median(1:functDim(3)));
                            end
                            % Form 3D image:
                            timesImg = funct3D;
                            timesImg.img = squeeze(funct.img(:,:,whichSlice,times));
                            timesImg.hdr.dime.dim(2:4) = size(timesImg.img);
                            % Plot Mosaic:
                            if isempty(parsed.output_dir) % did user specify a path for output?
                                [handles] = neuroimage_editor_mosaic('background',timesImg,'slices',1:length(times),'dimension',3,'slice_label_pos',2,'colorbar',1,...
                                    'colormap','jet','title',['functional_TRs','_',orientationName],'slice_locator',0,'print_res',parsed.printRes,'custom_slice_labels',timesLabels,...
                                    'figure_pos',[.1,.038,.83,.91],'axes_dim',axesDim,'print',fullfile(allSubDir{i},[fname,'_','functional_TRs','_',orientationName,'.png']));
                            elseif exist(parsed.output_dir,'dir')==7 % did user specify full path?
                                [handles] = neuroimage_editor_mosaic('background',timesImg,'slices',1:length(times),'dimension',3,'slice_label_pos',2,'colorbar',1,...
                                    'colormap','jet','title',['functional_TRs','_',orientationName],'slice_locator',0,'print_res',parsed.printRes,'custom_slice_labels',timesLabels,...
                                    'figure_pos',[.1,.038,.83,.91],'axes_dim',axesDim,'print',fullfile(parsed.output_dir,[fname,'_','functional_TRs','_',orientationName,'.png']));
                            else % did user specify a relative path to be created?
                                fpath = fullfile(allSubDir{i},parsed.output_dir); mkdir(fpath);
                                [handles] = neuroimage_editor_mosaic('background',timesImg,'slices',1:length(times),'dimension',3,'slice_label_pos',2,'colorbar',1,...
                                    'colormap','jet','title',['functional_TRs','_',orientationName],'slice_locator',0,'print_res',parsed.printRes,'custom_slice_labels',timesLabels,...
                                    'figure_pos',[.1,.038,.83,.91],'axes_dim',axesDim,'print',fullfile(fpath,[fname,'_','functional_TRs','_',orientationName,'.png']));
                            end
                            close(handles.figure); 
                        end
                        
                        if parsed.slice_by_time(2)
                            orientationName = 'Sagittal';
                            % Get Slice:
                            if ~isempty(parsed.which_slice_by_time) && parsed.which_slice_by_time(2)>0 && parsed.which_slice_by_time(2)<=functDim(2)
                                whichSlice = parsed.which_slice_by_time(2);
                            else
                                whichSlice = round(median(1:functDim(2)));
                            end
                            % Form 3D image:
                            timesImg = funct3D;
                            timesImg.img = squeeze(funct.img(whichSlice,:,:,times));
                            timesImg.hdr.dime.dim(2:4) = size(timesImg.img);
                            % Plot Mosaic:
                            if isempty(parsed.output_dir) % did user specify a path for output?
                                [handles] = neuroimage_editor_mosaic('background',timesImg,'slices',1:length(times),'dimension',3,'slice_label_pos',2,'colorbar',1,...
                                    'colormap','jet','title',['functional_TRs','_',orientationName],'slice_locator',0,'print_res',parsed.printRes,'custom_slice_labels',timesLabels,...
                                    'figure_pos',[.1,.038,.83,.91],'axes_dim',axesDim,'print',fullfile(allSubDir{i},[fname,'_','functional_TRs','_',orientationName,'.png']));
                            elseif exist(parsed.output_dir,'dir')==7 % did user specify full path?
                                [handles] = neuroimage_editor_mosaic('background',timesImg,'slices',1:length(times),'dimension',3,'slice_label_pos',2,'colorbar',1,...
                                    'colormap','jet','title',['functional_TRs','_',orientationName],'slice_locator',0,'print_res',parsed.printRes,'custom_slice_labels',timesLabels,...
                                    'figure_pos',[.1,.038,.83,.91],'axes_dim',axesDim,'print',fullfile(parsed.output_dir,[fname,'_','functional_TRs','_',orientationName,'.png']));
                            else % did user specify a relative path to be created?
                                fpath = fullfile(allSubDir{i},parsed.output_dir); mkdir(fpath);
                                [handles] = neuroimage_editor_mosaic('background',timesImg,'slices',1:length(times),'dimension',3,'slice_label_pos',2,'colorbar',1,...
                                    'colormap','jet','title',['functional_TRs','_',orientationName],'slice_locator',0,'print_res',parsed.printRes,'custom_slice_labels',timesLabels,...
                                    'figure_pos',[.1,.038,.83,.91],'axes_dim',axesDim,'print',fullfile(fpath,[fname,'_','functional_TRs','_',orientationName,'.png']));
                            end
                            close(handles.figure); 
                            
                        end
                        
                        if parsed.slice_by_time(3)
                            orientationName = 'Coronal';
                            % Get Slice:
                            if ~isempty(parsed.which_slice_by_time) && parsed.which_slice_by_time(3)>0 && parsed.which_slice_by_time(3)<=functDim(1)
                                whichSlice = parsed.which_slice_by_time(3);
                            else
                                whichSlice = round(median(1:functDim(1)));
                            end
                            % Form 3D image:
                            timesImg = funct3D;
                            timesImg.img = squeeze(funct.img(:,whichSlice,:,times));
                            timesImg.hdr.dime.dim(2:4) = size(timesImg.img);
                            % Plot Mosaic:
                            if isempty(parsed.output_dir) % did user specify a path for output?
                                [handles] = neuroimage_editor_mosaic('background',timesImg,'slices',1:length(times),'dimension',3,'slice_label_pos',2,'colorbar',1,...
                                    'colormap','jet','title',['functional_TRs','_',orientationName],'slice_locator',0,'print_res',parsed.printRes,'custom_slice_labels',timesLabels,...
                                    'figure_pos',[.1,.038,.83,.91],'axes_dim',axesDim,'print',fullfile(allSubDir{i},[fname,'_','functional_TRs','_',orientationName,'.png']));
                            elseif exist(parsed.output_dir,'dir')==7 % did user specify full path?
                                [handles] = neuroimage_editor_mosaic('background',timesImg,'slices',1:length(times),'dimension',3,'slice_label_pos',2,'colorbar',1,...
                                    'colormap','jet','title',['functional_TRs','_',orientationName],'slice_locator',0,'print_res',parsed.printRes,'custom_slice_labels',timesLabels,...
                                    'figure_pos',[.1,.038,.83,.91],'axes_dim',axesDim,'print',fullfile(parsed.output_dir,[fname,'_','functional_TRs','_',orientationName,'.png']));
                            else % did user specify a relative path to be created?
                                fpath = fullfile(allSubDir{i},parsed.output_dir); mkdir(fpath);
                                [handles] = neuroimage_editor_mosaic('background',timesImg,'slices',1:length(times),'dimension',3,'slice_label_pos',2,'colorbar',1,...
                                    'colormap','jet','title',['functional_TRs','_',orientationName],'slice_locator',0,'print_res',parsed.printRes,'custom_slice_labels',timesLabels,...
                                    'figure_pos',[.1,.038,.83,.91],'axes_dim',axesDim,'print',fullfile(fpath,[fname,'_','functional_TRs','_',orientationName,'.png']));
                            end
                            close(handles.figure); 
                        end
                        
                        % 3D QC Metric Volume Output:
                        if parsed.write_qc_volume
                            if writeType==1
                                save_nii(tSNR,fullfile(allSubDir{i},fname,'_',statName));
                            elseif writeType==2
                                save_untouch_nii(tSNR,fullfile(allSubDir{i},fname,'_',statName));
                            end
                        end
                        
                    catch
                        warning(['Failed to generate quality report for ',fullfile(allSubDir{i},fullfile(listing(j).name))])
                    end
                end
            end
        end
    end
end

    function shouldRun = checkExclude(inputStr)
        if isempty(parsed.exclude_str) 
            shouldRun = true;
        elseif ischar(parsed.exclude_str)
            shouldRun = isempty(strfind(inputStr,parsed.exclude_str)); % check for exclude_str
        elseif iscell(parsed.exclude_str)
            check = zeros(1,length(parsed.exclude_str));
            for ix = 1:length(parsed.exclude_str)
                check(ix) = ~isempty(strfind(inputStr,parsed.exclude_str{ix})); 
            end
            shouldRun = ~any(check);
        end
    end

end