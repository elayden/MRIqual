function MRIqual
% MRIqual  % A GUI for converting DICOM's, running visual checks on
% structural & functional neuroimaging data, and checking various quality
% measures such as tSNR, SFNR, 
% (file types: .nii, .nii.gz, .img/.hdr), checking temporal SNR (tSNR) of 
% functional images, average power within resting-state frequencies, 
% ghost ratio of structural images, visualization of local tSNR or SNR,
% functional ROI correlation checking, frequency spectrum visualization.
% Various options for statistics output include .csv, .txt, & command line.
%   
% Author: Elliot Layden, The University of Chicago, 2016-2019
% 
% Usage: 
% To begin, simply type "MRIqual" into the command line.
% Alternatively, simply right-click on the MRIqual.m
% script, and choose Run.
% 
% FEATURES:
% 
% I. DICOM Conversion
%
% II. View/Edit 3D Image
%   This functionality calls neuroimage_editor, a major external function
%   used throughout the modules below. It enables manually scanning through
%   2D slices with arrow keys, zooming with the mouse scrollwheel, drawing
%   ROIs via manual tracing or shapes (ovals,circles,rectangles) which can
%   be propogated through slices, displaying statistical or ROI-based
%   overlays, and many other functions. A more thorough description can be
%   found within the neuroimage_editor.m file.
% 
% III. Structural Quality
%   Signal to Noise Ratio (SNR)
%       1. SNR_SD = mean voxel intensity in specified signal ROI divided by 
%           SD of voxel intensities within a noise (background) ROI,
%           multiplied by a correction factor
%       2. SNR_mean = mean voxel intensity in specified signal ROI divided 
%           by mean of voxel intensities within a noise (background) ROI,
%           multiplied by a correction factor
%       *Reference:
%             Dietrich, O., Raya, J. G., Reeder, S. B., Reiser, M. F., &
%               Schoenberg, S. O. (2007). Measurement of signal?to?noise 
%               ratios in MR images: Influence of multichannel coils, 
%               parallel imaging, and reconstruction filters. Journal of 
%               Magnetic Resonance Imaging, 26(2), 375-385. (see A11 & A12)
%   Visualize local SNR
% 
% IV. Functional Quality
% 
%   1. Temporal Signal-to-Noise Ratio (tSNR)
%       = mean voxel signal across time divided by SD across time
%       *References: 
%                  Triantafyllou, C., Hoge, R. D., Krueger, G., Wiggins, 
%                   C. J., Potthast, A., Wiggins, G. C., & Wald, L. L. (2005). 
%                   Comparison of physiological noise at 1.5 T, 3 T and 7 T 
%                   and optimization of fMRI acquisition parameters. 
%                   Neuroimage, 26(1), 243-250.
%                  Murphy, K., Bodurka, J., & Bandettini, P. A. (2007). 
%                   How long to scan? The relationship between fMRI temporal 
%                   signal to noise ratio and necessary scan duration. 
%                   Neuroimage, 34(2), 565-574. 
%                   <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2223273/>
% 
%   2. Temporal Signal-to-Background Noise Ratio (tSBNR)
%       = spatial average across a signal ROI of temporal mean voxel
%       signals / spatial average across temporal standard deviation of 
%       noise/background ROI voxels
%       -whereas tSNR divides mean signal by the temporal SD of the same
%       voxels, this measure compares the mean signal to the average
%       temporal SD of background/noise voxels, which may provide a more
%       accurate measure of temporal fluctuations due specifically to noise
% 
%   3. Signal-to-Fluctuation-Noise Ratio (SFNR)
%       = mean voxel signal across time divided by temporal standard deviation
%           of residuals obtained by fitting second-order polynomial to voxel
%           signal data
%       *Reference: 
%           Friedman, L., & Glover, G. H. (2006). Report on a multicenter 
%               fMRI quality assurance protocol. Journal of Magnetic 
%               Resonance Imaging, 23(6), 827-839.
% 
%   4. Signal-to-Noise Ratio (SNR-functional)
%       = mean voxel signal of signal ROI (averaged across time) divided by
%           the spatial SD of a noise ROI signal (averaged across time),
%           multiplied by a correction factor
%           Reference: 
%               http://wagerlab.colorado.edu/wiki/doku.php/help/fmri_quality_control_overview
% 
%   5. Calculate average power within a given frequency range
% 
%   Visualize local tSNR, SFNR, local power, & power spectrum
% 
%   Power Spectrum & Power Spectral Density Estimation
%       MRIqual implements two methods for obtaining these
%       quantities: 
%       1. Welch's periodogram is used to calculate the power spectrum,
%       utilizing a Hamming window with 50% overlap. Specifically, time
%       courses are by default divided into the longest segments possible
%       so as to yield close to, but not exceeding, 8 segments w/ 50%
%       overlap. Power spectral density estimates are then scaled by the 
%       equivalent noise bandwidth (ENBW) of the window (in hertz). To 
%       obtain the ENBW using Matlab, e.g., "enbw(hamming(N), fs)", where 
%       fs = sampling frequency, N = # of samples.  
%       2. The Multitaper method (Matlab: pmtm.m) is used to calculate the
%       power spectral density estimate using default parameters. The power
%       spectrum, scaling by ENBW, is not easily obtainable by this method.
% 
% Miscellaneous
%   *Determine ROIs using rectangular selection tool, manual tracing,
%       loading a mask image, or specifying an intensity threshold.
%   *ROIs may be propogated through slices or through slices and time
%   *Generate Report: Calculates all metrics for either a functional or
%       structural image, and saves output as a text file with accompanying
%       .png images.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize GUI Window & Data:

% Identify Function Path and Add Helper Scripts:
script_fullpath = mfilename('fullpath');
[script_path,~,~] = fileparts(script_fullpath);
addpath(genpath(script_path)) % adds dicm2nii, neuroimage_editor, etc.

% Check for Dependencies: neuroimage_editor, dicm2nii, load_untouch_nii
if ~exist('neuroimage_editor_tbx.m','file')
    errordlg('ERROR: neuroimage_editor_tbx.m not found. Image viewing functionality disabled.',...
      'Missing Dependency','modal');
end
if ~exist('nifti_studio.m','file')
    errordlg('ERROR: nifti_studio.m not found. Image viewing functionality disabled.',...
      'Missing Dependency','modal');
end
if ~exist('dicm2nii.m','file')
    errordlg('ERROR: dicm2nii.m not found. DICOM conversion disabled.',...
      'Missing Dependency','modal')
end
if ~exist('load_untouch_nii.m','file')
    errordlg('ERROR: NIFTI_Tools not found. Image loading disabled.',...
      'Missing Dependency','modal')
end

% Initialize Figure (& auto-detect screen size):
screen_res = get(0, 'MonitorPositions'); % get(0,'ScreenSize');
figure_pos = [0.36948, 0.2541, 0.25, 0.55];
if any(figure_pos<=0) % avoid ScreenSize errors
    figure_pos = [342,50,700,680]; 
    screen_res = [1,1,(figure_pos(3:4)./[.5,.86])];
    figure_pos = [.37*screen_res(3), .255*screen_res(4), ...
        .25*screen_res(3), .55*screen_res(4)];
end

main_figure = figure('units','normalized','Position',[0.36948, 0.2541, 0.25, 0.55],'MenuBar','none',...
    'Name','MRIqual','NumberTitle','off','Color',ones(1,3),...
    'Visible','on','doublebuffer','on','CloseRequestFcn',@mainfig_closereq);
    
% Try loading background image:
try 
    background = imread(fullfile(script_path,'res','background_brain.png'));
    imagesc(background(:,:,1)); colormap('gray'); caxis([0,260]);
    set(gca,'Position',[0,0,1,1],'TickLength',[0,0],'XTickMode','manual',...
        'YTickMode','manual','XTick',[],'YTick',[]); grid off; 
catch
end

% Intialize GUI Data:
guidata1 = guidata(main_figure);
guidata1.script_path = script_path;
guidata1.screen_res = screen_res;
guidata1.main_figure = main_figure;
guidata1.figure_pos = figure_pos;
guidata1.background_color = [.8,.88,.98];
guidata1.extensions={'*.img';'*.nii';'*nii.gz'};
guidata1.apply_header = 0;
guidata1.viewer_figure = 28.1773; % initialize false handle
guidata1.font_color = [.873,.546,.347]; % copper
% Directory Memory Initialization:
guidata1.first_dir = pwd;
guidata1.last_dicom_dir = pwd;
guidata1.last_output_dir = pwd;
guidata1.last_struct_dir = pwd;
guidata1.last_funct_dir = pwd;
guidata1.funct_path = [];
% Main Figures Handle Initialization:
guidata1.dicom_figure = 1.48437;
guidata1.structurals_figure = 2.48482;
guidata1.funct_figure = 3.484727;
guidata1.h = 10.18384; % non-handle for view image 
numvox_total = [];
all_sig = [];
all_detrended = [];
guidata1.SFNR_all = [];
guidata1.signal_amp_all = [];
% Structural Specific Initialization:
guidata1.struct_name = '';
guidata1.struct_name_string = '';
guidata1.structurals_state = 0;
guidata1.struct = [];
guidata1.struct_ROI = [];
guidata1.struct_mean_signal = [];
guidata1.struct_ROI_type = '';
guidata1.struct_SNR_sd = [];
guidata1.struct_SNR_sd_local = [];
guidata1.struct_SNR_mean = [];
guidata1.struct_SNR_mean_local = [];
guidata1.struct_report_options = ones(1,5);
guidata(main_figure,guidata1)
% Functional Initialization
funct_initialize

%% Initalize Main Buttons:

% Create I. "Convert DICOMs" Button:
guidata1.dicom_button = uicontrol(guidata1.main_figure,'style',...
    'togglebutton','Units','normalized','Position',[.275,.78,.45,.1],...
    'BackgroundColor', [0,0,0],'FontName','Segoe UI Black','FontSize',14,'FontWeight','normal','String',...
    'Convert DICOMs','callback',@dicoms_button_callback,'ForegroundColor',guidata1.font_color); 

% Create II. "View/Edit 3D Images" Button:
guidata1.viewEdit_button = uicontrol(guidata1.main_figure,'style',...
    'togglebutton','Units','normalized','Position',[.275,.56,.45,.1],...
    'BackgroundColor', [0,0,0],'FontName','Segoe UI Black','FontSize',14,'FontWeight','normal','String',...
    'View/Edit 3D Images','callback',@viewEdit_button_callback,'ForegroundColor',guidata1.font_color); 

% Create III. "Structural Quality" Button:
guidata1.structurals_button = uicontrol(guidata1.main_figure,'style',...
    'togglebutton','Units','normalized','Position',[.275,.34,.45,.1],...
    'BackgroundColor', [0,0,0],'FontName','Segoe UI Black','FontSize',14,'FontWeight','normal','String',...
    'Structural Quality','callback',@structurals_button_callback,'ForegroundColor',guidata1.font_color); 

% Create IV. "Functional Quality" Button:
guidata1.functionals_button = uicontrol(guidata1.main_figure,'style',...
    'togglebutton','Units','normalized','Position',[.275,.12,.45,.1],...
    'BackgroundColor', [0,0,0],'FontName','Segoe UI Black','FontSize',14,'FontWeight','normal','String',...
    'Functional Quality','callback',@functionals_button_callback,'ForegroundColor',guidata1.font_color); 

guidata(main_figure,guidata1)

function mainfig_closereq(hObject, eventdata, overlay_num) %#ok
    selection = questdlg('Are you sure you want to exit?','Exit MRIqual',...
      'Yes','No','Yes');
    switch selection
        case 'Yes'
            guidata1 = guidata(main_figure);
            if isgraphics(guidata1.dicom_figure)
                delete(guidata1.dicom_figure)
            end
            if isgraphics(guidata1.structurals_figure)
                delete(guidata1.structurals_figure)
            end
            if isgraphics(guidata1.funct_figure)
                delete(guidata1.funct_figure)
            end
            if isgraphics(guidata1.main_figure)
                delete(guidata1.main_figure)
            end
        case 'No'
          return;
    end
end

function funct_initialize
    % Retrieve Data
    guidata1 = guidata(main_figure);
    guidata1.numvoxel_limit = 500; % limit for plotting power spectrum, signals
    guidata1.chunk_num = 50000; % number of voxels to chunk by
    guidata1.funct_name = '';
    guidata1.funct_name_string = '';
    guidata1.functionals_state = 0;
    guidata1.funct = []; % 4D
    guidata1.funct3D = []; % 3D after temporal averaging
    guidata1.numvox = 0;
    guidata1.funct_dim = [];
    guidata1.multi_select = false;
    guidata1.n_4D = [];
    guidata1.funct_mean_signal = [];
    guidata1.funct_ROI_type = [];
    guidata1.ROI_ind = [];
    guidata1.opts = false(1,5);
    guidata1.ax1 = 28.3928;
    guidata1.ROI_x = [];
    guidata1.ROI_y = [];
    guidata1.ROI_z = [];
    guidata1.ROI_ind_back = [];
    guidata1.roi_sig = [];
    guidata1.roi_background = [];
    guidata1.roi_sig_filt = [];
    guidata1.sig_detrended = [];
    guidata1.background_detrended = [];
    guidata1.roi_background_filt = [];
    guidata1.mean_signal_amp = [];
    guidata1.h_plots = [34.1242,33.8974,12.2349,13.3593,14.4859,65.6456];
    guidata1.legend_str = {'Signal (Mean Centered)','Detrended (Linear)',...
        'Detrended (Quadratic)','Linear Trend','Quadratic Trend',...
        'Band-Pass Filtered'};
    guidata1.Fs = [];
    guidata1.HighPass = [];
    guidata1.LowPass = [];
    guidata1.report_dir = [];
    guidata1.funct_tSNR = [];
    guidata1.funct_tSNR_local = [];
    guidata1.tSBNR = [];
    guidata1.funct_tSBNR = [];
    guidata1.funct_SFNR = [];
    guidata1.use_linear = [];
    guidata1.funct_SFNR_local = [];
    guidata1.funct_SNR_sd = [];
    guidata1.funct_SNR_sd_local = [];
    guidata1.funct_power = [];
    guidata1.funct_power_local = [];
    guidata1.mean_background_amp = [];
    guidata1.funct_report_options = ones(1,7);
    guidata(main_figure,guidata1)
end

%% I. CONVERT DICOMS 

function dicoms_button_callback(hObject, eventdata, overlay_num) %#ok
    % Retrieve Data
    guidata1 = guidata(main_figure);
    % Determine State
    guidata1.dicoms_state = get(hObject,'Value');
    if guidata1.dicoms_state
        set(guidata1.viewEdit_button,'Value',false)
        set(guidata1.structurals_button,'Value',false);
        set(guidata1.functionals_button,'Value',false);
        if isgraphics(guidata1.viewer_figure)
            delete(guidata1.viewer_figure)
        end
        if isgraphics(guidata1.structurals_figure)
            delete(guidata1.structurals_figure)
        end
        if isgraphics(guidata1.funct_figure)
            delete(guidata1.funct_figure)
        end
    else
        if isgraphics(guidata1.dicom_figure)
            delete(guidata1.dicom_figure)
        end
        return;
    end
    % Generate Dicom Figure 
%     figure_pos1 = [.225*guidata1.screen_res(3), .45*guidata1.screen_res(4),...
%         .55*guidata1.screen_res(3), .18*guidata1.screen_res(4)];
    guidata1.dicom_figure = figure('units','normalized','Position',[.2245,.4491,.55,.18],'MenuBar','none',...
        'Name','Convert DICOMs','NumberTitle','off','Color',[0,0,0],...
        'Visible','on','doublebuffer','on','CloseRequestFcn',@dicom_closereq);
    % Create "Select DICOMs Folder(s)" Button:
    guidata1.select_dicoms_button = uicontrol(guidata1.dicom_figure,'style',...
        'pushbutton','Units','normalized','Position',[.01,.65,.21,.23],...
        'BackgroundColor', [0,0,0],'ForegroundColor',guidata1.font_color,...
        'FontSize',10,'FontWeight','bold','String',...
        'Select DICOM Folder(s)','callback',@select_dicoms_callback); 
    % Create "Select Output Folder(s)" Text Box:
    guidata1.select_dicoms_text = uicontrol(guidata1.dicom_figure,'style',...
        'edit','Units','normalized','Position',[.23,.65,.75,.23],...
        'BackgroundColor', [1,1,1]);
    % Create "Select Output Folder(s)" Button:
    guidata1.select_output_button = uicontrol(guidata1.dicom_figure,'style',...
        'pushbutton','Units','normalized','Position',[.01,.35,.21,.23],...
        'BackgroundColor', [0,0,0],'ForegroundColor',guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'Select Output Folder(s)','callback',@select_output_callback); 
    % Create "Select Output Folder(s)" Text Box:
    guidata1.select_output_text = uicontrol(guidata1.dicom_figure,'style',...
        'edit','Units','normalized','Position',[.23,.35,.75,.23],...
        'BackgroundColor', [1,1,1]);
    % Create "File Type:" Text
    uicontrol('Parent',guidata1.dicom_figure,'Style','text','String','File Type:',...
        'Units', 'normalized','FontSize',10,'FontName','Helvetica','Position',...
        [.02,.025,.1,.2],'BackgroundColor',[0,0,0],'ForegroundColor',guidata1.font_color,...
        'FontWeight','Bold','HorizontalAlignment','left');
    % Create File Type Selection Drop-Down:
    guidata1.filetype_dropdown = uicontrol('Parent',guidata1.dicom_figure,...
        'Style','popupmenu','String','.img/hdr|.nii|.nii.gz|3D.nii|3D.nii.gz',...
        'Units', 'normalized','FontSize',9,'FontName','Helvetica','Position',...
        [.12,.05,.1,.2],'BackgroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','center');
    % Create "Output File Name:" Text
    uicontrol('Parent',guidata1.dicom_figure,'Style','text','String','Output File Name:',...
        'Units', 'normalized','FontSize',10,'FontName','Helvetica','Position',...
        [.35,.025,.2,.2],'BackgroundColor',[0,0,0],'ForegroundColor',guidata1.font_color,...
        'FontWeight','Bold','HorizontalAlignment','left');
    % Create Output File Name Entry Box:
    guidata1.output_name_text = uicontrol(guidata1.dicom_figure,'style',...
        'edit','Units','normalized','Position',[.52,.08,.2,.18],...
        'BackgroundColor', [1,1,1]);  
    % Create "Convert" Button:
    guidata1.run_dicoms_button = uicontrol(guidata1.dicom_figure,'style',...
        'pushbutton','Units','normalized','Position',[.85,.08,.13,.18],...
        'BackgroundColor', guidata1.font_color,'ForegroundColor',[0,0,0] ,...
        'FontSize',10,'FontWeight','bold','String',...
        'Convert','callback',@run_dicoms_callback); 
    % Save Data
    guidata(main_figure,guidata1)
end

function dicom_closereq(hObject, eventdata, overlay_num) %#ok
    guidata1 = guidata(main_figure);
    guidata1.dicom_button.Value=0;
    if isgraphics(guidata1.dicom_figure)
        delete(guidata1.dicom_figure)
    end
    guidata(main_figure,guidata1)
end

function select_dicoms_callback(hObject, eventdata, overlay_num) %#ok
    guidata1 = guidata(main_figure);
    [guidata1.dicom_dir] = uigetdir(guidata1.last_dicom_dir,'Select DICOM Folder:');
    if guidata1.dicom_dir~=0
        [guidata1.last_dicom_dir,~,~] = fileparts(guidata1.dicom_dir);
        guidata1.select_dicoms_text.String = guidata1.dicom_dir;
    end
    guidata(main_figure,guidata1)
end

function select_output_callback(hObject, eventdata, overlay_num) %#ok
    guidata1 = guidata(main_figure);
    [guidata1.output_dir] = uigetdir(guidata1.last_output_dir,'Select Output Folder:');
    if guidata1.output_dir~=0
        [guidata1.last_output_dir,~,~] = fileparts(guidata1.output_dir);
        guidata1.select_output_text.String = guidata1.output_dir;
    end
    guidata(main_figure,guidata1)
end

function run_dicoms_callback(hObject, eventdata, overlay_num) %#ok
    guidata1 = guidata(main_figure);
    if isdir(guidata1.select_dicoms_text.String) %#ok
        if ~isdir(guidata1.select_output_text.String) %#ok
            mkdir(guidata1.select_output_text.String)
        end
        switch guidata1.filetype_dropdown.Value
            case 1
                out_type = 2; % .img/hdr 
            case 2
                out_type = 0; % .nii
            case 3
                out_type = 1; % .nii.gz
            case 4
                out_type = 4; % 3D.nii
            case 5
                out_type = 5; % 3D.nii.gz
        end
        dicm2nii(guidata1.select_dicoms_text.String, ...
            guidata1.select_output_text.String, out_type) 
        % Rename if Output Name Specified:
        S = load(fullfile(guidata1.select_output_text.String,'dcmHeaders.mat'));
        field_names = fieldnames(S.h); h = S.h; %#ok
        if ~isempty(guidata1.output_name_text.String)
            % Rename dcmHeaders.mat
            movefile(fullfile(guidata1.select_output_text.String,'dcmHeaders.mat'),...
                fullfile(guidata1.select_output_text.String,...
                [guidata1.output_name_text.String,'dcmHeaders.mat']))
            % Rename Image Files:
            listing = dir(fullfile(guidata1.select_output_text.String,[field_names{1},'*']));
            [~,fname,ext] = fileparts(listing(1).name);
            nseries = numel(listing);
            if nseries<100
                numzero = '%02g';
            elseif nseries<1000
                numzero = '%03g';
            elseif nseries<10000
                numzero = '%04g';
            elseif nseries<100000
                numzero = '%05g';
            end     
            if out_type==2 % .img/.hdr
                movefile(fullfile(guidata1.select_output_text.String,[fname,'.img']),...
                    fullfile(guidata1.select_output_text.String,...
                    [guidata1.output_name_text.String,'.img']))
                movefile(fullfile(guidata1.select_output_text.String,[fname,'.hdr']),...
                    fullfile(guidata1.select_output_text.String,...
                    [guidata1.output_name_text.String,'.hdr']))
            elseif out_type==0
                movefile(fullfile(guidata1.select_output_text.String,[fname,ext]),...
                    fullfile(guidata1.select_output_text.String,...
                    [guidata1.output_name_text.String,ext]))
            elseif out_type==1
                movefile(fullfile(guidata1.select_output_text.String,[fname,ext]),...
                    fullfile(guidata1.select_output_text.String,...
                    [guidata1.output_name_text.String,'.nii.gz']))
            elseif out_type==4
                for ix = 1:nseries
                    movefile(fullfile(guidata1.select_output_text.String,listing(ix).name),...
                        fullfile(guidata1.select_output_text.String,...
                        [guidata1.output_name_text.String,sprintf(numzero,ix),ext]))
                end
            elseif out_type==5
                for ix = 1:nseries
                    movefile(fullfile(guidata1.select_output_text.String,listing(ix).name),...
                        fullfile(guidata1.select_output_text.String,...
                        [guidata1.output_name_text.String,sprintf(numzero,ix),'.nii.gz']))
                end
            end  
        else
            movefile(fullfile(guidata1.select_output_text.String,'dcmHeaders.mat'),...
                fullfile(guidata1.select_output_text.String,...
                [field_names{1},'_dcmHeaders.mat']))
        end
    else
        errordlg('ERROR: One or more invalid directories. Respecify.')
    end
    guidata(main_figure,guidata1)
end

%% II. View/Edit 3D Image

% Simply call neuroimage_editor
function viewEdit_button_callback(hObject, ~) 
    % Retrieve Data
    guidata1 = guidata(main_figure);
    % Determine State
    guidata1.viewEdit_state = get(hObject,'Value');
    if isgraphics(guidata1.h); delete(guidata1.h); end
    if guidata1.viewEdit_state
        set(guidata1.dicom_button,'Value',false);
        set(guidata1.structurals_button,'Value',false);
        set(guidata1.functionals_button,'Value',false);
        if isgraphics(guidata1.dicom_figure)
            delete(guidata1.dicom_figure)
        end
        if isgraphics(guidata1.structurals_figure)
            delete(guidata1.structurals_figure)
        end
        if isgraphics(guidata1.funct_figure)
            delete(guidata1.funct_figure)
        end
    else
        if isgraphics(guidata1.viewer_figure)
            delete(guidata1.viewer_figure)
        end
        return;
    end
    guidata1.viewer_figure = nifti_studio;
    if ~isempty(guidata1.viewer_figure.figure)
        set(guidata1.viewer_figure.figure,'DeleteFcn',@closed_viewer_figure) 
    else
        set(guidata1.viewEdit_button,'Value',false)
    end
    % Save Data
    guidata(main_figure,guidata1)
end

% If Viewer/Editor is closed, Return Button Value to 'false'
function closed_viewer_figure(~,~)
    guidata1 = guidata(main_figure);
    guidata1.viewEdit_button.Value=false;
    guidata(main_figure,guidata1)   
end

%% III. STRUCTURALS

function structurals_button_callback(hObject, eventdata, overlay_num) %#ok
    % Retrieve Data
    guidata1 = guidata(main_figure);
    % Determine State
    guidata1.structurals_state = get(hObject,'Value');
    if isgraphics(guidata1.h); delete(guidata1.h); end
    if guidata1.structurals_state
        set(guidata1.viewEdit_button,'Value',false)
        set(guidata1.dicom_button,'Value',false);
        set(guidata1.functionals_button,'Value',false);
        if isgraphics(guidata1.viewer_figure)
            delete(guidata1.viewer_figure)
        end
        if isgraphics(guidata1.dicom_figure)
            delete(guidata1.dicom_figure)
        end
        if isgraphics(guidata1.funct_figure)
            delete(guidata1.funct_figure)
        end
    else
        if isgraphics(guidata1.structurals_figure)
            delete(guidata1.structurals_figure)
        end
        return;
    end
    % Figure: Structurals 
%     figure_pos1 = [.27*guidata1.screen_res(3), .4*guidata1.screen_res(4), ...
%         .46*guidata1.screen_res(3), .32*guidata1.screen_res(4)];
    guidata1.structurals_figure = figure('units','normalized','Position',...
        [.2695,.3991,.46,.32],'MenuBar','none','Name','Structural Quality',...
        'NumberTitle','off','Color',[0,0,0],'Visible','on','doublebuffer',...
        'on','CloseRequestFcn',@structurals_closereq);
    % Button: "Select Structural"
    guidata1.select_struct_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.01,.85,.2,.12],...
        'BackgroundColor', [0,0,0],'ForegroundColor',[1,1,1],...
        'FontSize',10,'FontWeight','bold','String',...
        'Select Structural','callback',@select_struct_callback,'Interruptible','off');  
    % Input Box: Select Structural
    guidata1.select_struct_text = uicontrol(guidata1.structurals_figure,'style',...
        'edit','Units','normalized','Position',[.22,.85,.7,.12],...
        'BackgroundColor', [1,1,1],'callback',@struct_input_box);
    if ~isempty(guidata1.struct_name_string)
        guidata1.select_struct_text.String = guidata1.struct_name_string;
    end
    % Button: "Load"
    guidata1.load_struct_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.93,.85,.06,.12],...
        'BackgroundColor', [0,0,0],'ForegroundColor',[1,1,1],...
        'FontSize',10,'FontWeight','bold','String',...
        'Load','callback',@load_struct_callback,'Interruptible','off'); 
    % ButtonGroup: "ROIs" Section
    guidata1.struct_rois_bg = uibuttongroup('Title','ROIs','Position',...
        [.01,.02,.48,.8],'BackgroundColor',[0,0,0],'ForegroundColor',guidata1.font_color,...
        'FontName','Helvetica','FontSize',11,'FontWeight','bold','TitlePosition',...
        'centertop');
    % Text: "Automatically specify..."
    uicontrol('Parent',guidata1.structurals_figure,'Style','text','String',...
        'Automatically specify ROIs based on...','Units','normalized',...
        'FontSize',10,'FontName','Helvetica','Position',[.018,.63,.39,.08],...
        'BackgroundColor',[0,0,0],'ForegroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','left');
    % Button: "Image Region"
    guidata1.struct_ROI_auto_region_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.025,.49,.21,.12],...
        'BackgroundColor',[0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Image Region','ForegroundColor',guidata1.font_color,'callback',{@struct_ROI_spec,1},'Interruptible','off');  
    % Button: "Signal Intensity"
    guidata1.struct_ROI_auto_intensity_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.265,.49,.21,.12],...
        'BackgroundColor',[0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Signal Intensity','ForegroundColor',guidata1.font_color,'callback',{@struct_ROI_spec,2},'Interruptible','off'); 
    % Text: "Manually specify..." 
    uicontrol('Parent',guidata1.structurals_figure,'Style','text','String',...
        'Manually specify ROIs based on...','Units','normalized',...
        'FontSize',10,'FontName','Helvetica','Position',[.018,.36,.35,.08],...
        'BackgroundColor',[0,0,0],'ForegroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','left');
    % Button: "Intensity Threshold"
    guidata1.struct_ROI_manual_intensity_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.025,.22,.21,.12],...
        'BackgroundColor',[0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Intensity Threshold','ForegroundColor',guidata1.font_color,'callback',{@struct_ROI_spec,3},'Interruptible','off');  
    % Button: "Manual Tracing"
    guidata1.struct_ROI_manual_trace_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.265,.22,.21,.12],...
        'BackgroundColor',[0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Manual Tracing','ForegroundColor',guidata1.font_color,'callback',{@struct_ROI_spec,4},'Interruptible','off'); 
    % Button: "View Structural"
    guidata1.struct_view_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.025,.06,.21,.12],...
        'BackgroundColor',guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'View Structural','ForegroundColor','k','callback',@struct_view_struct,'Interruptible','off');      
    % Button: "View ROIs"
    guidata1.struct_view_ROIs_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.265,.06,.21,.12],...
        'BackgroundColor',guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'View ROIs','ForegroundColor','k','callback',@struct_view_ROIs,'Interruptible','off');  
    % ButtonGroup: "Calculate" Section
    guidata1.struct_rois_bg = uibuttongroup('Title','Calculate','Position',...
        [.51,.02,.48,.8],'BackgroundColor',[0,0,0],'FontName',...
        'Helvetica','FontSize',11,'FontWeight','bold','TitlePosition',...
        'centertop','ForegroundColor',guidata1.font_color);
    % Button: "SNR (Noise SD)"
    guidata1.struct_SNR_sd_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.525,.58,.21,.12],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'SNR (Noise SD)','ForegroundColor',guidata1.font_color,'callback',{@struct_SNR_callback,1,1},'Interruptible','off');  
    % Button: "Display Image" 
    guidata1.struct_display_SNR_sd_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.765,.58,.21,.12],...
        'BackgroundColor',guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'Display Image','ForegroundColor','k','callback',{@struct_display_SNR,1},'Interruptible','off'); 
    % Button: "SNR (Noise Mean)"
    guidata1.struct_SNR_mean_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.525,.4,.21,.12],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'SNR (Noise Mean)','ForegroundColor',guidata1.font_color,'callback',{@struct_SNR_callback,2,1},'Interruptible','off');
    % Button: "Display Image" 
    guidata1.struct_display_SNR_mean_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.765,.4,.21,.12],...
        'BackgroundColor',guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'Display Image','ForegroundColor','k','callback',{@struct_display_SNR,2},'Interruptible','off'); 
    % Button: "Generate Report:" 
    guidata1.struct_report_button = uicontrol(guidata1.structurals_figure,'style',...
        'pushbutton','Units','normalized','Position',[.645,.24,.21,.12],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Generate Report:','ForegroundColor','w','callback',@struct_report,'Interruptible','off'); 
    box_x = linspace(.57,.91,5);
    % Checkbox 1: Stats.txt
    guidata1.struct_box1 = uicontrol(guidata1.structurals_figure,'style',...
        'checkbox','Units','normalized','Position',[box_x(1),.105,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',1); 
    % Checkbox 2: SNR img Mosaics
    guidata1.struct_box2 = uicontrol(guidata1.structurals_figure,'style',...
        'checkbox','Units','normalized','Position',[box_x(2),.105,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',1); 
    % Checkbox 3: ROIs img 
    guidata1.struct_box3 = uicontrol(guidata1.structurals_figure,'style',...
        'checkbox','Units','normalized','Position',[box_x(3),.105,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',0); 
    % Checkbox 4: SNR (SD) img
    guidata1.struct_box4 = uicontrol(guidata1.structurals_figure,'style',...
        'checkbox','Units','normalized','Position',[box_x(4),.105,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',0); 
    % Checkbox 5: SNR (Mean) img
    guidata1.struct_box5 = uicontrol(guidata1.structurals_figure,'style',...
        'checkbox','Units','normalized','Position',[box_x(5),.105,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',0);
    % Add Checkbox Labels:
    box_labels = {'Stats.txt','SNR Mosaics','ROIs img','SNR (SD) img','SNR (Mean) img'};
    box_pos = [.515,.05,.13,.06]; 
    box_x_pos = linspace(.52,.855,5); 
    box_y_pos = [.16,.035,.16,.035,.16];
    for ix = 1:5
        box_pos(1) = box_x_pos(ix);
        box_pos(2) = box_y_pos(ix);
        uicontrol('Parent',guidata1.structurals_figure,'Style','text','String',...
          box_labels{ix},'Units','normalized','FontSize',8,'FontName',...
          'Helvetica','Position',box_pos,'BackgroundColor',...
          [0,0,0],'ForegroundColor','w','FontWeight','Normal',...
          'HorizontalAlignment','center');
    end
    % Save Data
    guidata(main_figure,guidata1)
end

function select_struct_callback(~,~,~) 
    guidata1 = guidata(main_figure);
    cd(guidata1.last_struct_dir)
    [guidata1.struct_name, new_path] = uigetfile(guidata1.extensions,...
        'Select Structural Image:','MultiSelect','off');
    cd(guidata1.first_dir)
    if guidata1.struct_name~=0
        guidata1.last_struct_dir = new_path;
        [~,~,ext] = fileparts(guidata1.struct_name);
        for ix = 1:3
            if ~isempty(strfind(guidata1.extensions{ix},ext))
                type = ix;
                others = setdiff(1:3,ix);
                break;
            end
        end
        guidata1.extensions = guidata1.extensions([type,others]);
        if guidata1.last_struct_dir~=0
            guidata1.select_struct_text.String = fullfile(guidata1.last_struct_dir,guidata1.struct_name);
            guidata1.struct_name_string = guidata1.select_struct_text.String;
        end
    end
    guidata(main_figure,guidata1)
end

function struct_input_box(~, ~, ~)
    guidata1 = guidata(main_figure);
    if ~strcmp(guidata1.select_struct_text.String,'')
        guidata1.struct_name_string = guidata1.select_struct_text.String;
        guidata(main_figure,guidata1)
    end
end

function load_struct_callback(~, ~, ~)
    guidata1 = guidata(main_figure);
    if ~isempty(guidata1.struct_name_string)
        if exist(guidata1.struct_name_string,'file')
            try
                guidata1.struct = load_nii(guidata1.struct_name_string);
                disp('Structural image successfully loaded.')
                guidata1.apply_header = 1;
            catch
                guidata1.struct = load_untouch_nii(guidata1.struct_name_string);
                disp('Non-orthogonal shearing detected in affine matrix. Successfully loaded raw image data.') 
                guidata1.apply_header = 0;
            end
            guidata1.struct.img = single(guidata1.struct.img); % change data type for later functions
        else
            errordlg('Structural image not found. Please respecify.'); return;
        end
    else
        errordlg('Please specify a structural image.'); return;
    end
    guidata(main_figure,guidata1)
end

function struct_ROI_spec(~, ~, ROI_type) 
    guidata1 = guidata(main_figure);
    if ~isempty(guidata1.struct)
        switch ROI_type
            case 1 %'Image Region' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                guidata1.struct_ROI_type = 'Automatic - Image Region';
                disp('Automatically calculating ROI locations based on standard image regions.')
                % Initialize struct_ROI:
                guidata1.struct_dim = guidata1.struct.hdr.dime.dim(2:4);
                guidata1.struct_ROI = guidata1.struct;
                guidata1.struct_ROI.img = zeros(guidata1.struct_dim);
                % Signal ROI:
                center_x = round(guidata1.struct_dim(1)/2);
                center_y = round(guidata1.struct_dim(2)/2);
                center_z = round(guidata1.struct_dim(3)/2);
                central_x = round(.1*guidata1.struct_dim(1));
                central_y = round(.1*guidata1.struct_dim(2));
                central_z = round(.15*guidata1.struct_dim(3));
                guidata1.struct_ROI.img(center_x-central_x:center_x+central_x,...
                    center_y-central_y:center_y+central_y,center_z-central_z:center_z+central_z) = 1;
                % Noise ROIs
                border_x = round(.08*guidata1.struct_dim(1));
                border_y = round(.75*guidata1.struct_dim(2));
                guidata1.struct_ROI.img(3:border_x,3:border_y,2:guidata1.struct_dim(3)-1) = 2;
                guidata1.struct_ROI.img(guidata1.struct_dim(1)-border_x:guidata1.struct_dim(1)-2,...
                    3:border_y,2:guidata1.struct_dim(3)-1) = 2;
            case 2 %'Signal Intensity' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                guidata1.struct_ROI_type = 'Automatic - Signal Intensity';
%                 disp('Automatically calculting ROI locations based on signal intensities.')
                % Initialize struct_ROI:
                guidata1.struct_dim = guidata1.struct.hdr.dime.dim(2:4);
                guidata1.struct_ROI = guidata1.struct;
                guidata1.struct_ROI.img = zeros(guidata1.struct_dim);
                % Determine Thresholds:
                intensities = sort(guidata1.struct.img(:));
                sig_min_perc = .75; sig_max_perc = .995;
                sig_min = intensities(round(sig_min_perc*length(intensities)));
                sig_max = intensities(round(sig_max_perc*length(intensities)));
                noise_perc = .7; 
                noise_thresh = intensities(round(noise_perc*length(intensities)));
                % Signal ROI:
                guidata1.struct_ROI.img(((guidata1.struct.img>=sig_min) + (guidata1.struct.img<=sig_max))==2) = 1;
                % Noise ROI:
                guidata1.struct_ROI.img(guidata1.struct.img<=noise_thresh) = 2;
                fprintf(['ROIs:  Intensity thresholds for signal ROI were ',...
                    '%.02f - %.02f. Maximum intensity for background ROI was %.02f.'],...
                    [sig_min,sig_max,noise_thresh])
            case 3 %'Intensity Threshold' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                guidata1.struct_ROI_type = 'Manual - Intensity Threshold';
                % Initialize struct_ROI:
                guidata1.struct_dim = guidata1.struct.hdr.dime.dim(2:4);
                guidata1.struct_ROI = guidata1.struct;
                guidata1.struct_ROI.img = zeros(guidata1.struct_dim);
                max_sig = max(guidata1.struct.img(:));
                % Specify thresholds for inclusive mask:
                prompt = {sprintf('Enter intensity thresholds for SIGNAL ROI: \n\nMin:'),'Max:'};
                dlg_title = 'Signal ROI'; num_lines = [1,50;1,50]; defaultans = {'0',num2str(max_sig)};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                if isempty(answer)
                    disp('User cancelled intensity threshold specification');
                    return;
                end
                sig_min = str2double(answer(1)); sig_max = str2double(answer(2));
                prompt = {'Enter maximum intensity for BACKGROUND ROI: '};
                dlg_title = 'Background ROI'; num_lines = [1,50;]; defaultans = {num2str(max_sig)};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                if isempty(answer)
                    disp('User cancelled intensity threshold specification');
                    return;
                end
                noise_thresh = str2double(answer(1));
                % Signal ROI:
                guidata1.struct_ROI.img(((guidata1.struct.img>=sig_min) + (guidata1.struct.img<=sig_max))==2) = 1;
                % Noise ROI:
                guidata1.struct_ROI.img(guidata1.struct.img<=noise_thresh) = 2;
                disp('Calculting ROI locations based on specified signal intensities.')
            case 4 %'Manual Tracing' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                guidata1.struct_ROI_type = 'Manual - Tracing';
%                 if isgraphics(guidata1.h); close(guidata1.h); end
                instructions_text = sprintf(['First, draw a SIGNAL ROI (green). ',...
                    'Second, draw a BACKGROUND/NOISE ROI (red). ',...
                    'To accomplish this, select "Draw" menu -> "Draw Color". ',...
                    'Once a color is selected, simply trace the outline of a shape on the image, ',...
                    'and this shape will automatically fill with color. Alternatively, select ',...
                    '"Draw" menu -> "Shapes" to draw using a specified shape. ',...
                    'Click and drag on the image to create shapes. To navigate between slices, ',...
                    'use the up and down arrow keys. An ROI drawn on one slice can be propagated through to other ',...
                    'slices by selecting "Draw" menu -> "Edit Drawing" -> "Propagate Through Slices" -> choose slice range. ',...
                    'Note that this will only apply to the most recent drawing. When you are finished drawing ',...
                    'select "Confirm ROIs" -> "Confirm". The ROIs will be automatically saved.']);
                h_mbox = msgbox(instructions_text,'ROI Tracing Instructions');
                uiwait(h_mbox,60)
                [guidata1.h,~,drawing_idx] = neuroimage_editor_tbx('apply_header',...
                    guidata1.apply_header,'background',guidata1.struct,...
                    'colorbar',1,'title_on',1,'axis_tick_on',1); %#ok
                uiwait(guidata1.h.figure,120) % wait until neuroimage_editor closes, or 60 sec
                drawing_idx = evalin('base','drawing_idx'); % pull last_draw in from 'base' workspace
                drawing_idx = permute(drawing_idx,[2,1,3]);
                guidata1.struct_ROI = guidata1.struct;
                guidata1.struct_ROI.img = zeros(guidata1.struct_dim);
                guidata1.struct_ROI.img(drawing_idx(:)==2) = 1;
                guidata1.struct_ROI.img(drawing_idx(:)==3) = 2;
        end
        if ROI_type~=4
            guidata1.struct_mean_signal = nanmean(guidata1.struct.img(guidata1.struct_ROI.img(:)==1));
        end
    else
        errordlg('Please load a structural image.'); return;
    end
    guidata(main_figure,guidata1)        
end

function struct_view_struct(hObject, eventdata, overlay_num) %#ok
    guidata1 = guidata(main_figure);
    if isgraphics(guidata1.h); close(guidata1.h); end
    if ~isempty(guidata1.struct) 
        handles = nifti_studio('apply_header',guidata1.apply_header,...
            'background',guidata1.struct,'colormap','gray',...
            'colorbar_on',1,'title_on',1,'axis_tick_on',1);
        guidata1.h = handles.figure;
        guidata(main_figure,guidata1) 
    else
        errordlg('Please load a structural image.'); return;
    end
end

function struct_view_ROIs(hObject, eventdata, overlay_num) %#ok
    guidata1 = guidata(main_figure);
    if isgraphics(guidata1.h); close(guidata1.h); end
    if ~isempty(guidata1.struct) 
        if ~isempty(guidata1.struct_ROI)
            handles = nifti_studio('apply_header',guidata1.apply_header,...
                'background',guidata1.struct,'colormap','gray',...
                'overlay',guidata1.struct_ROI,'colorbar_on',0,'title_on',0,'axis_tick_on',0);
            guidata1.h = handles.figure;
            guidata(main_figure,guidata1) 
        else
            errordlg('Please Specify Signal & Noise ROIs.'); return;
        end
    else
        errordlg('Please load a structural image.'); return;
    end
end

function struct_SNR_callback(~, ~, SNR_type, verbose) 
    guidata1 = guidata(main_figure);
    % Retrieve ROIs:
    if ~isempty(guidata1.struct_ROI) 
%         if isgraphics(guidata1.h); close(guidata1.h); end
    else
        errordlg('ERROR: No ROIs found.'); return;
    end
    % Calculate SNR:
    switch SNR_type
        case 1 % SNR (SD)
            % Calculate Stat:
            noise_SD = nanstd(guidata1.struct.img(guidata1.struct_ROI.img(:)==2));
            denominator = noise_SD; % sqrt(2/(4-pi))*noise_SD;
            if noise_SD == 0
                warning('The signal standard deviation within the noise ROI was 0. Has this image been thresholded or filtered?')
            end
            guidata1.struct_SNR_sd = guidata1.struct_mean_signal / denominator; % scalar
            if verbose
                disp(['Structural SNR (SD) = ',num2str(guidata1.struct_SNR_sd)])
            end
            % Calculate Local:
            guidata1.struct_SNR_sd_local = guidata1.struct;
            guidata1.struct_SNR_sd_local.img(guidata1.struct_ROI.img(:)==1) = guidata1.struct_SNR_sd_local.img(guidata1.struct_ROI.img(:)==1)./denominator;
            guidata1.struct_SNR_sd_local.img(guidata1.struct_ROI.img(:)~=1) = 0;
            % Calculate Global:
            guidata1.struct_SNR_sd_global = guidata1.struct;
            guidata1.struct_SNR_sd_global.img = guidata1.struct_SNR_sd_global.img./denominator;
        case 2 % SNR (MEAN)
            % Calculate Stat:
            noise_mean = nanmean(guidata1.struct.img(guidata1.struct_ROI.img(:)==2));
            denominator = noise_mean; % sqrt(2/pi)*noise_mean;
            if noise_mean == 0
                warning('The mean signal within the noise ROI was 0. Has this image been thresholded or filtered?')
            end
            guidata1.struct_SNR_mean = guidata1.struct_mean_signal / denominator;
            if verbose
                disp(['Structural SNR (Mean) = ',num2str(guidata1.struct_SNR_mean)])
            end
            % Calculate Local:
            guidata1.struct_SNR_mean_local = guidata1.struct;
            guidata1.struct_SNR_mean_local.img(guidata1.struct_ROI.img(:)==1) = guidata1.struct_SNR_mean_local.img(guidata1.struct_ROI.img(:)==1)./denominator;
            guidata1.struct_SNR_mean_local.img(guidata1.struct_ROI.img(:)~=1) = 0;
            % Calculate Global:
            guidata1.struct_SNR_mean_global = guidata1.struct;
            guidata1.struct_SNR_mean_global.img = guidata1.struct_SNR_mean_global.img./denominator;
    end
    guidata(main_figure,guidata1) 
end

function struct_display_SNR(hObject, eventdata, SNR_type) %#ok
    guidata1 = guidata(main_figure);
%     if isgraphics(guidata1.h); close(guidata1.h); end
    switch SNR_type
        case 1 % SNR (SD)
            if isempty(guidata1.struct_SNR_sd_local) || isempty(guidata1.struct_SNR_sd_global)
                errordlg('ERROR: First, calculate SNR (Noise SD)'); return;
            end
             % Check whether to display signal ROI only or Full Image:
            ans1 = 'Signal ROI'; ans2 = 'Full Image';
            answer = questdlg('Display which?','Display',ans1,ans2,ans2);
            switch answer
                case ans1
                    handles = nifti_studio('apply_header',guidata1.apply_header,...
                        'background',guidata1.struct_SNR_sd_local,...
                        'colorbar_on',1,'title_on',1,'axis_tick_on',1,...
                        'colormap','jet','title','SNR (Noise SD)',...
                        'background_caxis',[min(guidata1.struct_SNR_sd_local.img(:)),max(guidata1.struct_SNR_sd_local.img(:))]);
                    guidata1.h = handles.figure;
                case ans2
                    handles = nifti_studio('apply_header',guidata1.apply_header,...
                        'background',guidata1.struct,...
                        'overlay',guidata1.struct_SNR_sd_global,...
                        'colorbar_on',1,'title_on',1,'axis_tick_on',1,...
                        'colormap','jet','title','SNR (Noise SD)',...
                        'overlay_caxis',[0,max(guidata1.struct_SNR_sd_local.img(:))]);
                    guidata1.h = handles.figure;
            end
        case 2 % SNR (MEAN)
            if isempty(guidata1.struct_SNR_mean_local) || isempty(guidata1.struct_SNR_mean_global)
                errordlg('ERROR: First, calculate SNR (Noise Mean)'); return;
            end
            % Check whether to display signal ROI only or Full Image:
            ans1 = 'Signal ROI'; ans2 = 'Full Image';
            answer = questdlg('Display which?','Display',ans1,ans2,ans2);
            switch answer
                case ans1
                    handles = nifti_studio('apply_header',guidata1.apply_header,...
                        'background',guidata1.struct_SNR_mean_local,...
                        'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet','title','SNR (Noise Mean)',...
                        'background_caxis',[min(guidata1.struct_SNR_mean_local.img(:)),max(guidata1.struct_SNR_mean_local.img(:))]);
                    guidata1.h = handles.figure;
                case ans2
                    handles = nifti_studio('apply_header',guidata1.apply_header,...
                        'background',guidata1.struct,...
                        'overlay',guidata1.struct_SNR_mean_global,...
                        'colorbar_on',1,'title_on',1,'axis_tick_on',1,...
                        'colormap','jet','title','SNR (Noise Mean)',...
                        'overlay_caxis',[0,max(guidata1.struct_SNR_mean_local.img(:))]);
                    guidata1.h = handles.figure;
            end
    end
    guidata(main_figure,guidata1) 
end

function struct_report(~, ~, ~) 
    guidata1 = guidata(main_figure);
    if isgraphics(guidata1.h); close(guidata1.h); end
    if isempty(guidata1.struct)
        errordlg('Structural image not loaded.','ERROR'); return;
    end
    [guidata1.report_dir] = uigetdir(guidata1.last_struct_dir,'Select Output Folder:');
    if guidata1.report_dir == 0
       disp('Generate report canceled.') 
       return
    end
    % Get Checkbox Options:
    guidata1.struct_report_options = [guidata1.struct_box1.Value,...
        guidata1.struct_box2.Value,guidata1.struct_box3.Value,...
        guidata1.struct_box4.Value,guidata1.struct_box5.Value];
    guidata(main_figure,guidata1) 
    if isempty(guidata1.struct_ROI)
        disp('WARNING: No ROIs specified. Automatically calculating based on Signal Intensity.')
        struct_ROI_spec([], [], 2)
    end
    % Recalculate these in case anything has changed:
    struct_SNR_callback([], [], 1, 0);
    struct_SNR_callback([], [], 2, 0);
    % Outputs:
    guidata1 = guidata(main_figure);
    if guidata1.struct_report_options(1) % Stats.txt
        try [~,fname1,~] = fileparts(guidata1.struct.fileprefix);
        catch; fname1 = 'Structural';
        end
        fname2 = ['Stats_',fname1,'.txt'];
        % Write Text File:
        fileID = fopen(fullfile(guidata1.report_dir,fname2),'w');
        fprintf(fileID,'%s\r\n',['Structural Name:  ',fname1]);
        fprintf(fileID,'%s\r\n',['ROI Type:  ',guidata1.struct_ROI_type]);
        fprintf(fileID,'%s\r\n',['Signal ROI # Voxels:  ',num2str(sum(guidata1.struct_ROI.img(:)==1))]);
        fprintf(fileID,'%s\r\n',['Background ROI # Voxels:  ',num2str(sum(guidata1.struct_ROI.img(:)==2))]);
        fprintf(fileID,'%s\r\n',['SNR (SD) = ',num2str(guidata1.struct_SNR_sd)]);
        fprintf(fileID,'%s\r\n',['SNR (Mean) = ',num2str(guidata1.struct_SNR_mean)]);
        fclose(fileID);
    end
    % Determine slices and dimension for mosaic:
    if guidata1.struct_report_options(2) 
        dim3 = size(guidata1.struct.img,3);
        if (dim3>=15) 
            slices = round(linspace(3,.9*dim3,15));
            axes_dim = [3,5];
        elseif (dim3>=10) 
            slices = round(linspace(3,.9*dim3,9));
            axes_dim = [2,5];
        elseif (dim3>=5) 
            slices = round(linspace(3,.9*dim3,5));
            axes_dim = [1,5];
        else
            slices = dim3;
            axes_dim = [1,slices+1];
        end
    end
    if guidata1.struct_report_options(2) % SNR Mosaics
%         handles = nifti_studio('apply_header',guidata1.apply_header,...
%             'background',guidata1.struct,...
%             'overlay',guidata1.struct_SNR_sd_global,...
%             'colorbar_on',1,'title_on',1,'axis_tick_on',1,...
%             'colormap','gray','title','SNR (Noise SD)',...
%             'overlay_caxis',[0,max(guidata1.struct_SNR_sd_local.img(:))],...
%             'print',fullfile(guidata1.report_dir,'struct_SNR_sd_mosaic.png'));
        [handles] = neuroimage_editor_mosaic('background',guidata1.struct,...
            'overlay',guidata1.struct_SNR_sd_global,'slices',slices,'colormap','gray',...
            'slice_locator',0,'colorbar',1,'axes_dim',axes_dim,'figure_pos',[.1593,.03796,.71563,.908],...
            'print',fullfile(guidata1.report_dir,'struct_SNR_sd_mosaic.png'));
        close(handles.figure);
%         handles = nifti_studio('apply_header',guidata1.apply_header,...
%             'background',guidata1.struct,...
%             'overlay',guidata1.struct_SNR_mean_global,...
%             'colorbar_on',1,'title_on',1,'axis_tick_on',1,...
%             'colormap','gray','title','SNR (Noise Mean)',...
%             'overlay_caxis',[0,max(guidata1.struct_SNR_mean_local.img(:))],...
%             'print',fullfile(guidata1.report_dir,'struct_SNR_mean_mosaic.png'));
        [handles] = neuroimage_editor_mosaic('background',guidata1.struct,...
            'overlay',guidata1.struct_SNR_mean_global,'slices',slices,'colormap','gray',...
            'slice_locator',0,'colorbar',1,'axes_dim',axes_dim,'figure_pos',[.1593,.03796,.71563,.908],...
            'print',fullfile(guidata1.report_dir,'struct_SNR_mean_mosaic.png'));
            close(handles.figure);
    end
    if guidata1.struct_report_options(3) % ROIs Image
        if guidata1.apply_header
            try save_nii(guidata1.struct_ROI,fullfile(guidata1.report_dir,'struct_ROIs'))
            catch
                save_untouch_nii(guidata1.struct_ROI,fullfile(guidata1.report_dir,'struct_ROIs'))
            end
        else
            save_untouch_nii(guidata1.struct_ROI,fullfile(guidata1.report_dir,'struct_ROIs'))
        end
    end
    if guidata1.struct_report_options(4) % SNR (SD) img
        if guidata1.apply_header
            try save_nii(guidata1.struct_SNR_sd_global,fullfile(guidata1.report_dir,'struct_SNR_sd'))
            catch
                save_untouch_nii(guidata1.struct_SNR_sd_global,fullfile(guidata1.report_dir,'struct_SNR_sd'))
            end
        else
            save_untouch_nii(guidata1.struct_SNR_sd_global,fullfile(guidata1.report_dir,'struct_SNR_sd'))
        end
    end
    if guidata1.struct_report_options(5) % SNR (Mean) img
        if guidata1.apply_header
            try save_nii(guidata1.struct_SNR_mean_global,fullfile(guidata1.report_dir,'struct_SNR_mean'))
            catch
                save_untouch_nii(guidata1.struct_SNR_mean_global,fullfile(guidata1.report_dir,'struct_SNR_mean'))
            end
        else
            save_untouch_nii(guidata1.struct_SNR_mean_global,fullfile(guidata1.report_dir,'struct_SNR_mean'))
        end
    end
    guidata(main_figure,guidata1) 
    disp('Report completed.')
end

function structurals_closereq(~, ~, ~)
    guidata1 = guidata(main_figure);
    guidata1.structurals_button.Value=0;
    if isgraphics(guidata1.structurals_figure)
        delete(guidata1.structurals_figure)
    end
    if isgraphics(guidata1.h); close(guidata1.h); end
    guidata(main_figure,guidata1)
end

%% IV. Functionals

function functionals_button_callback(hObject, ~, ~)
    % Retrieve Data
    funct_initialize % initialize data or clear old data
    guidata1 = guidata(main_figure);
    % Determine State
    guidata1.functionals_state = get(hObject,'Value');
    if isgraphics(guidata1.h); close(guidata1.h); end
    if guidata1.functionals_state
        set(guidata1.viewEdit_button,'Value',false);
        set(guidata1.dicom_button,'Value',false);
        set(guidata1.structurals_button,'Value',false);
        if isgraphics(guidata1.viewer_figure)
            delete(guidata1.viewer_figure)
        end
        if isgraphics(guidata1.dicom_figure)
            delete(guidata1.dicom_figure)
        end
        if isgraphics(guidata1.structurals_figure)
            delete(guidata1.structurals_figure)
        end
    else
        if isgraphics(guidata1.funct_figure)
            delete(guidata1.funct_figure)
        end
        return;
    end
    % Figure: functionals 
%     figure_pos1 = [.28*guidata1.screen_res(3), .266*guidata1.screen_res(4),...
%         .46*guidata1.screen_res(3), .525*guidata1.screen_res(4)];
    guidata1.funct_figure = figure('units','normalized','Position',...
        [.2795,.2651,.46,.525],'MenuBar','none','Name','Functional Quality',...
        'NumberTitle','off','Color',[0,0,0],'Visible','on','doublebuffer',...
        'on','CloseRequestFcn',@functionals_closereq);
    % Button: "Select functional"
    guidata1.select_funct_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.01,.9,.2,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Select Functional','callback',@select_funct_callback,...
        'Interruptible','off','ForegroundColor','w');  
    % Input Box: Select functional
    guidata1.select_funct_text = uicontrol(guidata1.funct_figure,'style',...
        'edit','Units','normalized','Position',[.22,.9,.7,.075],...
        'BackgroundColor', [1,1,1],'callback',@funct_input_box);
    if ~isempty(guidata1.funct_name_string)
        guidata1.select_funct_text.String = guidata1.funct_name_string;
    end
    % Button: "Load"
    guidata1.load_funct_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.93,.9,.06,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Load','callback',@load_funct_callback,'Interruptible','off',...
        'ForegroundColor','w'); 
    % ButtonGroup: "ROIs" Section
    guidata1.funct_rois_bg = uibuttongroup('Title','ROIs','Position',...
        [.01,.285,.48,.6],'BackgroundColor',[0,0,0],'FontName',...
        'Helvetica','FontSize',11,'FontWeight','bold','TitlePosition',...
        'centertop','ForegroundColor','w');
    % Text: "Automatically specify..."
    uicontrol('Parent',guidata1.funct_figure,'Style','text','String',...
        'Automatically specify ROIs based on...','Units','normalized',...
        'FontSize',10,'FontName','Helvetica','Position',[.018,.75,.39,.08],...
        'BackgroundColor',[0,0,0],'ForegroundColor','w','FontWeight','Bold',...
        'HorizontalAlignment','left');
    % Button: "Image Region"
    guidata1.funct_ROI_auto_region_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.025,.7,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Image Region','ForegroundColor',guidata1.font_color,'callback',{@funct_ROI_spec,1},...
        'Interruptible','off');  
    % Button: "Signal Intensity"
    guidata1.funct_ROI_auto_intensity_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.265,.7,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Signal Intensity','ForegroundColor',guidata1.font_color,'callback',{@funct_ROI_spec,2},...
        'Interruptible','off'); 
    % Text: "Manually specify..." 
    uicontrol('Parent',guidata1.funct_figure,'Style','text','String',...
        'Manually specify ROIs based on...','Units','normalized',...
        'FontSize',10,'FontName','Helvetica','Position',[.018,.61,.35,.075],...
        'BackgroundColor',[0,0,0],'ForegroundColor','w','FontWeight','Bold',...
        'HorizontalAlignment','left');
    % Button: "Intensity Threshold"
    guidata1.funct_ROI_manual_intensity_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.025,.55,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Intensity Threshold','ForegroundColor',guidata1.font_color,'callback',{@funct_ROI_spec,3},...
        'Interruptible','off');  
    % Button: "Manual Tracing"
    guidata1.funct_ROI_manual_trace_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.265,.55,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Manual Tracing','ForegroundColor',guidata1.font_color,'callback',{@funct_ROI_spec,4},...
        'Interruptible','off'); 
    % Text: "Load brain mask..." 
    uicontrol('Parent',guidata1.funct_figure,'Style','text','String',...
        'Load brain mask...','Units','normalized',...
        'FontSize',10,'FontName','Helvetica','Position',[.018,.475,.35,.06],...
        'BackgroundColor',[0,0,0],'ForegroundColor','w','FontWeight','Bold',...
        'HorizontalAlignment','left');
    % Button: "Inclusive Mask"
    guidata1.funct_ROI_inclusive_mask = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.025,.41,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Inclusive Mask','ForegroundColor',guidata1.font_color,'callback',{@funct_ROI_spec,5},...
        'Interruptible','off');  
    % Button: "Exclusive Mask"
    guidata1.funct_ROI_exclusive_mask = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.265,.41,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Exclusive Mask','ForegroundColor',guidata1.font_color,'callback',{@funct_ROI_spec,6},...
        'Interruptible','off');  
    % Button: "View Functional"
    guidata1.funct_view_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.025,.31,.21,.075],...
        'BackgroundColor', guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'View Functional','ForegroundColor','k','callback',@funct_view_funct,...
        'Interruptible','off');      
    % Button: "View ROIs"
    guidata1.funct_view_ROIs_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.265,.31,.21,.075],...
        'BackgroundColor', guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'View ROIs','ForegroundColor','k','callback',@funct_view_ROIs,...
        'Interruptible','off');  
    % ButtonGroup: "Preprocess" Section
    guidata1.funct_filter_bg = uibuttongroup('Title','Preprocess','Position',...
        [.01,.02,.48,.26],'BackgroundColor',[0,0,0],'FontName',...
        'Helvetica','FontSize',11,'FontWeight','bold','TitlePosition',...
        'centertop','ForegroundColor','w');
    % Button: "Detrend"
    guidata1.funct_detrend_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.025,.145,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Detrend','ForegroundColor',guidata1.font_color,'callback',@funct_detrend,...
        'Interruptible','off');
    % Button: "Band-Pass"
    guidata1.funct_bandpass_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.025,.045,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Band-Pass','ForegroundColor',guidata1.font_color,'callback',@funct_bandpass,...
        'Interruptible','off');
    % Button: "View Signals"
    guidata1.funct_signals = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.265,.145,.21,.075],...
        'BackgroundColor', guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'View Signals','ForegroundColor','k','callback',@funct_view_signals,...
        'Interruptible','off');  
    % Button: "Power Spectrum"
    guidata1.funct_power_spec = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.265,.045,.21,.075],...
        'BackgroundColor',guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'Power Spectrum','ForegroundColor','k','callback',@funct_view_spectrum,...
        'Interruptible','off');  
    % ButtonGroup: "Calculate" Section
    guidata1.funct_rois_bg = uibuttongroup('Title','Calculate','Position',...
        [.505,.02,.485,.865],'BackgroundColor',[0,0,0],'FontName',...
        'Helvetica','FontSize',11,'FontWeight','bold','TitlePosition',...
        'centertop','ForegroundColor','w');
    % Button: "tSNR"
    guidata1.funct_tSNR_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.525,.747,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'tSNR','ForegroundColor',guidata1.font_color,'callback',{@funct_SNR_callback,1,1},...
        'Interruptible','off');  
    % Button: "Display tSNR" 
    guidata1.funct_display_tSNR_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.765,.747,.21,.075],...
        'BackgroundColor', guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'Display tSNR','ForegroundColor','k','callback',{@funct_display_SNR,1},...
        'Interruptible','off'); 
    % Button: "tSBNR"
    guidata1.funct_tSBNR_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.525,.6402,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'tSBNR','ForegroundColor',guidata1.font_color,'callback',{@funct_SNR_callback,2,1},...
        'Interruptible','off');  
    % Button: "Display tSBNR" 
    guidata1.funct_display_tSBNR_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.765,.6402,.21,.075],...
        'BackgroundColor', guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'Display tSBNR','ForegroundColor','k','callback',{@funct_display_SNR,2},...
        'Interruptible','off'); 
    % Button: "SFNR"
    guidata1.funct_SFNR_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.525,.5334,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'SFNR','ForegroundColor',guidata1.font_color,'callback',{@funct_SNR_callback,3,1},...
        'Interruptible','off');
    % Button: "Display SFNR" 
    guidata1.funct_display_SFNR_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.765,.5334,.21,.075],...
        'BackgroundColor',guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'Display SFNR','ForegroundColor','k','callback',{@funct_display_SNR,3},...
        'Interruptible','off'); 
    % Button: "SNR-Funct"
    guidata1.funct_SNR_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.525,.4266,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'SNR-Funct','ForegroundColor',guidata1.font_color,'callback',{@funct_SNR_callback,4,1},...
        'Interruptible','off');
    % Button: "Display SNR-Funct" 
    guidata1.funct_display_SNR_mean_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.765,.4266,.21,.075],...
        'BackgroundColor', guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'Display SNR-Funct','ForegroundColor','k','callback',{@funct_display_SNR,4},...
        'Interruptible','off'); 
    % Button: "Mean Power"
    guidata1.funct_mean_power_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.525,.32,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Mean Power','ForegroundColor',guidata1.font_color,'callback',{@funct_SNR_callback,5,1},...
        'Interruptible','off');
    % Button: "Display Power" 
    guidata1.funct_display_mean_power_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.765,.32,.21,.075],...
        'BackgroundColor',guidata1.font_color,'FontSize',10,'FontWeight','bold','String',...
        'Display Power','ForegroundColor','k','callback',{@funct_display_SNR,5},...
        'Interruptible','off'); 
    % Button: "Generate Report:" 
    guidata1.funct_report_button = uicontrol(guidata1.funct_figure,'style',...
        'pushbutton','Units','normalized','Position',[.645,.205,.21,.075],...
        'BackgroundColor', [0,0,0],'FontSize',10,'FontWeight','bold','String',...
        'Generate Report:','ForegroundColor','w','callback',@funct_report,...
        'Interruptible','off'); 
%     orig = [.5444,.5989,.6533,.7078,.7622,.8167,.8711,.9256];
%     box_xpos_orig = linspace(.5444,.9256,8);
    box_xpos = linspace(.535,.935,9);
    % Checkbox 1: Stats.txt
    guidata1.funct_box1 = uicontrol(guidata1.funct_figure,'style',...
        'checkbox','Units','normalized','Position',[box_xpos(1),.08,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',1); 
    % Checkbox 2: Img Mosaics
    guidata1.funct_box2 = uicontrol(guidata1.funct_figure,'style',...
        'checkbox','Units','normalized','Position',[box_xpos(2),.08,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',1); 
    % Checkbox 3: Preprocess Data
    guidata1.funct_box3 = uicontrol(guidata1.funct_figure,'style',...
        'checkbox','Units','normalized','Position',[box_xpos(3),.08,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',0);   
    % Checkbox 4: ROIs Img
    guidata1.funct_box4 = uicontrol(guidata1.funct_figure,'style',...
        'checkbox','Units','normalized','Position',[box_xpos(4),.08,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',0); 
    % Checkbox 5: tSNR Img
    guidata1.funct_box5 = uicontrol(guidata1.funct_figure,'style',...
        'checkbox','Units','normalized','Position',[box_xpos(5),.08,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',0); 
    % Checkbox 6: tSBNR Img
    guidata1.funct_box6 = uicontrol(guidata1.funct_figure,'style',...
        'checkbox','Units','normalized','Position',[box_xpos(6),.08,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',0); 
    % Checkbox 7: SFNR Img
    guidata1.funct_box7 = uicontrol(guidata1.funct_figure,'style',...
        'checkbox','Units','normalized','Position',[box_xpos(7),.08,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',0); 
    % Checkbox 8: SNR Img
    guidata1.funct_box8 = uicontrol(guidata1.funct_figure,'style',...
        'checkbox','Units','normalized','Position',[box_xpos(8),.08,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',0); 
    % Checkbox 9: Power Img
    guidata1.funct_box9 = uicontrol(guidata1.funct_figure,'style',...
        'checkbox','Units','normalized','Position',[box_xpos(9),.08,.05,.05],...
        'BackgroundColor',[0,0,0],'FontSize',8,...
        'FontWeight','Normal','ForegroundColor','w','Value',0); 
    % Add Checkbox Labels:
    box_labels = {'Stats.txt','Mosaics','Preprocessed','ROIs Img','tSNR Img','tSBNR Img','SFNR Img','SNR Img','Power Img'};
    box_pos = [.515,.05,.095,.04]; 
    box_x_pos = linspace(.512,.9,9);
    box_y_pos = [.12,.03,.12,.03,.12,.03,.12,.03,.12];
    box_widths = [.07,.095,.115,repmat(.095,1,5),.079];
    for ix = 1:9
        if ix==2 
            box_pos(1) = box_x_pos(ix)-.01;
        elseif ix==3
            box_pos(1) = box_x_pos(ix)-.018;
        elseif ix==4
            box_pos(1) = box_x_pos(ix)-.008;
        elseif ix==9
            box_pos(1) = box_x_pos(ix)+.008;
        else
            box_pos(1) = box_x_pos(ix);
        end
        box_pos(2) = box_y_pos(ix);
        box_pos(3) = box_widths(ix);
        uicontrol('Parent',guidata1.funct_figure,'Style','text','String',...
          box_labels{ix},'Units','normalized','FontSize',7,'FontName',...
          'Helvetica','Position',box_pos,'BackgroundColor',...
          [0,0,0],'FontWeight','Normal','ForegroundColor','w',...
          'HorizontalAlignment','center');
    end
    % Save Data
    guidata(main_figure,guidata1)
end

function select_funct_callback(~,~,~)
    guidata1 = guidata(main_figure);
    cd(guidata1.last_funct_dir)
    [funct_name, funct_path] = uigetfile(guidata1.extensions,...
        'Select One 4D Functional or Series of 3D Functionals:','MultiSelect','on');
    cd(guidata1.first_dir)
    if funct_path==0 % cancel
        return; 
    else
        funct_initialize
        guidata1 = guidata(main_figure);
        guidata1.funct_name = funct_name;
        guidata1.funct_path = funct_path;
        guidata1.last_funct_dir = guidata1.funct_path; 
    end
    if ischar(guidata1.funct_name); guidata1.funct_name = {guidata1.funct_name}; end
    check_num = numel(guidata1.funct_name);
    if check_num>1
        guidata1.multi_select = true; 
        guidata1.n_4D = check_num; 
    else
        guidata1.multi_select = false; 
    end
    [~,~,ext] = fileparts(guidata1.funct_name{1});
    for ix = 1:3
        if ~isempty(strfind(guidata1.extensions{ix},ext))
            type = ix;
            others = setdiff(1:3,ix);
            break;
        end
    end
    guidata1.extensions = guidata1.extensions([type,others]);
    if guidata1.last_funct_dir~=0
        guidata1.select_funct_text.String = fullfile(guidata1.last_funct_dir,guidata1.funct_name{1});
        guidata1.funct_name_string = guidata1.select_funct_text.String;
    end
    guidata(main_figure,guidata1)
end

function funct_input_box(~,~,~)
    guidata1 = guidata(main_figure);
    if ~strcmp(guidata1.select_funct_text.String,'')
        guidata1.funct_name_string = guidata1.select_funct_text.String;
        guidata(main_figure,guidata1)
    end
end

function load_funct_callback(~,~,~)
    guidata1 = guidata(main_figure);
    % Clear old data:
    guidata1.ROI_ind = [];
    guidata1.ROI_ind_back = [];
    guidata1.numvox = [];
    guidata1.numvox_background = [];
    guidata1.roi_sig = [];
    guidata1.roi_sig_filt = [];
    guidata1.roi_background = [];
    guidata1.roi_background_filt = [];
    guidata1.sig_detrended = [];
    guidata1.background_detrended = [];
    guidata1.funct_SNR_sd = [];
    guidata1.SFNR = [];
    guidata1.tSNR = [];
    guidata1.tSBNR = [];
    guidata1.SFNR_all = [];
    guidata1.signal_amp_all = [];
    guidata1.funct_SNR_sd_all = [];
    guidata1.funct_SNR_sd = [];
    % Begin:
    if ~isempty(guidata1.funct_name)
        if ~exist(fullfile(guidata1.funct_path,guidata1.funct_name{1}),'file')
        	errordlg('Functional image not found. Please respecify.'); return;
        end
    else
        errordlg('Please specify a functional image.'); return;
    end
    try
        disp('Loading functional series...')
        guidata1.funct = load_nii(fullfile(guidata1.funct_path,guidata1.funct_name{1}));
        guidata1.apply_header = 1;
    catch
        guidata1.funct = load_untouch_nii(fullfile(guidata1.funct_path,guidata1.funct_name{1}));
        warning('Non-orthogonal shearing detected in affine matrix. Successfully loaded raw image data.') 
        guidata1.apply_header = 0;
    end
    guidata1.funct_dim = guidata1.funct.hdr.dime.dim(2:4);
%     guidata1.funct_TR = guidata1.funct.hdr.dime.pixdim(5);
    if guidata1.multi_select % if 3D image series
        funct_4D = zeros([guidata1.funct_dim,guidata1.n_4D]);
        funct_4D(:,:,:,1) = guidata1.funct.img;
        hWait = waitbar(0,'Loading functional series...');
        try
            if guidata1.apply_header
                for i = 2:guidata1.n_4D
                    waitbar(i/guidata1.n_4D,hWait);
                    img = load_nii(fullfile(guidata1.funct_path,guidata1.funct_name{i}));
                    funct_4D(:,:,:,i) = img.img;
                end
            else
                for i = 2:guidata1.n_4D
                    waitbar(i/guidata1.n_4D,hWait);
                    img = load_untouch_nii(fullfile(guidata1.funct_path,guidata1.funct_name{i}));
                    funct_4D(:,:,:,i) = img.img;
                end
            end; close(hWait);
            guidata1.funct.img = funct_4D;
            clear funct_4D
        catch; close(hWait); errordlg('Error concatanating 3D functional series.'); return;
        end
    else
        guidata1.n_4D = size(guidata1.funct.img,4);
    end
    guidata1.funct.img = single(guidata1.funct.img); % change data type for later functions
    % Initialize 3D Mean Funct:
    guidata1.funct3D = guidata1.funct;
    guidata1.funct3D.img = mean(guidata1.funct.img,4);
    guidata1.funct3D.hdr.dime.dim(1) = 3;
    guidata1.funct3D.hdr.dime.dim(5) = 1;
    guidata(main_figure,guidata1)
end

function funct_ROI_spec(~, ~, ROI_type) 
    guidata1 = guidata(main_figure);
    % Clear Old Data:
    guidata1.ROI_ind = [];
    guidata1.ROI_ind_back = [];
    guidata1.numvox = [];
    guidata1.numvox_background = [];
    guidata1.roi_sig = [];
    guidata1.roi_sig_filt = [];
    guidata1.roi_background = [];
    guidata1.roi_background_filt = [];
    guidata1.sig_detrended = [];
    guidata1.background_detrended = [];
    guidata1.funct_SNR_sd = [];
    guidata1.SFNR = [];
    % Begin:
    if ~isempty(guidata1.funct)
        hwait = waitbar(.1,'Creating ROIs...');
        switch ROI_type
            case 1 %'Image Region' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp('Automatically calculating ROI locations based on standard image regions.')
                guidata1.funct_ROI_type = 'Automatic - Image Region';
                if isgraphics(guidata1.h); close(guidata1.h); end
                % Initialize funct_ROI:
                funct_ROI = zeros(guidata1.funct_dim);
                % Signal ROI:
                center_x = round(guidata1.funct_dim(1)/2);
                center_y = round(guidata1.funct_dim(2)/2);
                center_z = round(guidata1.funct_dim(3)/2);
                central_x = round(.1*guidata1.funct_dim(1));
                central_y_upper = round(.2*guidata1.funct_dim(2));
                central_y_lower = round(.06*guidata1.funct_dim(2));
                central_z = round(.15*guidata1.funct_dim(3));
                funct_ROI(center_x-central_x:center_x+central_x,...
                    center_y-central_y_lower:center_y+central_y_upper,...
                    center_z-central_z:center_z+central_z) = 1;
                % Noise ROIs
                border_x = round(.08*guidata1.funct_dim(1));
                border_y = round(.75*guidata1.funct_dim(2));
                funct_ROI(3:border_x,3:border_y,2:guidata1.funct_dim(3)-1) = 2;
                funct_ROI(guidata1.funct_dim(1)-border_x:guidata1.funct_dim(1)-2,...
                    3:border_y,2:guidata1.funct_dim(3)-1) = 2;
            case 2 %'Signal Intensity' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                guidata1.funct_ROI_type = 'Automatic - Signal Intensity';
                if isgraphics(guidata1.h); close(guidata1.h); end
                % Initialize funct_ROI:
                funct_ROI = zeros(guidata1.funct_dim);
                % Determine Thresholds:
                intensities = sort(guidata1.funct3D.img(:));
                sig_min_perc = .75; sig_max_perc = .995;
                sig_min = intensities(round(sig_min_perc*length(intensities)));
                sig_max = intensities(round(sig_max_perc*length(intensities)));
                noise_perc = .7; 
                noise_thresh = intensities(round(noise_perc*length(intensities)));
                % Signal ROI:
                funct_ROI(((guidata1.funct3D.img>=sig_min) + (guidata1.funct3D.img<=sig_max))==2) = 1;
                % Noise ROI:
                funct_ROI(guidata1.funct3D.img<=noise_thresh) = 2;
                disp(['Creating ROIs - Intensity thresholds for signal ROI: ',num2str(sig_min),'-',num2str(sig_max),...
                    '. Maximum intensity for background ROI: ',num2str(noise_thresh)])
            case 3 %'Intensity Threshold' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                guidata1.funct_ROI_type = 'Manual - Intensity Threshold';
                % Initialize funct_ROI:
                funct_ROI = zeros(guidata1.funct_dim);
                max_sig = max(guidata1.funct.img(:));
                % Specify thresholds for inclusive mask:
                prompt = {sprintf('Enter intensity thresholds for SIGNAL ROI: \n\nMin:'),'Max:'};
                dlg_title = 'Signal ROI'; num_lines = [1,50;1,50]; defaultans = {'0',num2str(max_sig)};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                if isempty(answer)
                    disp('User cancelled intensity threshold specification');
                    delete(hwait);
                    return;
                end
                sig_min = str2double(answer(1)); sig_max = str2double(answer(2));
                prompt = {'Enter maximum intensity for BACKGROUND ROI: '};
                dlg_title = 'Background ROI'; num_lines = [1,50;]; defaultans = {num2str(max_sig)};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                if isempty(answer)
                    disp('User cancelled intensity threshold specification');
                    delete(hwait);
                    return;
                end
                noise_thresh = str2double(answer(1));
                % Signal ROI:
                funct_ROI(((guidata1.funct3D.img>=sig_min) + (guidata1.funct3D.img<=sig_max))==2) = 1;
                % Noise ROI:
                funct_ROI(guidata1.funct3D.img<=noise_thresh) = 2;
                disp('Calculting ROI locations based on specified signal intensities.')
            case 4 %'Manual Tracing' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                guidata1.funct_ROI_type = 'Manual - Tracing';
                if isgraphics(guidata1.h); close(guidata1.h); end
                instructions_text = sprintf(['First, draw a SIGNAL ROI (green). ',...
                    'Second, draw a BACKGROUND/NOISE ROI (red). ',...
                    'To accomplish this, select "Draw" menu -> "Draw Color". ',...
                    'Once a color is selected, simply trace the outline of a shape on the image, ',...
                    'and this shape will automatically fill with color. Alternatively, select ',...
                    '"Draw" menu -> "Shapes" to draw using a specified shape. ',...
                    'Click and drag on the image to create shapes. To navigate between slices, ',...
                    'use the up and down arrow keys. An ROI drawn on one slice can be propagated through to other ',...
                    'slices by selecting "Draw" menu -> "Edit Drawing" -> "Propagate Through Slices" -> choose slice range. ',...
                    'Note that this will only apply to the most recent drawing. When you are finished drawing ',...
                    'select "Confirm ROIs" -> "Confirm". The ROIs will be automatically saved.']);
                h_mbox = msgbox(instructions_text,'ROI Tracing Instructions');
                uiwait(h_mbox,60)
                [guidata1.h,~,drawing_idx] = neuroimage_editor_tbx('apply_header',guidata1.apply_header,...
                    'background',guidata1.funct3D,'colorbar',1,'title_on',1,'axis_tick_on',1); %#ok
                uiwait(guidata1.h.figure,120) % wait until neuroimage_editor closes, or 60 sec
                drawing_idx = evalin('base','drawing_idx'); % pull last_draw in from 'base' workspace
                drawing_idx = permute(drawing_idx,[2,1,3]);
                funct_ROI = zeros(guidata1.funct_dim);
                funct_ROI(drawing_idx(:)==2) = 1;
                funct_ROI(drawing_idx(:)==3) = 2;
            case 5 %'Inclusive Mask' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Get input filename:
                cd(guidata1.last_funct_dir);
                [mask_name, mask_path] = uigetfile({'*.img';'*.nii';'*nii.gz'},...
                    'Select Inclusive Mask:','MultiSelect','off');
                cd(guidata1.first_dir);
                if mask_path ~= 0 % if not 'cancel'
                    guidata1.funct_ROI_type = 'Inclusive Mask';
                    try
                        if guidata1.apply_header
                            ROI_img = load_nii(fullfile(mask_path,mask_name));
                        else
                            ROI_img = load_untouch_nii(fullfile(mask_path,mask_name));
                        end
                    catch
                        errordlg('ERROR: Mask failed to load.'); return;
                    end
                    % Initialize funct_ROI (sig = 1, background = 2):
                    funct_ROI = ones(guidata1.funct_dim)*2;
                    if islogical(ROI_img.img)
                        funct_ROI(ROI_img.img) = 1;
                    else
                        funct_ROI(ROI_img.img==1) = 1;
                    end
                    disp('Successfully loaded inclusive mask.')
                else
                    disp('User cancelled action.'); delete(hwait); return;
                end
            case 6 %'Exclusive Mask' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                % Get input filename:
                cd(guidata1.last_funct_dir);
                [mask_name, mask_path] = uigetfile({'*.img';'*.nii';'*nii.gz'},...
                    'Select Exclusive Mask:','MultiSelect','off');
                cd(guidata1.first_dir);
                if mask_path ~= 0 % if not 'cancel'
                    guidata1.funct_ROI_type = 'Exclusive Mask';
                    try
                        if guidata1.apply_header
                            ROI_img = load_nii(fullfile(mask_path,mask_name));
                        else
                            ROI_img = load_untouch_nii(fullfile(mask_path,mask_name));
                        end
                    catch
                        errordlg('ERROR: Mask failed to load.'); return;
                    end
                    % Initialize funct_ROI (sig = 1, background = 2):
                    funct_ROI = ones(guidata1.funct_dim)*1;
                    if islogical(ROI_img.img)
                        funct_ROI(ROI_img.img) = 2;
                    else
                        funct_ROI(ROI_img.img==1) = 2;
                    end
                    disp('Successfully loaded exclusive mask.')
                else
                    disp('User cancelled action.'); delete(hwait); return;
                end
        end
        waitbar(.3,hwait);
        % Parse ROI Data:
        guidata1.ROI_ind = find(funct_ROI(:)==1);
        guidata1.ROI_ind_back = find(funct_ROI(:)==2);
        guidata1.funct_mean_signal = mean(guidata1.funct3D.img(guidata1.ROI_ind));
        % Form Signal ROI voxel x time matrix:
        if ~isempty(guidata1.ROI_ind)
            [guidata1.ROI_x,guidata1.ROI_y,guidata1.ROI_z] = ind2sub(guidata1.funct_dim,guidata1.ROI_ind);
            guidata1.numvox = length(guidata1.ROI_ind); 
            guidata1.roi_sig = zeros(guidata1.n_4D,guidata1.numvox);
            for i = 1:guidata1.numvox
                guidata1.roi_sig(:,i) = guidata1.funct.img(guidata1.ROI_x(i),guidata1.ROI_y(i),guidata1.ROI_z(i),:);
            end
            waitbar(.5,hwait);
        end
        % Form Background ROI voxel x time matrix (if exists):
        if ~isempty(guidata1.ROI_ind_back)
            [ROI_x_back,ROI_y_back,ROI_z_back] = ind2sub(guidata1.funct_dim,guidata1.ROI_ind_back);
            guidata1.numvox_background = numel(guidata1.ROI_ind_back); 
            guidata1.roi_background = zeros(guidata1.n_4D,guidata1.numvox_background);
            waitbar(.8,hwait);
            for i = 1:guidata1.numvox_background
                guidata1.roi_background(:,i) = guidata1.funct.img(ROI_x_back(i),...
                    ROI_y_back(i),ROI_z_back(i),:);
            end
        end
        waitbar(1,hwait);
    end
    disp('Finished creating ROIs.')
    delete(hwait)
    guidata(main_figure,guidata1)        
end

function funct_view_funct(~,~,~)
    guidata1 = guidata(main_figure);
    if isgraphics(guidata1.h); close(guidata1.h); end
    if ~isempty(guidata1.funct3D) 
        handles = nifti_studio('apply_header',guidata1.apply_header,...
            'background',guidata1.funct3D,...
            'colorbar_on',1,'title_on',1,'axis_tick_on',1);
        guidata1.h = handles.figure;
        guidata(main_figure,guidata1) 
    else
        errordlg('Please load a functional image.'); return;
    end
end

function funct_view_ROIs(~,~,~)
    guidata1 = guidata(main_figure);
    if isgraphics(guidata1.h); close(guidata1.h); end
    if ~isempty(guidata1.funct3D) 
        if ~isempty(guidata1.ROI_ind)
            funct_ROI = reconstruct_image(ones(1,guidata1.numvox),guidata1.ROI_ind);
            funct_ROI.img(guidata1.ROI_ind_back) = 2;
            handles = nifti_studio('apply_header',guidata1.apply_header,...
                'background',guidata1.funct3D,...
                'overlay',funct_ROI,'colorbar_on',0,'title_on',0,'axis_tick_on',...
                0,'colormap','gray');
            guidata1.h = handles.figure;
            guidata(main_figure,guidata1) 
        else
            errordlg('Please Specify Signal & Noise ROIs.'); return;
        end
    else
        errordlg('Please load a functional image.'); return;
    end
end

function funct_detrend(~,~,~)
    guidata1 = guidata(main_figure);
    if isempty(guidata1.funct); errordlg('First load a functional image.'); return; end
    if isempty(guidata1.roi_sig); errordlg('First specify ROIs.'); return; end
    % Takes ~36.12 s per 50,000 voxels with 150 time points
    [guidata1.use_linear, cancelled] = waitTimeUserInput(guidata1.numvox);
    if cancelled; return; end
    % Detrend:
    if guidata1.use_linear
        disp('Performing linear detrending.')
        guidata1.sig_detrended = detrend(guidata1.roi_sig);
        if ~isempty(guidata1.roi_background)
            guidata1.background_detrended = detrend(guidata1.roi_background);
        end
    else
        disp('Performing quadratic detrending.')
        hwait = waitbar(0,'Detrending...');
        guidata1.sig_detrended = guidata1.roi_sig; % initialize
        x = (1:guidata1.n_4D)';
        for i = 1:guidata1.numvox
            if rem(i,1000)==0; waitbar(i/guidata1.numvox,hwait); end
            p = polyfit(x,guidata1.roi_sig(:,i),2);
            predicted = polyval(p,x);
            guidata1.sig_detrended(:,i) = guidata1.roi_sig(:,i)-predicted;
        end
        waitbar(1,hwait,'Done!'); pause(.1); delete(hwait);
        if ~isempty(guidata1.roi_background)
            hwait = waitbar(0,'Detrending background...');
            guidata1.background_detrended = guidata1.roi_background; % initialize
            x = (1:guidata1.n_4D)';
            for i = 1:guidata1.numvox_background
                if rem(i,1000)==0; waitbar(i/guidata1.numvox_background,hwait); end
                p = polyfit(x,guidata1.roi_background(:,i),2);
                predicted = polyval(p,x);
                guidata1.background_detrended(:,i) = guidata1.roi_background(:,i)-predicted;
            end
            waitbar(1,hwait,'Done!'); pause(.1); delete(hwait);
        end
    end
    disp('Finished detrending.')
    guidata(main_figure,guidata1) 
end

function funct_bandpass(~,~,~,Fs,HighPass,LowPass)
    guidata1 = guidata(main_figure);
    if ~isempty(guidata1.funct)
        if ~isempty(guidata1.ROI_ind)
            if nargin<4
                % Detect TR:
                try
                TR_detect = guidata1.funct.hdr.dime.pixdim(5);
                % Detect TR Units:
                switch bitand(guidata1.funct.hdr.dime.xyzt_units, 56)	
                    case 8
                      t_units = 's';
                    case 16
                      t_units = 'ms';
                    case 24
                      t_units = 'microseconds';	% microsecond
                    otherwise
                      t_units = 'unknown_units';
                end
                catch; TR_detect = [];
                end 
                % Prompt User: effective TR; lower bound (high-pass); upper-bound (low-pass)
                if ~isempty(TR_detect) && TR_detect ~= 0
                    prompt = {sprintf(['Detected TR = ',num2str(TR_detect),' ',num2str(t_units),...
                        '. Multiply by # repetitions per 3D image to obtain sampling time. \n\nSampling Time:']),...
                        'Min Frequency (Hz):','Max Frequency (Hz):'};
                else
                    prompt = {'Sampling Time (s):','Min Frequency (Hz):','Max Frequency (Hz):'};
                end
                dlg_title = 'Band-Pass Settings'; num_lines = [1,40;1,40;1,40]; defaultans = {'2','.01','.1'};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                if isempty(answer)
                    disp('User cancelled band-pass filtration.'); return;
                end
                TR = str2double(answer(1)); 
                guidata1.HighPass = str2double(answer(2));
                guidata1.LowPass = str2double(answer(3));
                frange = guidata1.LowPass-guidata1.HighPass;
                guidata1.fmax_display = guidata1.LowPass+frange*.3;  
                guidata1.Fs = 1/TR; % samples per second
            else
                guidata1.Fs = Fs;
                guidata1.HighPass = HighPass;
                guidata1.LowPass = LowPass;
            end
        else
            errordlg('First, specify ROIs.'); return;
        end
    else
        errordlg('First, load a functional image.'); return;
    end
    hwait = waitbar(.3,'Band pass filtering...');
    % Infinite Impulse Response (IIR) filter: 
    % See for more info: https://www.minidsp.com/applications/dsp-basics/fir-vs-iir-filtering
    if ~isempty(guidata1.sig_detrended)
        [guidata1.roi_sig_filt,~,~] = bst_bandpass_filtfilt(guidata1.sig_detrended',...
            guidata1.Fs, guidata1.HighPass, guidata1.LowPass, 0, 'iir');
    else
        [guidata1.roi_sig_filt,~,~] = bst_bandpass_filtfilt(guidata1.roi_sig',...
            guidata1.Fs, guidata1.HighPass, guidata1.LowPass, 0, 'iir');
    end
    waitbar(.7,hwait);
    if ~isempty(guidata1.background_detrended)
        [guidata1.roi_background_filt,~,~] = bst_bandpass_filtfilt(guidata1.background_detrended',...
            guidata1.Fs, guidata1.HighPass, guidata1.LowPass, 0, 'iir');
    else
        [guidata1.roi_background_filt,~,~] = bst_bandpass_filtfilt(guidata1.roi_background',...
            guidata1.Fs, guidata1.HighPass, guidata1.LowPass, 0, 'iir');
    end
    guidata1.roi_sig_filt = guidata1.roi_sig_filt';
    guidata1.roi_background_filt = guidata1.roi_background_filt';
    waitbar(1,hwait);
    guidata(main_figure,guidata1)
    if ~isempty(guidata1.sig_detrended)
        disp(['Successfully applied ',num2str(guidata1.HighPass),' - ',...
            num2str(guidata1.LowPass),' Hz band-pass filter to detrended time series.'])
    else
        disp(['Successfully applied ',num2str(guidata1.HighPass),' - ',...
            num2str(guidata1.LowPass),' Hz band-pass filter.'])
    end
    delete(hwait)
end

function funct_view_signals(~,~,~)
    guidata1 = guidata(main_figure);
    if isempty(guidata1.funct); errordlg('First load functional image.'); return; end
    if isempty(guidata1.roi_sig); errordlg('First specify ROIs.'); return; end
%     figure_pos = [.007*guidata1.screen_res(3),.25*guidata1.screen_res(4), ...
%     .988*guidata1.screen_res(3),.5*guidata1.screen_res(4)];
    h_signals = figure('units','normalized','Position',[.0065,.2491,.988,.5],...
        'Name','Randomly Selected Voxel Time Series',...
        'NumberTitle','off','Color',zeros(1,3),'MenuBar','none');
    guidata1.signum = randi(guidata1.numvox);
    guidata1.mean_centered = guidata1.roi_sig(:,guidata1.signum)-mean(guidata1.roi_sig(:,guidata1.signum));
    guidata1.opts = false(1,6); guidata1.opts(1) = true; 
    if ~isempty(guidata1.sig_detrended) 
        guidata1.opts(2) = true; 
    else
        guidata1.opts(2) = false; 
    end
    if ~isempty(guidata1.roi_sig_filt) 
        guidata1.opts(6) = true; 
    else
        guidata1.opts(6) = false; 
    end
    if ~isempty(guidata1.use_linear)
        guidata1.opts(4) = true; 
    else
        guidata1.opts(5) = true;
    end
    file_menu = uimenu(h_signals,'Label','File');
        uimenu(file_menu,'Label','Choose New Signal','Callback',@new_signal);
        uimenu(file_menu,'Label','Print','Callback',@print_signal);
    select_menu = uimenu(h_signals,'Label','Display Selected');
    if guidata1.opts(1)
        select_menu1 = uimenu(select_menu,'Label','Signal (Mean Centered) *','Callback',{@select_display,1}); 
    else
        select_menu1 = uimenu(select_menu,'Label','Signal (Mean Centered)','Callback',{@select_display,1});
    end
    if guidata1.opts(2)
        select_menu2 = uimenu(select_menu,'Label','Detrended (Linear) *','Callback',{@select_display,2}); 
    else
        select_menu2 = uimenu(select_menu,'Label','Detrended (Linear)','Callback',{@select_display,2});
    end
    if guidata1.opts(3)
        select_menu3 = uimenu(select_menu,'Label','Detrended (Quadratic) *','Callback',{@select_display,3}); 
    else
        select_menu3 = uimenu(select_menu,'Label','Detrended (Quadratic)','Callback',{@select_display,3});
    end
    if guidata1.opts(4)
        select_menu4 = uimenu(select_menu,'Label','Linear Trend *','Callback',{@select_display,4});    
    else
        select_menu4 = uimenu(select_menu,'Label','Linear Trend','Callback',{@select_display,4});
    end
    if guidata1.opts(5)
        select_menu5 = uimenu(select_menu,'Label','Quadratic Trend *','Callback',{@select_display,5});
    else
        select_menu5 = uimenu(select_menu,'Label','Quadratic Trend','Callback',{@select_display,5});
    end
    if ~isempty(guidata1.roi_sig_filt)
        if guidata1.opts(6)
            select_menu6 = uimenu(select_menu,'Label','Band-Pass Filtered *','Callback',{@select_display,6});
        else
            select_menu6 = uimenu(select_menu,'Label','Band-Pass Filtered','Callback',{@select_display,6});
        end
    end
    guidata(main_figure,guidata1)
    update_plot(1)
    function new_signal(~,~,~)
        guidata1 = guidata(main_figure);
        guidata1.signum = randi(guidata1.numvox);
        guidata1.mean_centered = guidata1.roi_sig(:,guidata1.signum)-mean(guidata1.roi_sig(:,guidata1.signum));
        guidata(main_figure,guidata1)
        update_plot(1)
    end
    function print_signal(~,~,~)
        % Identify File Extension:
        ext_opts = {'.png',' ';'.tiff',' ';'.bmp',' '};
        [title1,path1] = uiputfile(ext_opts,'Specify filename:');
        figure_title = fullfile(path1,title1);
        [~,~,file_ext] = fileparts(figure_title);
        % Identify File Extension:
        exts = {'-dpng','-dtiff','-dbmp'};
        for ix = 1:3
            if ~isempty(strfind(exts{ix},file_ext(2:end)))
                ext_ind = ix; break;
            end
        end
        % Print Figure:
        set(h_signals,'InvertHardcopy','off','PaperPositionMode','auto')
        print(h_signals,exts{ext_ind},'-r300','-loose',figure_title)
        disp(['Printed figure: ',figure_title])
    end
    function select_display(~,~,opt_num)
        guidata1 = guidata(main_figure);
        guidata1.opts(opt_num) = ~guidata1.opts(opt_num);
        if guidata1.opts(1)
            set(select_menu1,'Label','Signal (Mean Centered) *') 
        else
            set(select_menu1,'Label','Signal (Mean Centered)')
        end
        if guidata1.opts(2)
            set(select_menu2,'Label','Detrended (Linear) *'); 
        else
            set(select_menu2,'Label','Detrended (Linear)');
        end
        if guidata1.opts(3)
            set(select_menu3,'Label','Detrended (Quadratic) *'); 
        else
            set(select_menu3,'Label','Detrended (Quadratic)');
        end
        if guidata1.opts(4)
            set(select_menu4,'Label','Linear Trend *');    
        else
            set(select_menu4,'Label','Linear Trend');
        end
        if guidata1.opts(5)
            set(select_menu5,'Label','Quadratic Trend *');
        else
            set(select_menu5,'Label','Quadratic Trend');
        end
        if ~isempty(guidata1.roi_sig_filt)
            if guidata1.opts(6)
                set(select_menu6,'Label','Band-Pass Filtered *');
            else
                set(select_menu6,'Label','Band-Pass Filtered');
            end
        end
        guidata(main_figure,guidata1)
        update_plot
    end
    function update_plot(~)
        guidata1 = guidata(main_figure);
        if nargin==1
            if isgraphics(guidata1.ax1,'axes'); cla(guidata1.ax1); end
            guidata1.h_plots(1) = plot(guidata1.mean_centered,'r');
            guidata1.ax1 = gca; hold(guidata1.ax1,'on')
            guidata1.opts(1) = true;
            set(guidata1.ax1,'Position',[.035,.07,.955,.9]);
        end
        if guidata1.opts(1) && nargin==0
            if isgraphics(guidata1.h); close(guidata1.h); end
            guidata1.h_plots(1) = plot(guidata1.mean_centered,'r');
        elseif ~guidata1.opts(1)
            if ishandle(guidata1.h_plots(1)); set(guidata1.h_plots(1),'Visible','off'); end
        end
        if guidata1.opts(2)
            if ishandle(guidata1.h_plots(2)); delete(guidata1.h_plots(2)); end
            x = (1:1:guidata1.n_4D)';
            p = polyfit(x,guidata1.mean_centered,1);
            predicted = polyval(p,x);
            detrended_linear = guidata1.mean_centered-predicted;
            guidata1.h_plots(2) = plot(guidata1.ax1,detrended_linear,'c');
        else
            if ishandle(guidata1.h_plots(2)); set(guidata1.h_plots(2),'Visible','off'); end
        end
        if guidata1.opts(3)
            if ishandle(guidata1.h_plots(3)); delete(guidata1.h_plots(3)); end
            x = (1:1:guidata1.n_4D)';
            p = polyfit(x,guidata1.mean_centered,2);
            predicted = polyval(p,x);
            detrended_linear = guidata1.mean_centered-predicted;
            guidata1.h_plots(3) = plot(guidata1.ax1,detrended_linear,'m');
        else
            if ishandle(guidata1.h_plots(3)); set(guidata1.h_plots(3),'Visible','off'); end
        end
        if guidata1.opts(4)
            x = (1:1:guidata1.n_4D)';
            p = polyfit(x,guidata1.mean_centered,1);
            predicted = polyval(p,x);
            if ishandle(guidata1.h_plots(4)); delete(guidata1.h_plots(4)); end
            guidata1.h_plots(4) = plot(guidata1.ax1,predicted,'g');
        else
            if ishandle(guidata1.h_plots(4)); set(guidata1.h_plots(4),'Visible','off'); end
        end
        if guidata1.opts(5)
            x = (1:1:guidata1.n_4D)';
            p = polyfit(x,guidata1.mean_centered,2);
            predicted = polyval(p,x);
            if ishandle(guidata1.h_plots(5)); delete(guidata1.h_plots(5)); end
            guidata1.h_plots(5) = plot(guidata1.ax1,predicted,'b');
        else
            if ishandle(guidata1.h_plots(5)); set(guidata1.h_plots(5),'Visible','off'); end
        end
        if ~isempty(guidata1.roi_sig_filt)
            if guidata1.opts(6)
                if ishandle(guidata1.h_plots(6)); delete(guidata1.h_plots(6)); end
                guidata1.h_plots(6) = plot(guidata1.ax1,guidata1.roi_sig_filt(:,guidata1.signum),'k');
            else
                if ishandle(guidata1.h_plots(6)); set(guidata1.h_plots(6),'Visible','off'); end
            end
        end
        ax_color = [.96,.96,.96];
        if isgraphics(guidata1.ax1,'axes')
            set(guidata1.ax1,'color',ax_color,'box','on','XColor',ax_color,...
                'YColor',ax_color,...
                'BoxStyle','full','XMinorGrid','on','YMinorGrid','on',...
                'XColorMode','manual','YColorMode','manual',...
                'GridColor',ax_color,'GridColorMode','manual',...
                'MinorGridColor',ax_color,'MinorGridColorMode','manual',...
                'XGrid','on','YGrid','on','TickDir','out');
        end
        objs = []; for i = 1:6; if guidata1.opts(i); objs = [objs,guidata1.h_plots(i)]; end; end; %#ok
        h_leg = legend(guidata1.ax1,objs,guidata1.legend_str(guidata1.opts),'location','northeast','FontSize',8);
        leg_pos = get(h_leg,'Position');
        leg_pos(1) = .99-leg_pos(3); leg_pos(2) = .97-leg_pos(4);
        set(h_leg,'Position',leg_pos); % ax_pos = [.035,.07,.955,.9]
        guidata(main_figure,guidata1);
    end
end

function funct_view_spectrum(~,~,~)
    % Load Data and Check
    guidata1 = guidata(main_figure);
    if isempty(guidata1.funct); errordlg('First load functional image.'); return; end
    if isempty(guidata1.roi_sig); errordlg('First specify ROIs.'); return; end
    if isempty(guidata1.Fs)
        prompt = 'Sampling Time (s):';
        dlg_title = 'Sample Time'; num_lines = [1,30]; defaultans = {'2'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        if isempty(answer); disp('User cancelled action.'); return; end
        guidata1.Fs = 1/str2double(answer(1));
        guidata(main_figure,guidata1);
    end
   spec_type = questdlg('Choose power spectral density estimation method: ',...
        'Power Spectrum Method','Welch''s Power Spectrum (faster)',...
        'Multitaper PSD (slower, more accurate)','Cancel',...
        'Welch''s Power Spectrum (faster)');
    switch spec_type
        case 'Welch''s Power Spectrum (faster)'
            spec_method=1;
        case 'Multitaper PSD (slower, more accurate)'
            spec_method=2;
        case 'Cancel'
            disp('User cancelled operation.'); return;
    end
    % Initialize Figure
%     power_figure_pos = [.006*guidata1.screen_res(3), .1*guidata1.screen_res(4),...
%         .989*guidata1.screen_res(3), .825*guidata1.screen_res(4)];
    h_power = figure('units','normalized','Position',[.0055,.0991,.989,.825],'MenuBar','none',...
        'Name','Power Spectrum','NumberTitle','off','Color',[0,0,0]);
    ax_pos4 = [.04,.567,.45,.4; .535,.567,.45,.4; .04,.065,.45,.4; .535,.065,.45,.4];
    ax_pos2 = [.04,.567,.94,.4; .04,.065,.94,.4];
    ax_pos1 = [.04,.07,.94,.897];
    ax = [12.3982,14.334,18.484,18.436];
    if guidata1.numvox > guidata1.numvoxel_limit 
        rand_selection = randperm(guidata1.numvox,guidata1.numvoxel_limit);
    else
        rand_selection = [];
    end
    % Check Processing Status, and Process Spectrums:
    plot_menus = false(4,23);
    plot_menus(1,8) = true; plot_menus(2,9) = true; 
    plot_menus(3,[16,20]) = true; plot_menus(4,[18,22]) = true;
    process_status = [~isempty(guidata1.roi_sig_filt),...
        ~isempty(guidata1.roi_background),~isempty(guidata1.roi_background_filt)];
    if sum(process_status)==3 % all present
        ind_use = 1:4; 
    elseif ~process_status(1) && process_status(2) && ~process_status(3) % background, no filter
        ind_use = [1,3]; plot_menus(2,9) = false; plot_menus(2,10) = true; % add background to plot 2
        plot_menus(4,[18,22]) = false; plot_menus(2,9) = false;  plot_menus(4,2) = true;
    elseif process_status(1) && ~process_status(2) && ~process_status(3) % no background
        ind_use = 1:2; plot_menus(3,20) = false; plot_menus(4,22) = false;
    elseif process_status(1) && process_status(2) && ~process_status(3) % bandpass sig, not background
        ind_use = 1:3; plot_menus(4,22) = false;
    end
    % Calculate Spectrums
    specs = calc_specs(process_status);
    % Initialize File Menus:
    file_menu = uimenu(h_power,'Label','File');
        uimenu(file_menu,'Label','Print','Callback',@print_spec);
        uimenu(file_menu,'Label','Draw ROI for PowerSpec','Callback',@drawROI4Power)
    h_plot_menus = zeros(4,23);
    plot_menu_labels = {'Subplot 1','Subplot 2','Subplot 3','Subplot 4',...
        'Disable','Power Spec','Mean Power Spec','Pre-Filter (Signal ROI)',...
        'Post-Filter (Signal ROI)','Pre-Filter (Background ROI)',...
        'Post-Filter (Background ROI)','Mean Pre-Filter (Signal ROI)',...
        'Mean Post-Filter (Signal ROI)','Mean Pre-Filter (Background ROI)',...
        'Mean Post-Filter (Background ROI)','CI (On)','CI (Off)'};
    legend_str = {'Mean Pre-Filter (Signal ROI)','95% CI (Pre-Filter, Signal ROI)',...
        'Mean Post-Filter (Signal ROI)','95% CI (Post-Filter, Signal ROI)',...
        'Mean Pre-Filter (Background ROI)','95% CI (Pre-Filter, Background ROI)',...
        'Mean Post-Filter (Background ROI)','95% CI (Post-Filter, Background ROI)'};
    xlimits = nan(2,4); ylimits = nan(2,4); kvec = [1,0];
    for i = 1:4
        h_plot_menus(i,1) = uimenu('Label',plot_menu_labels{i});
        if ~plot_menus(i,2)
            h_plot_menus(i,2) = uimenu(h_plot_menus(i,1),'Label',...
                plot_menu_labels{5},'Callback',@select_display);
        else
            h_plot_menus(i,2) = uimenu(h_plot_menus(i,1),'Label',...
                [plot_menu_labels{5},' *'],'Callback',@select_display);
        end
        h_plot_menus(i,3) = uimenu(h_plot_menus(i,1),'Label',plot_menu_labels{6});
        h_plot_menus(i,4) = uimenu(h_plot_menus(i,1),'Label',plot_menu_labels{7});
        for j = ind_use
            if plot_menus(i,j+7)
                h_plot_menus(i,j+7) = uimenu(h_plot_menus(i,3),'Label',...
                    [plot_menu_labels{j+7},' *'],'Callback',@select_display);
            else
                h_plot_menus(i,j+7) = uimenu(h_plot_menus(i,3),'Label',...
                    plot_menu_labels{j+7},'Callback',@select_display);
            end
        end
        for j = ind_use
            h_plot_menus(i,j+11) = uimenu(h_plot_menus(i,4),'Label',plot_menu_labels{j+11});
            for k = 1:2
                ix1 = 15+j*2-kvec(k);
                if plot_menus(i,ix1)
                    h_plot_menus(i,ix1) = uimenu(h_plot_menus(i,j+11),...
                        'Label',[plot_menu_labels{k+15},' *'],'Callback',@select_display);
                else
                    h_plot_menus(i,ix1) = uimenu(h_plot_menus(i,j+11),...
                        'Label',plot_menu_labels{k+15},'Callback',@select_display);
                end
            end
        end
        h_ax_lim = uimenu(h_plot_menus(i,1),'Label','Axis Limits');
        h_xlim = uimenu(h_ax_lim,'Label','X-Limits');
        h_ylim = uimenu(h_ax_lim,'Label','Y-Limits');
        uimenu(h_xlim,'Label','Autoscale','Callback',{@set_limits,i,'x'})
        uimenu(h_xlim,'Label','Specify','Callback',{@set_limits,i,'x'})
        uimenu(h_ylim,'Label','Autoscale','Callback',{@set_limits,i,'y'})
        uimenu(h_ylim,'Label','Specify','Callback',{@set_limits,i,'y'})
    end
    update_plot
    % Update Selection:
    function select_display(hObj,~)
        ind = h_plot_menus==hObj;
        col_ind = find(sum(ind,1)); row_ind = find(sum(ind,2));
        val = ~plot_menus(ind);
        plot_menus(ind) = val;
        str = get(h_plot_menus(ind),'Label');
        if val; str = [str,' *'];
        else
            str = str(1:end-2);
        end; set(h_plot_menus(ind),'Label',str)
        % Disable Others if Necessary
        if plot_menus(row_ind,2) && col_ind~=2 && val % Enable from Disable
            plot_menus(row_ind,2) = false;
            str = get(h_plot_menus(row_ind,2),'Label'); 
            str = str(1:end-2);
            set(h_plot_menus(row_ind,2),'Label',str)
        end
        if col_ind==2 && val % Disable Button
            for ix = (col_ind+1):23
                if plot_menus(row_ind,ix)
                    plot_menus(row_ind,ix) = false;
                    str = get(h_plot_menus(row_ind,ix),'Label'); 
                    str = str(1:end-2);
                    set(h_plot_menus(row_ind,ix),'Label',str)
                end
            end
        elseif val && col_ind>7 && col_ind<12 % plot a full power spec
            for ix = setdiff([8:11,16:23],col_ind)
                if plot_menus(row_ind,ix)
                    plot_menus(row_ind,ix) = false;
                    str = get(h_plot_menus(row_ind,ix),'Label'); 
                    str = str(1:end-2);
                    set(h_plot_menus(row_ind,ix),'Label',str)
                end
            end
        elseif val && col_ind>15 % plot a mean power spec
            for ix = 8:11
                if plot_menus(row_ind,ix)
                    plot_menus(row_ind,ix) = false;
                    str = get(h_plot_menus(row_ind,ix),'Label'); 
                    str = str(1:end-2);
                    set(h_plot_menus(row_ind,ix),'Label',str)
                end
            end
            if col_ind==16 
                plot_menus(row_ind,17) = false;
            elseif col_ind==17
                plot_menus(row_ind,16) = false;
            elseif col_ind==18
                plot_menus(row_ind,19) = false;           
            elseif col_ind==19
                plot_menus(row_ind,18) = false;
            elseif col_ind==20
                plot_menus(row_ind,21) = false;
            elseif col_ind==21
                plot_menus(row_ind,20) = false;
            elseif col_ind==22
                plot_menus(row_ind,23) = false;
            elseif col_ind==23
                plot_menus(row_ind,22) = false;
            end
        end 
        update_plot
    end
    % Update Plot:
    function update_plot
        if isgraphics(h_power,'figure')
            guidata1 = guidata(main_figure);
            % Determine *Which Axes & *Where to Place:
            for ix = 1:4; if isgraphics(ax(ix),'axes'); delete(ax(ix)); end; end
            ax = [12.3982,14.334,18.3843,16.4368];
            if ~plot_menus(1,2) && ~plot_menus(2,2) && (~plot_menus(3,2) || ~plot_menus(4,2))
                ax(1) = axes('parent',h_power); set(ax(1),'Position',ax_pos4(1,:)); hold(ax(1),'on');
                ax(2) = axes('parent',h_power); set(ax(2),'Position',ax_pos4(2,:)); hold(ax(2),'on')
            elseif ~plot_menus(1,2) && plot_menus(2,2) && (~plot_menus(3,2) || ~plot_menus(4,2))
                ax(1) = axes('parent',h_power); set(ax(1),'Position',ax_pos2(1,:)); hold(ax(1),'on')
            elseif plot_menus(1,2) && ~plot_menus(2,2) && (~plot_menus(3,2) || ~plot_menus(4,2))
                ax(2) = axes('parent',h_power); set(ax(2),'Position',ax_pos2(1,:)); hold(ax(2),'on')
            elseif ~plot_menus(1,2) && plot_menus(2,2) && plot_menus(3,2) && plot_menus(4,2)
                ax(1) = axes('parent',h_power); set(ax(1),'Position',ax_pos1); hold(ax(1),'on')
            elseif plot_menus(1,2) && ~plot_menus(2,2) && plot_menus(3,2) && plot_menus(4,2) 
                ax(2) = axes('parent',h_power); set(ax(2),'Position',ax_pos1); hold(ax(2),'on')
            elseif ~plot_menus(1,2) && ~plot_menus(2,2) && plot_menus(3,2) && plot_menus(4,2) 
                ax(1) = axes('parent',h_power); set(ax(1),'Position',ax_pos2(1,:)); hold(ax(1),'on')
                ax(2) = axes('parent',h_power); set(ax(2),'Position',ax_pos2(2,:)); hold(ax(2),'on')
            end
            if ~plot_menus(3,2) && ~plot_menus(4,2) && (~plot_menus(1,2) || ~plot_menus(2,2))
                ax(3) = axes('parent',h_power); set(ax(3),'Position',ax_pos4(3,:)); hold(ax(3),'on')
                ax(4) = axes('parent',h_power); set(ax(4),'Position',ax_pos4(4,:)); hold(ax(4),'on')
            elseif ~plot_menus(3,2) && plot_menus(4,2) && (~plot_menus(1,2) || ~plot_menus(2,2))
                ax(3) = axes('parent',h_power); set(ax(3),'Position',ax_pos2(2,:)); hold(ax(3),'on')
            elseif plot_menus(3,2) && ~plot_menus(4,2) && (~plot_menus(1,2) || ~plot_menus(2,2))
                ax(4) = axes('parent',h_power); set(ax(4),'Position',ax_pos2(2,:)); hold(ax(4),'on')
            elseif plot_menus(1,2) && plot_menus(2,2) && ~plot_menus(3,2) && plot_menus(4,2)
                ax(3) = axes('parent',h_power); set(ax(3),'Position',ax_pos1); hold(ax(3),'on')
            elseif plot_menus(1,2) && plot_menus(2,2) && plot_menus(3,2) && ~plot_menus(4,2) 
                ax(4) = axes('parent',h_power); set(ax(4),'Position',ax_pos1); hold(ax(4),'on')
            elseif plot_menus(1,2) && plot_menus(2,2) && ~plot_menus(3,2) && ~plot_menus(4,2) 
                ax(3) = axes('parent',h_power); set(ax(3),'Position',ax_pos2(1,:)); hold(ax(3),'on')
                ax(4) = axes('parent',h_power); set(ax(4),'Position',ax_pos2(2,:)); hold(ax(4),'on')
            end
            ax_color = [.96,.96,.96];
            for ix = 1:4
                if isgraphics(ax(ix),'axes')
                    set(ax(ix),'color',ax_color,'box','on','XColor',ax_color,...
                        'YColor',ax_color,...
                        'BoxStyle','full','XMinorGrid','on','YMinorGrid','on',...
                        'XColorMode','manual','YColorMode','manual',...
                        'GridColor',ax_color,'GridColorMode','manual',...
                        'MinorGridColor',ax_color,'MinorGridColorMode','manual',...
                        'XGrid','on','YGrid','on','TickDir','out');
                end
            end
            % Determine What to Plot:
            for ix = find(~plot_menus(:,2))'
                ind_plots = find(plot_menus(ix,:));
                % Check for Full Power Specs vs. Mean Specs:
                full_specs = ((ind_plots>7) + (ind_plots<12))==2;
                mean_specs = ind_plots>15;
                if any(full_specs)
                    which1 = ind_plots(full_specs)-7;
                    switch which1
                        case 1
                            if isempty(rand_selection)
                                plot(ax(ix),specs.f,specs.sig.prefilt);
                                if spec_method==1
                                    title(ax(ix),'Signal ROI Pre-filter Power Spectrum');
                                elseif spec_method==2
                                    title(ax(ix),'Signal ROI Pre-filter PSD');
                                end
                            else
                                plot(ax(ix),specs.f,specs.sig.prefilt(:,rand_selection));
                                if spec_method==1
                                    title(ax(ix),sprintf(['Signal ROI Pre-filter Power Spectrum (',...
                                        num2str(guidata1.numvoxel_limit),' random voxels)']))
                                elseif spec_method==2
                                    title(ax(ix),sprintf(['Signal ROI Pre-filter PSD (',...
                                        num2str(guidata1.numvoxel_limit),' random voxels)']))
                                end
                            end
                        case 2
                            if isempty(rand_selection)
                                plot(ax(ix),specs.f,specs.sig.postfilt);
                                if spec_method==1 
                                    title(ax(ix),'Signal ROI Post-filter Power Spectrum');
                                elseif spec_method==2 
                                    title(ax(ix),'Signal ROI Post-filter PSD');
                                end
                            else
                                plot(ax(ix),specs.f,specs.sig.postfilt(:,rand_selection));
                                if spec_method==1 
                                    title(ax(ix),sprintf(['Signal ROI Post-filter Power Spectrum (',...
                                        num2str(guidata1.numvoxel_limit),' random voxels)']))
                                elseif spec_method==2 
                                    title(ax(ix),sprintf(['Signal ROI Post-filter PSD (',...
                                        num2str(guidata1.numvoxel_limit),' random voxels)']))
                                end
                            end    
                        case 3
                            if isempty(rand_selection) || (guidata1.numvox_background <= guidata1.numvoxel_limit)
                                plot(ax(ix),specs.f,specs.background.prefilt);
                                if spec_method==1 
                                    title(ax(ix),'Background ROI Pre-filter Power Spectrum');
                                elseif spec_method==2 
                                    title(ax(ix),'Background ROI Pre-filter PSD');
                                end
                            else
                                plot(ax(ix),specs.f,specs.background.prefilt(:,rand_selection));
                                if spec_method==1
                                    title(ax(ix),sprintf(['Background ROI Pre-filter Power Spectrum (',...
                                        num2str(guidata1.numvoxel_limit),' random voxels)']))
                                elseif spec_method==2
                                    title(ax(ix),sprintf(['Background ROI Pre-filter PSD (',...
                                        num2str(guidata1.numvoxel_limit),' random voxels)']))
                                end
                            end
                        case 4
                            if isempty(rand_selection) || (guidata1.numvox_background <= guidata1.numvoxel_limit)
                                plot(ax(ix),specs.f,specs.background.postfilt);
                                if spec_method==1 
                                    title(ax(ix),'Background ROI Post-filter Power Spectrum');
                                elseif spec_method==2 
                                    title(ax(ix),'Background ROI Post-filter PSD');
                                end
                            else
                                plot(ax(ix),specs.f,specs.background.postfilt(:,rand_selection));
                                if spec_method==1
                                    title(ax(ix),sprintf(['Background ROI Post-filter Power Spectrum (',...
                                        num2str(guidata1.numvoxel_limit),' random voxels)']))
                                elseif spec_method==2
                                    title(ax(ix),sprintf(['Background ROI Post-filter PSD (',...
                                        num2str(guidata1.numvoxel_limit),' random voxels)']))
                                end
                            end      
                    end
                    h_xlab = xlabel(ax(ix),'Frequency (Hz)'); 
                    pos_xlab = get(h_xlab,'Position');
                    pos_xlab(2) = pos_xlab(2)+.012; set(h_xlab,'Position',pos_xlab)
                    if spec_method==1
                        ylabel(ax(ix),'Power (a.u.)')
                    elseif spec_method==2; ylabel(ax(ix),'Power/Hz')
                    end
                elseif any(mean_specs)
                     which2 = ind_plots(mean_specs)-15;
                     h1 = zeros(1,8); h_plots = [];
                     if any(which2==1)
                         h1(1) = plot(ax(ix),specs.f,specs.sig.mean_prefilt,'r');
                         h1(2) = plot(ax(ix),specs.f,specs.sig.CI_lower_prefilt,'r--');
                         plot(ax(ix),specs.f,specs.sig.CI_upper_prefilt,'r--');
                         h_plots = [h_plots,h1(1),h1(2)]; %#ok
                     elseif any(which2==2)
                         h1(1) = plot(ax(ix),specs.f,specs.sig.mean_prefilt,'r');
                         h_plots = [h_plots,h1(1)]; %#ok
                     end
                     if any(which2==3)
                         h1(3) = plot(ax(ix),specs.f,specs.sig.mean_postfilt,'m');
                         h1(4) = plot(ax(ix),specs.f,specs.sig.CI_lower_postfilt,'m--');
                         plot(ax(ix),specs.f,specs.sig.CI_upper_postfilt,'m--')
                         h_plots = [h_plots,h1(3),h1(4)]; %#ok
                     elseif any(which2==4)
                         disp(plot_menus)
                         h1(3) = plot(ax(ix),specs.f,specs.sig.mean_postfilt,'m');
                         h_plots = [h_plots,h1(3)]; %#ok
                     end
                     if any(which2==5)
                         h1(5) = plot(ax(ix),specs.f,specs.background.mean_prefilt,'b');
                         h1(6) = plot(ax(ix),specs.f,specs.background.CI_lower_prefilt,'b--');
                         plot(ax(ix),specs.f,specs.background.CI_upper_prefilt,'b--')
                         h_plots = [h_plots,h1(5),h1(6)]; %#ok
                     elseif any(which2==6)
                         h1(5) = plot(ax(ix),specs.f,specs.background.mean_prefilt,'b');
                         h_plots = [h_plots,h1(5)]; %#ok
                     end
                     if any(which2==7)
                         h1(7) = plot(ax(ix),specs.f,specs.background.mean_postfilt,'c');
                         h1(8) = plot(ax(ix),specs.f,specs.background.CI_lower_postfilt,'c--');
                         plot(ax(ix),specs.f,specs.background.CI_upper_postfilt,'c--')
                         h_plots = [h_plots,h1(7),h1(8)]; %#ok
                     elseif any(which2==8)
                         h1(7) = plot(ax(ix),specs.f,specs.background.mean_postfilt,'c');
                         h_plots = [h_plots,h1(7)]; %#ok
                     end
                    legend(h_plots,legend_str(h1~=0),'location','northeast');
                    h_xlab = xlabel(ax(ix),'Frequency (Hz)'); 
                    pos_xlab = get(h_xlab,'Position');
                    pos_xlab(2) = pos_xlab(2)+.008; set(h_xlab,'Position',pos_xlab)
                    if spec_method==1 
                        title(ax(ix),'Average Power Spectrum')
                        ylabel(ax(ix),'Power (a.u.)')
                    elseif spec_method==2
                        title(ax(ix),'Average Power Spectral Density')
                        ylabel(ax(ix),'Power/Hz')
                    end
                end
                % Set Limits:
                if ~isnan(xlimits(1,ix)); xlim(ax(ix),xlimits(:,ix)'); end
                if ~isnan(ylimits(1,ix)); ylim(ax(ix),ylimits(:,ix)'); end
            end
            guidata(main_figure,guidata1);
        end
    end
    % Set Axis Limits:
    function set_limits(hObj,~,plot_num,which_ax)
        type = get(hObj,'Label');
        switch type
            case 'Autoscale'
                if strcmp(which_ax,'x')
                   xlimits(:,plot_num) = nan; 
                elseif strcmp(which_ax,'y')
                   ylimits(:,plot_num) = nan; 
                end
            case 'Specify'
                if strcmp(which_ax,'x')
                    prompt = {'Min: ','Max: '};
                    dlg_title = 'X-Limits'; num_lines = [1,15;1,15]; 
                    answer = inputdlg(prompt,dlg_title,num_lines);
                    if isempty(answer); disp('User cancelled action.'); return; end
                    xlimits(:,plot_num) = [str2double(answer{1});str2double(answer{2})];
                elseif strcmp(which_ax,'y')
                    prompt = {'Min: ','Max: '};
                    dlg_title = 'Y-Limits'; num_lines = [1,15;1,15]; 
                    answer = inputdlg(prompt,dlg_title,num_lines);
                    if isempty(answer); disp('User cancelled action.'); return; end
                    ylimits(:,plot_num) = [str2double(answer{1});str2double(answer{2})];
                end
        end
        update_plot
    end
    % Print Function:
    function print_spec(~,~,~)
        % Identify File Extension:
        ext_opts = {'.png',' ';'.tiff',' ';'.bmp',' '};
        [title1,path1] = uiputfile(ext_opts,'Specify filename:');
        figure_title = fullfile(path1,title1);
        [~,~,file_ext] = fileparts(figure_title);
        % Identify File Extension:
        exts = {'-dpng','-dtiff','-dbmp'};
        for ix = 1:3
            if ~isempty(strfind(exts{ix},file_ext(2:end)))
                ext_ind = ix; break;
            end
        end
        % Print Figure:
        set(h_power,'InvertHardcopy','off','PaperPositionMode','auto')
        print(h_power,exts{ext_ind},'-r300','-loose',figure_title)
        disp(['Printed figure: ',figure_title])
    end
    % Draw ROI for PowerSpec Button:
    function drawROI4Power(~,~,~)
        guidata1 = guidata(main_figure);
        instructions = sprintf(['Simply draw an ROI in Red (1), propagate through slices if desired, and close the neuroimage_editor.',...
            ' The ROI will be automatically saved, the average signal extracted, and the power spectrum calculated.']);
        h_mbox = msgbox(instructions,'ROI Tracing Instructions');
        uiwait(h_mbox,60)
        [h_GUI,last_draw] = neuroimage_editor_tbx('apply_header',guidata1.apply_header,...
            'background',guidata1.funct3D); %#ok
        uiwait(h_GUI,60) % wait until neuroimage_editor closes, or 60 sec
        last_draw = evalin('base','last_draw'); % pull last_draw in from 'base' workspace
        % Check if Detrending previously performed:
        if guidata1.use_linear
            detrend_type = 1;
        elseif guidata1.use_linear==0
            detrend_type = 2;
        elseif isempty(guidata1.use_linear)
            detrend_type = 0;
        end
        if isempty(guidata1.HighPass); guidata1.HighPass=0; end
        if isempty(guidata1.LowPass); guidata1.LowPass=0; end
        [~,~] = batch_power_spectrum(guidata1.funct, guidata1.Fs,'ROI',last_draw,...
            'spec_method',2,'power_type',1,'detrend',detrend_type,...
            'bandpass',[guidata1.HighPass,guidata1.LowPass],'anatomical',guidata1.funct3D);
        guidata(main_figure,guidata1);
    end
    % Calculate Spectrums:
    function specs = calc_specs(process_status)
        guidata1 = guidata(main_figure);
        if ~isempty(guidata1.sig_detrended) && all(size(guidata1.sig_detrended)==size(guidata1.roi_sig))
            prefilt_sig = guidata1.sig_detrended;
        else
            prefilt_sig = guidata1.roi_sig;
        end
        if ~isempty(guidata1.background_detrended) && all(size(guidata1.background_detrended)==size(guidata1.roi_background))
            prefilt_background = guidata1.background_detrended;
        else
            prefilt_background = guidata1.roi_background;
        end
        % Output Structure:
        specs = struct('f',[],'sig',struct('prefilt',[],'mean_prefilt',[],...
            'CI_lower_prefilt',[],'CI_upper_prefilt',[],'postfilt',[],...
            'mean_postfilt',[],'CI_lower_postfilt',[],'CI_upper_postfilt',[]),...
            'background',struct('prefilt',[],'mean_prefilt',[],...
            'CI_lower_prefilt',[],'CI_upper_prefilt',[],'postfilt',[],...
            'mean_postfilt',[],'CI_lower_postfilt',[],'CI_upper_postfilt',[]));
        % Pre-Filter Signal ROI Spectrum
        if spec_method==1
            if guidata1.numvox < 10000
                [specs.sig.prefilt,f] = pwelch(prefilt_sig,[],[],[],guidata1.Fs,'power');
            else
                hwait = waitbar(0,'Calculating power spectral density...');
                nchunk = ceil(guidata1.numvox/10000);
                chunks = 1:10000:guidata1.numvox; chunks(end+1) = guidata1.numvox+1;
                for ix = 1:nchunk
                    waitbar(ix/nchunk,hwait,'Calculating power spectral density...');
                    [specs.sig.prefilt(:,chunks(ix):chunks(ix+1)-1),f] = pwelch(prefilt_sig(:,chunks(ix):chunks(ix+1)-1),[],[],[],guidata1.Fs,'power');
                end
                close(hwait)
            end
        elseif spec_method==2
            if guidata1.numvox < 10000
                [specs.sig.prefilt,f] = pmtm(prefilt_sig,[],[],guidata1.Fs);
            else
                hwait = waitbar(0,'Calculating power spectral density...');
                nchunk = ceil(guidata1.numvox/10000);
                chunks = 1:10000:guidata1.numvox; chunks(end+1) = guidata1.numvox+1;
                for ix = 1:nchunk
                    waitbar(ix/nchunk,hwait,'Calculating power spectral density...');
                    [specs.sig.prefilt(:,chunks(ix):chunks(ix+1)-1),f] = pmtm(prefilt_sig(:,chunks(ix):chunks(ix+1)-1),[],[],guidata1.Fs);
                end
                close(hwait)
            end
        end
%         specs.sig.prefilt = specs.sig.prefilt./repmat(sum(specs.sig.prefilt),length(f),1); % Normalize [0,1]
        specs.f = f;
        sortamps = sort(specs.sig.prefilt,2);
        specs.sig.mean_prefilt = mean(sortamps,2);
        specs.sig.CI_lower_prefilt = sortamps(:,round(.05*size(sortamps,2)));
        specs.sig.CI_upper_prefilt = sortamps(:,round(.95*size(sortamps,2)));
        % Post-Filter Signal ROI Spectrum
        if process_status(1)
           if spec_method==1
                if guidata1.numvox < 10000
                    [specs.sig.postfilt,~] = pwelch(guidata1.roi_sig_filt,[],[],[],guidata1.Fs,'power');
                else
                    hwait = waitbar(0,'Calculating power spectral density (post-filtered)...');
                    for ix = 1:nchunk
                        waitbar(ix/nchunk,hwait,'Calculating power spectral density (post-filtered)...');
                        [specs.sig.postfilt(:,chunks(ix):chunks(ix+1)-1),~] = pwelch(guidata1.roi_sig_filt(:,chunks(ix):chunks(ix+1)-1),[],[],[],guidata1.Fs,'power');
                    end
                    close(hwait)
                end
            elseif spec_method==2
                if guidata1.numvox < 10000
                    [specs.sig.postfilt,~] = pmtm(guidata1.roi_sig_filt,[],[],guidata1.Fs);
                else
                    hwait = waitbar(0,'Calculating power spectral density (post-filtered)...');
                    for ix = 1:nchunk
                        waitbar(ix/nchunk,hwait,'Calculating power spectral density (post-filtered)...');
                        [specs.sig.postfilt(:,chunks(ix):chunks(ix+1)-1),~] = pmtm(guidata1.roi_sig_filt(:,chunks(ix):chunks(ix+1)-1),[],[],guidata1.Fs);
                    end
                    close(hwait)
                end
           end
%             specs.sig.postfilt = specs.sig.postfilt./repmat(sum(specs.sig.postfilt),length(f),1); % Normalize [0,1]
            sortamps = sort(specs.sig.postfilt,2);
            specs.sig.mean_postfilt = mean(sortamps,2);
            specs.sig.CI_lower_postfilt = sortamps(:,round(.05*size(sortamps,2)));
            specs.sig.CI_upper_postfilt = sortamps(:,round(.95*size(sortamps,2)));
        end
        % Pre-Filter Background ROI
        if process_status(2)
           if spec_method==1
                if guidata1.numvox < 10000
                    [specs.background.prefilt,~] = pwelch(prefilt_background,[],[],[],guidata1.Fs,'power');
                else
                    hwait = waitbar(0,'Calculating power spectral density (background, pre-filtered)...');
                    nchunk = ceil(guidata1.numvox/10000);
                    chunks = 1:10000:guidata1.numvox; chunks(end+1) = guidata1.numvox+1;
                    for ix = 1:nchunk
                        waitbar(ix/nchunk,hwait,'Calculating power spectral density (background, pre-filtered)...');
                        [specs.background.prefilt(:,chunks(ix):chunks(ix+1)-1),~] = pwelch(prefilt_background(:,chunks(ix):chunks(ix+1)-1),[],[],[],guidata1.Fs,'power');
                    end
                    close(hwait)
                end
            elseif spec_method==2
                if guidata1.numvox < 10000
                    [specs.background.prefilt,~] = pmtm(prefilt_background,[],[],guidata1.Fs);
                else
                    hwait = waitbar(0,'Calculating power spectral density (background, pre-filtered)...');
                    nchunk = ceil(guidata1.numvox/10000);
                    chunks = 1:10000:guidata1.numvox; chunks(end+1) = guidata1.numvox+1;
                    for ix = 1:nchunk
                        waitbar(ix/nchunk,hwait,'Calculating power spectral density (background, pre-filtered)...');
                        [specs.background.prefilt(:,chunks(ix):chunks(ix+1)-1),~] = pmtm(prefilt_background(:,chunks(ix):chunks(ix+1)-1),[],[],guidata1.Fs);
                    end
                    close(hwait)
                end
            end
%             specs.background.prefilt = specs.background.prefilt./repmat(sum(specs.background.prefilt),length(f),1); % Normalize [0,1]
            sortamps = sort(specs.background.prefilt,2);
            specs.background.mean_prefilt = mean(sortamps,2);
            specs.background.CI_lower_prefilt = sortamps(:,round(.05*size(sortamps,2)));
            specs.background.CI_upper_prefilt = sortamps(:,round(.95*size(sortamps,2)));
        end
        % Post-Filter Background ROI
        if process_status(3)
           if spec_method==1
                if guidata1.numvox < 10000
                    [specs.background.postfilt,~] = pwelch(guidata1.roi_background_filt,[],[],[],guidata1.Fs,'power');
                else
                    hwait = waitbar(0,'Calculating power spectral density (background, post-filtered)...');
                    nchunk = ceil(guidata1.numvox/10000);
                    chunks = 1:10000:guidata1.numvox; chunks(end+1) = guidata1.numvox+1;
                    for ix = 1:nchunk
                        waitbar(ix/nchunk,hwait,'Calculating power spectral density (background, post-filtered)...');
                        [specs.background.postfilt(:,chunks(ix):chunks(ix+1)-1),~] = pwelch(guidata1.roi_background_filt(:,chunks(ix):chunks(ix+1)-1),[],[],[],guidata1.Fs,'power');
                    end
                    close(hwait)
                end
            elseif spec_method==2
                if guidata1.numvox < 10000
                    [specs.background.postfilt,~] = pmtm(guidata1.roi_background_filt,[],[],guidata1.Fs);
                else
                    hwait = waitbar(0,'Calculating power spectral density (background, post-filtered)...');
                    for ix = 1:nchunk
                        waitbar(ix/nchunk,hwait,'Calculating power spectral density (background, post-filtered)...');
                        [specs.background.postfilt(:,chunks(ix):chunks(ix+1)-1),~] = pmtm(guidata1.roi_background_filt(:,chunks(ix):chunks(ix+1)-1),[],[],guidata1.Fs);
                    end
                    close(hwait)
                end
           end
%             specs.background.postfilt = specs.background.postfilt./repmat(sum(specs.background.postfilt),length(f),1); % Normalize [0,1]
            sortamps = sort(specs.background.postfilt,2);
            specs.background.mean_postfilt = mean(sortamps,2);
            specs.background.CI_lower_postfilt = sortamps(:,round(.05*size(sortamps,2)));
            specs.background.CI_upper_postfilt = sortamps(:,round(.95*size(sortamps,2)));
        end
    end % end calc_specs
end % end funct_view_spectrum

function funct_SNR_callback(~, ~, SNR_type, verbose, ~) 
    guidata1 = guidata(main_figure);
    % Delete display figure if present:
    if isgraphics(guidata1.h); delete(guidata1.h); end
    % Check if Funct Loaded:
    if isempty(guidata1.funct); errordlg('Please load a functional image.'); return; end
    % If no ROIs found, exit:
    if isempty(guidata1.ROI_ind); errordlg('No ROIs found. Please specify.'); return; end
    % Determine Type of SNR Requested:
    switch SNR_type
        case 1 % tSNR = mean voxel signal across time divided by SD across time
            % Calculate Stat:
            noise_SD = std(guidata1.roi_sig);
            mean_sig = mean(guidata1.roi_sig);
            % Avoid Inf's or extreme outliers, in case signal ROI included any zero voxels:
            if any(noise_SD<.01)
                ind = noise_SD<.01;
                noise_SD(ind) = nan;
                mean_sig(ind) = nan;
            end
%             guidata1.tSNR = mean_sig ./ noise_SD;
%             guidata1.funct_tSNR = nanmean(guidata1.tSNR);
            guidata1.tSNR = nanmean(guidata1.funct.img,4) ./ nanstd(guidata1.funct.img,0,4);
            
            % Set NaN's, Inf's, and out of range values (negatives) to 0:
            guidata1.tSNR(isinf(guidata1.tSNR)) = 0;
            guidata1.tSNR(isnan(guidata1.tSNR)) = 0;
            guidata1.tSNR(guidata1.tSNR<0) = 0;
            guidata1.tSNR(guidata1.tSNR>9999) = 0;
            % Calculate for Signal ROI:
            guidata1.tSNR_ROI = mean_sig ./ noise_SD;
            guidata1.funct_tSNR = nanmean(guidata1.tSNR_ROI);
            if verbose
                disp(['tSNR = ',num2str(guidata1.funct_tSNR)])
            end
        case 2 % tSBNR
            % Calculate Stat:
            noise_SD = std(guidata1.roi_background);
            mean_sig = mean(guidata1.roi_sig);
            % Avoid Inf's or extreme outliers, in case signal ROI included any zero voxels:
            if any(noise_SD<.0001)
                ind = noise_SD<.0001;
                noise_SD(ind) = nan;
            end
            noise_SD = nanmean(noise_SD);
            guidata1.tSBNR = nanmean(guidata1.funct.img,4) ./ noise_SD;
            % Set NaN's, Inf's, and other out of range values to 0;
            guidata1.tSBNR(isinf(guidata1.tSBNR)) = 0;
            guidata1.tSBNR(isnan(guidata1.tSBNR)) = 0;
            guidata1.tSBNR(guidata1.tSBNR<0) = 0;
            guidata1.tSBNR(guidata1.tSBNR>9999) = 0;
            % Calculate for Signal ROI
            guidata1.tSBNR_ROI = mean_sig ./ noise_SD;
            guidata1.funct_tSBNR = nanmean(guidata1.tSBNR_ROI);
            if verbose
                disp(['tSBNR = ',num2str(guidata1.funct_tSBNR)])
            end
        case 3 % SFNR
            if isempty(guidata1.funct); errordlg('First load a functional image.'); return; end
            if isempty(guidata1.roi_sig); errordlg('First specify ROIs.'); return; end
            % Provide Warning Regarding Wait-Time for this Process:
            waitbar_on = false;
            if isempty(guidata1.sig_detrended) || (~all(size(guidata1.sig_detrended)==size(guidata1.roi_sig)))
                if guidata1.numvox>20000 % override to auto use linear instead of 2nd order
                    if nargin==5
                        guidata1.use_linear=true;
                    else
                        [guidata1.use_linear, cancelled] = waitTimeUserInput(guidata1.numvox);
                        if cancelled; return; end
                    end
                else
                    guidata1.use_linear = false;
                end
                % Detrend:
                if guidata1.use_linear
                    disp('Calculating SFNR using linear detrending.')
                    guidata1.sig_detrended = detrend(guidata1.roi_sig); 
                else
                    disp('Calculating SFNR using quadratic detrending.')
                    guidata1.sig_detrended = guidata1.roi_sig; % initialize
                    x = (1:guidata1.n_4D)';
                    hwait = waitbar(0,'Detrending...'); waitbar_on = true;
                    for i = 1:guidata1.numvox
                        if rem(i,1000)==0; waitbar(i/guidata1.numvox,hwait); end
                        p = polyfit(x,guidata1.roi_sig(:,i),2);
                        predicted = polyval(p,x);
                        guidata1.sig_detrended(:,i) = guidata1.roi_sig(:,1)-predicted;
                    end
                end
            end
            noise_SD = std(guidata1.sig_detrended);
            mean_sig = mean(guidata1.roi_sig,1);
            % Avoid Inf's or extreme outliers, in case signal ROI included any zero voxels:
            if any(noise_SD<.001)
                ind = noise_SD<.001;
                noise_SD(ind) = nan;
                mean_sig(ind) = nan;
            end
            % Calculate Local:
            guidata1.SFNR = mean_sig ./ noise_SD;
            guidata1.funct_SFNR = nanmean(guidata1.SFNR);
            if waitbar_on; delete(hwait); end
            if verbose
                disp(['SFNR = ',num2str(guidata1.funct_SFNR)])
            end
        case 4 % SNR-Funct
            if isempty(guidata1.ROI_ind_back); errordlg('SNR-Funct requires a background ROI.'); return; end
             % Calculate Stat:
            noise_SD = std(guidata1.funct3D.img(guidata1.ROI_ind_back));
            denominator = noise_SD; % sqrt(2/(4-pi))*noise_SD;
            guidata1.funct_SNR_sd_local = guidata1.funct3D.img(guidata1.ROI_ind)./denominator;
            guidata1.funct_SNR_sd_all = guidata1.funct3D.img./denominator;
            guidata1.funct_SNR_sd = nanmean(guidata1.funct_SNR_sd_local); % guidata1.funct_mean_signal / denominator;
            if verbose; disp(['SNR-Functional = ',num2str(guidata1.funct_SNR_sd)]); end
        case 5 % Mean Power
            check_prev = [isempty(guidata1.Fs), isempty(guidata1.HighPass), isempty(guidata1.LowPass)];
            if any(check_prev)
                power_settings = questdlg('Use previous bandpass filter settings?',...
                    'Mean Power Settings','Yes','No','Cancel','Yes');
                switch power_settings
                    case 'Yes'
                        dont_use = false;
                    case 'No'
                        dont_use = true;
                    case 'Cancel'
                        return;
                end
            else
                dont_use = false;
            end
            if any(check_prev) || dont_use
                % Detect TR:
                try
                TR_detect = guidata1.funct.hdr.dime.pixdim(5);
                % Detect TR Units:
                switch bitand(guidata1.funct.hdr.dime.xyzt_units, 56)	
                    case 8
                      t_units = 's';
                    case 16
                      t_units = 'ms';
                    case 24
                      t_units = 'microseconds';	% microsecond
                    otherwise
                      t_units = 'unknown_units';
                end
                catch; TR_detect = [];
                end 
                % Prompt User: effective TR; lower bound (high-pass); upper-bound (low-pass)
                if ~isempty(TR_detect) && TR_detect ~= 0
                    prompt = {sprintf(['Detected TR = ',num2str(TR_detect),' ',num2str(t_units),...
                        '. Multiply by # repetitions per 3D image to obtain sampling time. \n\nSampling Time:']),...
                        'Min Frequency (Hz):','Max Frequency (Hz):'};
                else
                    prompt = {'Sampling Time (s):','Min Frequency (Hz):','Max Frequency (Hz):'};
                end
                dlg_title = 'Mean Power Settings'; num_lines = [1,40;1,40;1,40]; defaultans = {'2','.01','.1'};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                if isempty(answer)
                    disp('User cancelled mean power calculation.'); return;
                end
                guidata1.Fs = str2double(answer(1));
                guidata1.HighPass = str2double(answer(2));
                guidata1.LowPass = str2double(answer(3));
            end
            if ~isempty(guidata1.sig_detrended) && all(size(guidata1.sig_detrended)==size(guidata1.roi_sig))
                prefilt_sig = guidata1.sig_detrended;
            else
                prefilt_sig = guidata1.roi_sig;
            end
            if ~isempty(guidata1.background_detrended) && all(size(guidata1.background_detrended)==size(guidata1.roi_background))
                prefilt_background = guidata1.background_detrended;
            else
                prefilt_background = guidata1.roi_background;
            end
            if guidata1.numvox < 10000
                [signal_amp,f] = pwelch(prefilt_sig,[],[],[],guidata1.Fs,'power');
            else
                hwait = waitbar(0,'Calculating power spectral density using Welch''s periodogram...');
                nchunk = ceil(guidata1.numvox/10000);
                chunks = 1:10000:guidata1.numvox; chunks(end+1) = guidata1.numvox+1;
                for ix = 1:nchunk
                    waitbar(ix/nchunk,hwait);
                    [signal_amp(:,chunks(ix):chunks(ix+1)-1),f] = pwelch(prefilt_sig(:,chunks(ix):chunks(ix+1)-1),[],[],[],guidata1.Fs,'power');
                end
                close(hwait)
            end
            ind = ((f>=guidata1.HighPass) + (f<=guidata1.LowPass))==2;
            signal_amp = mean(signal_amp(ind,:),1);
            guidata1.mean_signal_amp = mean(signal_amp);
            if ~isempty(prefilt_background) 
                if size(prefilt_background,2) < 10000
                    [background_amp,~] = pwelch(prefilt_background,[],[],[],guidata1.Fs,'power');
                else
                    hwait = waitbar(0,'Calculating power spectral density (background)...');
                    nchunk = ceil(guidata1.numvox_background/10000);
                    chunks = 1:10000:guidata1.numvox_background; chunks(end+1) = guidata1.numvox_background+1;
                    for ix = 1:nchunk
                        waitbar(ix/nchunk,hwait);
                        [background_amp(:,chunks(ix):chunks(ix+1)-1),~] = pwelch(prefilt_background(:,chunks(ix):chunks(ix+1)-1),[],[],[],guidata1.Fs,'power');
                    end
                    close(hwait)
                end
                ind = ((f>=guidata1.HighPass) + (f<=guidata1.LowPass))==2;
                background_amp = nanmean(background_amp(ind,:),1);
                guidata1.mean_background_amp = mean(background_amp);
                background_present = true;
            else
                background_present = false;
            end
            if verbose
                if background_present
                    [~,p,~,stats] = ttest2(signal_amp,background_amp);
                    disp(['Mean power in range ',num2str(guidata1.HighPass),...
                        '-',num2str(guidata1.LowPass),' Hz: ',...
                        sprintf('\n%.04f',guidata1.mean_signal_amp),' (Signal ROI) ',...
                        sprintf('\n%.04f',guidata1.mean_background_amp),' (Background ROI)  ',...
                        sprintf('\nt(%g) = %.04g, p = %.3g',[stats.df,stats.tstat,p])])
                else
                    disp(['Mean power in range ',num2str(guidata1.HighPass),...
                        '-',num2str(guidata1.LowPass),' Hz: ',...
                        sprintf('%.02f',guidata1.mean_signal_amp),' (Signal ROI).']);
                end
            end
            % Save Output:
            guidata1.funct_power = signal_amp;
    end
    guidata(main_figure,guidata1) 
end

function [img] = reconstruct_image(stat_vec,ind,time_series)
    guidata1 = guidata(main_figure);
    if nargin<3
        img = guidata1.funct3D; img.img(:) = 0;
        img.img(ind) = stat_vec;
    elseif time_series
        img = guidata1.funct; img.img(:) = 0;
        if ~isempty(guidata1.roi_sig_filt)
            for i = 1:guidata1.numvox
                img.img(guidata1.ROI_x(i),guidata1.ROI_y(i),guidata1.ROI_z(i),:) = guidata1.roi_sig_filt(:,i);
            end
        elseif ~isempty(guidata1.sig_detrended)
            for i = 1:guidata1.numvox
                img.img(guidata1.ROI_x(i),guidata1.ROI_y(i),guidata1.ROI_z(i),:) = guidata1.sig_detrended(:,i);
            end
        end
    end
end

function funct_display_SNR(hObject, eventdata, SNR_type) %#ok
    guidata1 = guidata(main_figure);
%     if isgraphics(guidata1.h.figure); delete(guidata1.h.figure); end
    if isempty(guidata1.funct); errordlg('Please load a functional image.'); return; end
    if isempty(guidata1.ROI_ind); errordlg('Please specify ROIs.'); return; end
    switch SNR_type
        case 1 % tSNR
            if isempty(guidata1.tSNR); errordlg('ERROR: First, calculate tSNR'); return;end
            img = guidata1.funct3D; img.img = guidata1.tSNR;
            handles = nifti_studio('apply_header',guidata1.apply_header,...
                'background',img,...
                'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet','title','tSNR',...
                'background_caxis',[min(img.img(:)),.7*max(img.img(:))]);
            guidata1.h = handles.figure;
        case 2 % tSBNR (temporal Signal-to-Background Noise Ratio)
            if isempty(guidata1.tSBNR); errordlg('ERROR: First, calculate tSBNR'); return; end
            img = guidata1.funct3D; img.img = guidata1.tSBNR;
            handles = nifti_studio('apply_header',guidata1.apply_header,...
                'background',img,...
                'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet','title','tSBNR',...
                'background_caxis',[min(img.img(:)),.85*max(img.img(:))]);
            guidata1.h = handles.figure;
        case 3 % SFNR
            if isempty(guidata1.SFNR); errordlg('ERROR: First, calculate SFNR'); return; end
            % Check whether to display signal ROI only or Full SFNR Image:
            ans1 = 'Signal ROI only'; ans2 = 'Full Image (additional processing)';
            answer = questdlg('Display which?','SFNR Display',ans1,ans2,ans1);
            switch answer
                case ans1
                    funct_SFNR_local = reconstruct_image(guidata1.SFNR,guidata1.ROI_ind);
                    handles = nifti_studio('apply_header',guidata1.apply_header,...
                        'background',funct_SFNR_local,...
                        'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet','title','SFNR',...
                        'background_caxis',[min(funct_SFNR_local.img(:)),.85*max(funct_SFNR_local.img(:))]);
                    guidata1.h = handles.figure;
                case ans2
                    if isempty(guidata1.SFNR_all)
                        % Prepare all data:
                        if isempty(numvox_total)
                            numvox_total = numel(guidata1.funct3D.img);
                        end
                        if numvox_total>20000 % override to auto use linear instead of 2nd order
                            [guidata1.use_linear, cancelled] = waitTimeUserInput(guidata1.numvox);
                            if cancelled; return; end
                        else
                            guidata1.use_linear = false;
                        end
                        % Detrend:
                        disp('Preparing data...');
                        all_sig = zeros(guidata1.n_4D, numvox_total); all_detrended = all_sig;
                        [ROI_x_all,ROI_y_all,ROI_z_all] = ind2sub(guidata1.funct_dim,1:numvox_total);
                        for i = 1:numvox_total
                            all_sig(:,i) = guidata1.funct.img(ROI_x_all(i),...
                                ROI_y_all(i),ROI_z_all(i),:);
                        end
                        if guidata1.use_linear
                            disp('Calculating SFNR using linear detrending...')
                            all_detrended = detrend(all_sig); 
                        else
                            disp('Calculating SFNR using quadratic detrending...')
                            hwait = waitbar(0,'Detrending...'); 
                            x = (1:guidata1.n_4D)';
                            for i = 1:numvox_total
                                if rem(i,10000)==0; waitbar(i/numvox_total,hwait,'Detrending...'); end
                                p = polyfit(x,all_sig(:,i),2);
                                predicted = polyval(p,x);
                                all_detrended(:,i) = all_sig(:,1)-predicted;
                            end
                            delete(hwait);
                        end
                        noise_SD = std(all_detrended);
                        mean_sig = mean(all_sig,1);
                        % Avoid Inf's or extreme outliers, in case signal ROI included any zero voxels:
                        if any(noise_SD<.001)
                            ind = noise_SD<.001;
                            noise_SD(ind) = 0;
                            mean_sig(ind) = 0;
                        end
                        % Calculate Local:
                        guidata1.SFNR_all = mean_sig ./ noise_SD;
                        % Set NaN's, Inf's, and out of range values to 0:
                        guidata1.SFNR_all(isinf(guidata1.SFNR_all)) = 0;
                        guidata1.SFNR_all(isnan(guidata1.SFNR_all)) = 0;
                        guidata1.SFNR_all(guidata1.SFNR_all<0) = 0;
                        guidata1.SFNR_all(guidata1.SFNR_all>9999) = 0;
                        guidata(main_figure,guidata1) % SFNR_all must be saved within switch statement before call to neuroimage_editor
                    end
                    % Display Image:
                    funct_SFNR_all = reconstruct_image(guidata1.SFNR_all,1:numvox_total);
                    handles = nifti_studio('apply_header',guidata1.apply_header,...
                        'background',funct_SFNR_all,...
                        'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet','title','SFNR',...
                        'background_caxis',[min(funct_SFNR_all.img(:)),.7*max(funct_SFNR_all.img(:))]);
                    guidata1.h = handles.figure;
            end
        case 4 % SNR-Funct
            if isempty(guidata1.funct_SNR_sd); errordlg('ERROR: First, calculate SNR-Functional'); return; end
            img = guidata1.funct3D; img.img = guidata1.funct_SNR_sd_all;
            handles = nifti_studio('apply_header',guidata1.apply_header,...
                'background',img,...
                'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet','title','SNR-Funct',...
                'background_caxis',[min(img.img(:)),.85*max(img.img(:))]);%'background_caxis',[0,max(guidata1.funct_SNR_sd_local)*1.3]);
            guidata1.h = handles.figure;
        case 5 % Mean Power
            % Check that "Mean Power" calc has run:
            if isempty(guidata1.funct_power); errordlg('ERROR: First, calculate Mean Power'); return; end
            % Check whether to display signal ROI only or Full SFNR Image:
            ans1 = 'Signal ROI only'; ans2 = 'Full Image (additional processing)';
            answer = questdlg('Display which?','Power Display',ans1,ans2,ans1);
            switch answer
                case ans1
                    funct_power = reconstruct_image(guidata1.funct_power,guidata1.ROI_ind);
                    usePercentile = .99;
                    sortIntense = sort(funct_power.img(funct_power.img>0),'ascend');
                    maxColor = sortIntense(round(usePercentile*numel(funct_power.img(funct_power.img>0))));
                    handles = nifti_studio('apply_header',guidata1.apply_header,...
                        'background',funct_power,...
                        'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet','title','Power',...
                        'background_caxis',[min(funct_power.img(:)),maxColor]);
                    guidata1.h = handles.figure;
                case ans2
                    if isempty(guidata1.signal_amp_all)
                        % Check if all_sig data for all voxels has been collected:
                        if isempty(numvox_total)
                            numvox_total = numel(guidata1.funct3D.img);
                        end
                        if isempty(all_sig)
                            disp('Retrieving data...')
                            all_sig = zeros(guidata1.n_4D, numvox_total); all_detrended = all_sig;
                            [ROI_x_all,ROI_y_all,ROI_z_all] = ind2sub(guidata1.funct_dim,1:numvox_total);
                            for i = 1:numvox_total
                                all_sig(:,i) = guidata1.funct.img(ROI_x_all(i),...
                                    ROI_y_all(i),ROI_z_all(i),:);
                            end
                        end
                        if numvox_total < 10000
                            if isempty(all_detrended)
                                [signal_amp,f] = pwelch(all_sig,[],[],[],guidata1.Fs,'power');
                            else
                                [signal_amp,f] = pwelch(all_detrended,[],[],[],guidata1.Fs,'power');
                            end
                        else
                            hwait = waitbar(0,'Calculating power spectral density using Welch''s periodogram...');
                            nchunk = ceil(numvox_total/10000);
                            chunks = 1:10000:numvox_total; chunks(end+1) = numvox_total+1;
                            if isempty(all_detrended)
                                for ix = 1:nchunk
                                    waitbar(ix/nchunk,hwait,'Calculating power spectral density using Welch''s periodogram...');
                                    [signal_amp(:,chunks(ix):chunks(ix+1)-1),f] = pwelch(all_sig(:,chunks(ix):chunks(ix+1)-1),[],[],[],guidata1.Fs,'power');
                                end
                            else
                                for ix = 1:nchunk
                                    waitbar(ix/nchunk,hwait,'Calculating power spectral density using Welch''s periodogram...');
                                    [signal_amp(:,chunks(ix):chunks(ix+1)-1),f] = pwelch(all_detrended(:,chunks(ix):chunks(ix+1)-1),[],[],[],guidata1.Fs,'power');
                                end
                            end
                            close(hwait)
                        end
                        ind = ((f>=guidata1.HighPass) + (f<=guidata1.LowPass))==2;
                        guidata1.signal_amp_all = mean(signal_amp(ind,:),1);
                        guidata(main_figure,guidata1);
                    end
                    funct_power_all = reconstruct_image(guidata1.signal_amp_all,1:numvox_total);
                    usePercentile = .99;
                    sortIntense = sort(funct_power_all.img(funct_power_all.img>0),'ascend');
                    idx_use = round(usePercentile*numel(funct_power_all.img(funct_power_all.img>0)));
                    if idx_use<=0
                        handles = nifti_studio('apply_header',guidata1.apply_header,...
                            'background',funct_power_all,...
                            'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet','title','Power');
                    else
                        maxColor = sortIntense(idx_use);
                        handles = nifti_studio('apply_header',guidata1.apply_header,...
                            'background',funct_power_all,...
                            'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet','title','Power',...
                            'background_caxis',[min(funct_power_all.img(:)),maxColor]);
                    end
                    guidata1.h = handles.figure;
            end
    end
    guidata(main_figure,guidata1) 
end

function funct_report(~, ~, ~) 
    %% Select Output Dir and Recalculate if Necessary:
    guidata1 = guidata(main_figure);
    if isgraphics(guidata1.h); delete(guidata1.h); end
    if isempty(guidata1.funct); errordlg('ERROR: Functional image not loaded.'); return; end
    if ~isempty(guidata1.report_dir) && isdir(guidata1.report_dir) %#ok
        [guidata1.report_dir] = uigetdir(guidata1.report_dir,'Select Output Folder:');
    else
        [guidata1.report_dir] = uigetdir(guidata1.last_funct_dir,'Select Output Folder:');
    end
    if guidata1.report_dir==0; disp('Generate report canceled.'); return; end
    % Get Checkbox Options:
    guidata1.funct_report_options = [guidata1.funct_box1.Value,...
        guidata1.funct_box2.Value,guidata1.funct_box3.Value,...
        guidata1.funct_box4.Value,guidata1.funct_box5.Value,...
        guidata1.funct_box6.Value,guidata1.funct_box7.Value,...
        guidata1.funct_box8.Value,guidata1.funct_box9.Value];
    guidata(main_figure,guidata1) % re-save data before callback functions
    % Check if ROIs Present:
    if isempty(guidata1.ROI_ind)
        disp('WARNING: No ROIs specified. Automatically calculating ROIs based on Signal Intensity.')
        funct_ROI_spec([], [], 2) 
    end
    % Calculate stats if empty:
    if isempty(guidata1.funct_tSNR)
        funct_SNR_callback([], [], 1, 0); % tSNR (4th input = verbose)
    end
    if isempty(guidata1.funct_tSBNR)
        funct_SNR_callback([], [], 2, 0, 1); % SFNR, auto-use linear detrend (5th input)
    end
    if isempty(guidata1.funct_SFNR)
        funct_SNR_callback([], [], 3, 0, 1); % SFNR, auto-use linear detrend (5th input)
    end
    if isempty(guidata1.funct_SNR_sd) && ~isempty(guidata1.ROI_ind_back)
        funct_SNR_callback([], [], 4, 0);
    end
    % Band-pass if necessary:
    if ~isempty(guidata1.Fs) % band-pass has been performed; use prev settings
       bandpass_present = true;
%        funct_bandpass([],[],[],guidata1.Fs,guidata1.HighPass,guidata1.LowPass) % re-run (need to add input to print spectrum)
       funct_SNR_callback([], [], 5, 0);
    else
        bandpass_present = false;
    end
% ############################# Outputs ###################################
    %% WRITE Stats.txt:
    if guidata1.funct_report_options(1) % Stats.txt
        try [~,fname1,~] = fileparts(guidata1.funct.fileprefix);
        catch; fname1 = 'Functional';
        end
        fname2 = ['Stats_',fname1,'.txt'];
        % Write Text File:
        fileID = fopen(fullfile(guidata1.report_dir,fname2),'w');
        fprintf(fileID,'%s\r\n',['Functional Name:  ',fname1]);
        fprintf(fileID,'%s\r\n',['ROI Type:  ',guidata1.funct_ROI_type]);
        fprintf(fileID,'%s\r\n',['Signal ROI # Voxels:  ',num2str(length(guidata1.ROI_ind))]);
        fprintf(fileID,'%s\r\n',['Background ROI # Voxels:  ',num2str(length(guidata1.ROI_ind_back))]);
        fprintf(fileID,'%s\r\n',' ');
        fprintf(fileID,'%s\r\n',['tSNR = ',num2str(guidata1.funct_tSNR)]);
        fprintf(fileID,'%s\r\n',['tSBNR = ',num2str(guidata1.funct_tSBNR)]);
        if guidata1.use_linear
            fprintf(fileID,'%s\r\n',['SFNR (linear detrend) = ',num2str(guidata1.funct_SFNR)]);
        else
            fprintf(fileID,'%s\r\n',['SFNR (quadratic detrend) = ',num2str(guidata1.funct_SFNR)]);
        end
        fprintf(fileID,'%s\r\n',['SNR-Funct = ',num2str(guidata1.funct_SNR_sd)]);
        fprintf(fileID,'%s\r\n',' ');
        if bandpass_present
            fprintf(fileID,'%s\r\n','Band-pass filter settings: ');
            fprintf(fileID,'%s\r\n',['Sampling Frequency (Fs) = ',num2str(guidata1.Fs),...
                '; Frequency Range: ',num2str(guidata1.HighPass),'-',num2str(guidata1.LowPass),' Hz']);
            fprintf(fileID,'%s\r\n',' ');
            if ~isempty(guidata1.mean_background_amp)
                fprintf(fileID,'%s\r\n',['Mean Power in Freq Range: ',...
                    num2str(guidata1.mean_signal_amp),' (Signal ROI); ',...
                    num2str(guidata1.mean_background_amp),' (Background ROI)']);
                fprintf(fileID,'%s\r\n',['Power Ratio (Signal ROI/Background ROI): ',...
                    num2str(guidata1.mean_signal_amp/guidata1.mean_background_amp)]);
            else
                fprintf(fileID,'%s\r\n',['Mean Power in Freq Range: ',...
                    num2str(guidata1.mean_signal_amp),' (Signal ROI)']);
            end
            % Add correlation matrix of metrics:
            if ~isempty(guidata1.ROI_ind_back)
                mat2corr = [guidata1.tSNR_ROI',guidata1.tSBNR_ROI',guidata1.SFNR',guidata1.funct_SNR_sd_local,guidata1.funct_power'];
                mat2corr = mat2corr(sum(isnan(mat2corr),2)==0,:);
                r = corr(mat2corr);
                fprintf(fileID,'%s\r\n',' ');
                fprintf(fileID,'%s\r\n','Correlation Among Metrics within Signal ROI:');
                fprintf(fileID,'%s\r\n','r = ');
                fprintf(fileID,'%s\r\n',sprintf('           tSNR, tSBNR, SFNR, SNR-Funct, Power'));
                fprintf(fileID,'%s\r\n',sprintf('tSNR:      %.3g,    %.3g,  %.3g,  %.3g,  %.3g',r(1,:)));
                fprintf(fileID,'%s\r\n',sprintf('tSBNR:     %.3g,  %.3g,    %.3g,  %.3g,  %.3g',r(2,:)));
                fprintf(fileID,'%s\r\n',sprintf('SFNR:      %.3g,  %.3g,  %.3g,    %.3g,  %.3g',r(3,:)));
                fprintf(fileID,'%s\r\n',sprintf('SNR-Funct: %.3g,  %.3g,  %.3g,  %.3g,    %.3g',r(4,:)));
                fprintf(fileID,'%s\r\n',sprintf('Power:     %.3g,  %.3g,  %.3g,  %.3g,  %.3g',r(5,:)));
            end
        else
            fprintf(fileID,'%s\r\n','No band-pass filtration performed.');
            % Add correlation matrix of metrics:
            mat2corr = [guidata1.tSNR_ROI',guidata1.tSBNR_ROI',guidata1.SFNR',guidata1.funct_SNR_sd_local];
            mat2corr = mat2corr(sum(isnan(mat2corr),2)==0,:);
            r = corr(mat2corr);
            fprintf(fileID,'%s\r\n',' ');
            fprintf(fileID,'%s\r\n','Correlation Among Metrics within Signal ROI:');
            fprintf(fileID,'%s\r\n','r = ');
            fprintf(fileID,'%s\r\n',sprintf('           tSNR, tSBNR, SFNR, SNR-Funct'));
            fprintf(fileID,'%s\r\n',sprintf('tSNR:      %.3g,    %.3g,  %.3g,  %.3g',r(1,:)));
            fprintf(fileID,'%s\r\n',sprintf('tSBNR:     %.3g,  %.3g,    %.3g,  %.3g',r(2,:)));
            fprintf(fileID,'%s\r\n',sprintf('SFNR:      %.3g,  %.3g,  %.3g,    %.3g',r(3,:)));
            fprintf(fileID,'%s\r\n',sprintf('SNR-Funct: %.3g,  %.3g,  %.3g,  %.3g  ',r(4,:)));
        end
        fclose(fileID);
    end
%% 2D Img Mosaics and Save 3D Images:

    % Determine slices and dimension for mosaic:
    if guidata1.funct_report_options(2) 
        dim3 = size(guidata1.funct3D.img,3);
        if (dim3>=15) 
            slices = round(linspace(3,dim3,15));
            axes_dim = [3,5];
        elseif (dim3>=10) 
            slices = round(linspace(3,dim3,9));
            axes_dim = [2,5];
        elseif (dim3>=5) 
            slices = round(linspace(3,dim3,5));
            axes_dim = [1,5];
        else
            slices = dim3;
            axes_dim = [1,slices+1];
        end
    end
    
    % ROIs:
    if guidata1.funct_report_options(2) || guidata1.funct_report_options(4)
        % Reconstruct 3D Image (saves RAM for ROIs, but not whole images):
        funct_ROI = reconstruct_image(ones(1,guidata1.numvox),guidata1.ROI_ind);
        if ~isempty(guidata1.ROI_ind_back)
            funct_ROI.img(guidata1.ROI_ind_back) = 2;
        end
        % 2D Mosaic:
        if guidata1.funct_report_options(2)
%             handles = nifti_studio('apply_header',guidata1.apply_header,'background',funct_ROI,...
%                 'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet',...
%                 'print',fullfile(guidata1.report_dir,[fname1,'_ROIs_mosaic.png']));
            [handles] = neuroimage_editor_mosaic('background',guidata1.funct3D,'overlay',funct_ROI,'slices',slices,...
                'colormap','jet','title','Signal and Background ROIs','slice_locator',0,...
                'axes_dim',axes_dim,'print',fullfile(guidata1.report_dir,[fname1,'_ROIs_mosaic.png']));
            close(handles.figure); 
        end
        % 3D Save:
        if guidata1.funct_report_options(4)
            save_img(guidata1.apply_header,funct_ROI,fullfile(guidata1.report_dir,[fname1,'_ROIs']))
        end
    end
    % tSNR:
    if guidata1.funct_report_options(2) || guidata1.funct_report_options(5)
        % Reconstruct:
        disp_tSNR = guidata1.funct3D; disp_tSNR.img = guidata1.tSNR;
        % 2D Mosaic:
        if guidata1.funct_report_options(2)
%             handles = nifti_studio('apply_header',guidata1.apply_header,'background',disp_tSNR,...
%                 'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet',...
%                 'print',fullfile(guidata1.report_dir,[fname1,'_tSNR_mosaic.png']));
            [handles] = neuroimage_editor_mosaic('background',guidata1.funct3D,'overlay',disp_tSNR,'slices',slices,...
                'colormap','jet','slice_locator',0,'axes_dim',axes_dim,'colorbar',1,...
                'print',fullfile(guidata1.report_dir,[fname1,'_tSNR_mosaic.png']));
            close(handles.figure);
        end
        % Save 3D Image:
        if guidata1.funct_report_options(5)
            save_img(guidata1.apply_header,disp_tSNR,fullfile(guidata1.report_dir,[fname1,'_tSNR']))
        end
        clear disp_tSNR
    end
    % tSBNR
    if guidata1.funct_report_options(2) || guidata1.funct_report_options(6)
        % Reconstruct:
        disp_tSBNR = guidata1.funct3D; disp_tSBNR.img = guidata1.tSBNR;
        % 2D Mosaic:
        if guidata1.funct_report_options(2)
%             handles = nifti_studio('apply_header',guidata1.apply_header,'background',disp_tSBNR,...
%                 'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet',...
%                 'print',fullfile(guidata1.report_dir,[fname1,'_tSBNR_mosaic.png']));
            [handles] = neuroimage_editor_mosaic('background',guidata1.funct3D,'overlay',disp_tSBNR,'slices',slices,...
                'colormap','jet','slice_locator',0,'axes_dim',axes_dim,'colorbar',1,...
                'print',fullfile(guidata1.report_dir,[fname1,'_tSBNR_mosaic.png']));
            close(handles.figure); 
        end
        % Save 3D Image:
        if guidata1.funct_report_options(6)
            save_img(guidata1.apply_header,disp_tSBNR,fullfile(guidata1.report_dir,[fname1,'_tSBNR']))
        end
        clear disp_tSBNR
    end
    % SFNR 
    if guidata1.funct_report_options(2) || guidata1.funct_report_options(7)
        % Reconstruct:
        if ~isempty(guidata1.SFNR_all)
            disp_SFNR = reconstruct_image(guidata1.SFNR_all,1:numvox_total);
        else
            disp_SFNR = reconstruct_image(guidata1.SFNR,guidata1.ROI_ind);
        end
         % 2D Mosaic:
        if guidata1.funct_report_options(2)
%             handles = nifti_studio('apply_header',guidata1.apply_header,'background',disp_SFNR,...
%                 'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet',...
%                 'print',fullfile(guidata1.report_dir,[fname1,'_SFNR_mosaic.png']));
            [handles] = neuroimage_editor_mosaic('background',guidata1.funct3D,'overlay',disp_SFNR,'slices',slices,...
                'colormap','jet','slice_locator',0,'axes_dim',axes_dim,'colorbar',1,...
                'print',fullfile(guidata1.report_dir,[fname1,'_SFNR_mosaic.png']));
            close(handles.figure);
        end
        if guidata1.funct_report_options(7)
            save_img(guidata1.apply_header,disp_SFNR,fullfile(guidata1.report_dir,[fname1,'_SFNR']))
        end
        clear disp_SFNR
    end
    % SNR-Funct 
    if ~isempty(guidata1.ROI_ind_back) &&  (guidata1.funct_report_options(2) || guidata1.funct_report_options(8))
        % Reconstruct:
        disp_SNR_funct = guidata1.funct3D; disp_SNR_funct.img = guidata1.funct_SNR_sd_all;
        % 2D Mosaic:
        if guidata1.funct_report_options(2)
%             handles = nifti_studio('apply_header',guidata1.apply_header,'background',disp_SNR_funct,...
%                 'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet',...
%                 'print',fullfile(guidata1.report_dir,[fname1,'_SNR_mosaic.png']));
            [handles] = neuroimage_editor_mosaic('background',guidata1.funct3D,'overlay',disp_SNR_funct,'slices',slices,...
                'colormap','jet','slice_locator',0,'axes_dim',axes_dim,'colorbar',1,...
                'print',fullfile(guidata1.report_dir,[fname1,'_SFNR_mosaic.png']));
            close(handles.figure); 
        end
        if guidata1.funct_report_options(8)
            save_img(guidata1.apply_header,disp_SNR_funct,fullfile(guidata1.report_dir,[fname1,'_funct_SNR']))
        end
        clear disp_SNR_funct
    end
    % Mean Power 
    if guidata1.funct_report_options(2) || guidata1.funct_report_options(9)
        % Reconstruct:
        if ~isempty(guidata1.signal_amp_all)
            disp_power = reconstruct_image(guidata1.signal_amp_all,1:numvox_total);
        elseif ~isempty(guidata1.funct_power)
            disp_power = reconstruct_image(guidata1.funct_power,guidata1.ROI_ind);
        else
            disp_power = [];
        end
        % 2D Mosaic:
        if guidata1.funct_report_options(2)
%             handles = nifti_studio('apply_header',guidata1.apply_header,'background',disp_power,...
%                 'colorbar_on',1,'title_on',1,'axis_tick_on',1,'colormap','jet',...
%                 'print',fullfile(guidata1.report_dir,[fname1,'_Power_mosaic.png']));
              [handles] = neuroimage_editor_mosaic('background',guidata1.funct3D,...
                  'overlay',disp_power,'slices',slices,'colormap','jet',...
                  'slice_locator',0,'axes_dim',axes_dim,'colorbar',1,...
                  'print',fullfile(guidata1.report_dir,[fname1,'_Power_mosaic.png']));
            close(handles.figure); 
        end
        % Save 3D Image:
        if guidata1.funct_report_options(9)
            save_img(guidata1.apply_header,disp_power,fullfile(guidata1.report_dir,[fname1,'_power']))
        end
        clear disp_power
    end
    % Save Preprocessed Time Series:
    if guidata1.funct_report_options(3) % preprocessed
        funct_preprocess = reconstruct_image([],[],1);
        if guidata1.use_linear; detrended_type = 'linear_detrend'; 
        else
            detrended_type = 'quadratic_detrend';
        end
        if ~isempty(guidata1.sig_detrended) && ~isempty(guidata1.roi_sig_filt)
            suffix = ['_',detrended_type,'_bandpass'];
        elseif ~isempty(guidata1.sig_detrended) && isempty(guidata1.roi_sig_filt)
            suffix = ['_',detrended_type,];
        elseif ~isempty(guidata1.sig_detrended) && isempty(guidata1.roi_sig_filt)   
            suffix = '_bandpass';
        end
        save_img(guidata1.apply_header,funct_preprocess,fullfile(guidata1.report_dir,[fname1,'preprocess_',suffix]))
    end
    disp('Report Finished.')
    guidata(main_figure,guidata1) 
end

function [use_linear, cancelled] = waitTimeUserInput(numvox)
    cancelled = false;
    wait_time = questdlg(['Quadratic detrending will take approximately ',...
        sprintf('%.02f',numvox/50000*36.12/60),' minutes. Continue or ',...
        'use linear detrending (nearly instantaneous)?'],'Detrend Wait-Time',...
        'Linear','Quadratic','Cancel','Linear');
    switch wait_time
        case 'Linear'
            use_linear=true;
        case 'Quadratic'
            use_linear=false;
        case 'Cancel'
            cancelled = true; use_linear = true;
            disp('User cancelled operation.');
    end
end

function functionals_closereq(~, ~, ~)
    guidata1 = guidata(main_figure);
    guidata1.functionals_button.Value=0;
    if isgraphics(guidata1.funct_figure)
        delete(guidata1.funct_figure)
    end
    if isgraphics(guidata1.h); close(guidata1.h); end
    guidata(main_figure,guidata1)
end

end % End: MRIqual.m