function [handles,last_draw,drawing_idx] = neuroimage_editor_tbx(apply_header,varargin)
% NEUROIMAGE_EDITOR  % A simple GUI for navigating, visualizing, and
% editing 3D images (file types: .nii, .nii.gz, .img/.hdr)
%
% Author: Elliot Layden, The University of Chicago, 2017-2019
% 
% Usage: 
% To begin, simply type "neuroimage_editor" into the command line
% (specify neuroimage_editor(1) if you want to apply header info (affine)
% or specify neuroimage_editor(0) if you do not want to apply this info).
% If neither is specified, a dialog box will pop-up asking you to select
% one or the other option. Alternatively, simply right-click 
% "neuroimage_editor.m" -> Run. Next, a file selection menu will be 
% displayed; select your desired 3D image file, and begin viewing/editing.
% 
% Optional Output:
%   'handles',          Handles structure for GUI figure, axes, etc.
%   'last_draw',        If declared as an output, a variable named
%                       'last_draw' will be assigned in the 'base'
%                       workspace. 'last_draw' is a #voxel x 3 matrix of
%                       indices for the last drawing or 'propagate drawing'
%                       created (i.e., [ix,iy,iz]); written to 'base' upon
%                       closing of the GUI
%   'drawing_idx',      a 3D matrix of size equal to the background image
%                       with integers denoting the current drawing; written 
%                       to 'base' upon closing of the GUI
% 
% Optional Inputs (Name-Value Pair Arguments):
%   (Note: these can also be loaded within the GUI)
%   'background',       filename or path-filename to a background image to
%                       be loaded; can also be loaded image structure 
%   'overlay1',         filename or path-filename to an image to be loaded
%                       as overlay #1; can also be loaded image structure
%   'overlay2',         filename or path-filename to an image to be loaded
%                       as overlay #2; can also be loaded image structure
%   'colorbar',         customization: true/false
%   'title_on',         customization: true/false (refers to slice numbers)
%   'axis_tick_on',     customization: true/false
%   'colormap',         colormap for main figure (background), specified as
%                       string, e.g. 'colormap','jet'
%   'print',            if valid filepath specified, will print the current
%                       loaded image after all other name-value pairs are 
%                       applied
%   'axes',             an axes graphics object (handle); used to embed a 
%                       neuroimage_editor window within an already present
%                       figure/axes
%   'background_caxis', 1x2 vector [cmin, cmax] to scale colors of
%                       background image
%   'overlay_caxis',    1x2 vector [cmin, cmax] to scale colors of
%                       overlay image
%   'title',            String specifying custom figure window title;
%                       overides default filename title
% 
% [Note: if you encounter error messages when loading a file that includes 
% "non-orthogonal shearing", this means that applying the affine to the 
% image is introducing distortion beyond the tolerance range of NIFTI Tools. 
% In this case, try loading the image without applying the affine/header 
% info (neuroimage_editor(0); equivalent to using load_untouch_nii.m in 
% NIFTI Tools)]
%
% General Features:
% 'apply_header' = true/false (apply affine contained in header?)
% -use up and down arrow keys to navigate slices (z-dimension)
% -use the mouse scroll wheel to zoom in/out
% -click and drag mouse to draw shapes; for instance, drawing an arc will 
%     cause the interior of the arc to become filled with color; any closed 
%     figure will also be filled with color (particularly useful for erasing 
%     undesirable parts of images like artifacts, or for creating ROIs)
%     -can draw to edit underlying image data, to draw new ROIs which can
%     be saved as a separate file from the main image, or to edit loaded
%     overlay images (such as a set of ROIs saved from an external program)
% -specify draw color based on image's intensity units (default = 0, i.e. 
% erase) for the underlying data or for an overlay, or based on a set of 
% pre-defined colors for ROI drawing
% -the transparency of overlays or ROI drawings can be easily adjusted
% -auto-detects screen resolution for positioning of figure window
% -figure can be resized if desired, with objects positioned in normalized 
%     units
% 
% Buttons/Actions (detailed):
% -Click once: draws current color at voxel (default = 0)
% -Click and Drag: upon release, draws line or fills in shape with color
% -Mouse Scroll Wheel (toward computer): zooms in (10 zoom settings: 100%, 
%     90%, 80%,...10% of slice area in display window) at current mouse 
%     location; (away from computer): zooming out at  different mouse 
%     location will also adjust the zoom location

% Menu's and Options:
% "FILE" Menu:
%   -"Open" (hotkey: 'o'): select new file to open; be careful to save first
%   -"Save Image" (hotkey: 's'): save current image and any edits with the 
%       same image filename as was loaded
%   -"Save Image As": specify new file name to save image as
%   -"Save Drawing": save current ROI drawing using predefined colors
%   -"Save Drawing As": specify new filename with which to save ROI drawing
%   -"Print Current 2D Image": Print current image display (including any 
%       overlays or drawings present) as a .png, .tiff, or .bmp in 300 dpi 
%       (can change dpi setting at line 343 "-r<dpi>"; auto-crops to 
%       image-only)
%   -"Exit" (hotkey: 'e'): Exits the program, first prompting user to
%       verify the desire to exit
% "EDIT" Menu:
%   -"Undo" (hotkey: 'z'): undo last action (whether edit to underying
%       data, drawing, or overlay); stores 2 actions back in RAM
%   -"Redo" (hotkey: 'r'): redo undone action (up to 2 forward)
%   -"Go to Slice..." (hotkey: 'g'): navigate to specified slice
% "DISPLAY" Menu:
%   -"Orientation": change orientation between coronal, sagittal, or axial
%   -"Colormap": change colormap of the underyling data/background image
%       -note: overlay colormaps can be changed in the Overlay menu, see
%       below
%   -"Autoscale Slice Color": if yes, autoscales color axis of each slice
%       (this can slightly reduce performance when navigating quickly
%       through slices)
%       -if "No", can either use the max and min color/intensity of the
%       whole 3D image to scale each slice, or can specify custom color
%       axis limits
%   -"Slice # (On/Off)": turn on/off title which displays slice #
%   -"Colorbar (On/Off)": turn on/off colorbar
%   -"Axis Tick (On/Off)": turn on/off axis ticks
%       -note: if all of the above 3 are turned off, the image display is
%       maximized in the figure window
% "DRAW" Menu:
%   -"Edit Underlying Data": this menu allows user to edit the background
%       image's underlying data
%       -"Select Draw Color": specify an image intensity to use for drawing 
%           (default: 0)
%       -"Crop Border" (hotkey: 'c'): specify a # of voxels to crop from
%           each edge of current slice
%       -"Fill Slice": fills entire slice with current draw color intensity
%   -"Edit Drawing": this menu allows user to create an ROI drawing
%       overlayed on top of the background image
%       -"Select Draw Color": select from a set of predefined colors (note
%       	that the numbers in parentheses are the values which will be
%       	written to image voxels if the drawing is saved); if "Disable" 
%           is selected, then drawing control will return to the background
%           image/underlying data
%       -"Select Draw Opacity": select the opacity/transparency of the
%           drawing overlay; this can be chosen separately for any loaded
%           overlays (see Overlay menu below)
%       -"Clear Drawing": this clears any current drawing on the 3D image
%           (but not any edits made to the background image or to any 
%           loaded overlays); user is prompted to confirm action
% "OVERLAY" Menu:
%   -"Add/Edit Overlays": this creates a new GUI window specific to actions
%       involving any external overlay images that are loaded
%       -Options include, "Open", "Close", "Save", change colormap, change
%           opacity, change color axis limits (color range), toggle
%           colorbar, and "Edit"
%       -Note: toggling "Colorbar" on will replace the background image
%           colorbar with one that corresponds to the loaded overlay;
%           "Colorbar (On/Off)" in the "Display" menu must be toggled on
%           for this to have an effect
%       -Note 2: if "Edit" is toggled on, user is prompted to enter an
%           intensity value for drawing which corresponds to the overlay 
%           image (and not the background image); if "Edit" is toggled off,
%           drawing control will be returned to the previous mode (either
%           background/underlying data editing or ROI drawing)
% 
% Dependencies: NIFTI_tools (Shen, 2005) must be on Matlab's path (download
% link: http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
% Note that the necessary functions have been included in the 
% neuroimage_editor download folder, so no further action is required.
% 
% Note: colorbars and overlays will not function properly for Matlab
% versions prior to 2014b, due to the major graphics update which came in
% 2014b
% 
%% Load Image and Initialize to Display a Middle Slice:

figure_color = [.2,.2,.2]; % zeros(1,3);
axes_color = ones(1,3);

% Identify Function Path and Add Helper Scripts:
script_fullpath = mfilename('fullpath');
[script_path,~,~] = fileparts(script_fullpath);
addpath(genpath(script_path))

% Determine whether to apply header:
apply_header_dlg = false;
if nargin > 0 
    try
        if apply_header
            untouch_nii = false;
        elseif ~apply_header
            untouch_nii = true;
        end
    catch
        apply_header_dlg = true;
    end
else
    apply_header_dlg = true;
end
if apply_header_dlg
    hdr_answer = questdlg('Apply header information to image?','Header Action');
    if strcmp(hdr_answer,'Yes')
        untouch_nii = false;
    elseif strcmp(hdr_answer,'No')
        untouch_nii = true;
    else
        handles = [];
        return;
    end
end

% Determine Inputs:
inputs = varargin;
parsed_inputs = struct('background',[],'overlay1',[],'overlay2',[],...
    'colorbar',1,'title_on',1,'axis_tick_on',1,'colormap','gray','print',[],...
    'axes',[],'background_caxis',[],'overlay_caxis',[],'title','',...
    'overlay_alpha',.9);
poss_input = {'background','overlay1','overlay2','colorbar','title_on',...
    'axis_tick_on','colormap','print','axes','background_caxis','overlay_caxis',...
    'title','overlay_alpha'};
input_ind = zeros(1,length(poss_input));
for ixx = 1:length(poss_input)
    ixx1 = find(strcmp(poss_input{ixx},inputs));
    if ~isempty(ixx1)
        input_ind(ixx) = ixx1;
        input1 = inputs{input_ind(ixx)+1};
        if ischar(input1) && exist(input1,'file')==2
            parsed_inputs.(poss_input{ixx}) = input1;
        elseif isstruct(input1)
            parsed_inputs.(poss_input{ixx}) = input1;
        elseif ixx>3 && islogical(input1) || isnumeric(input1)
            parsed_inputs.(poss_input{ixx}) = input1;
        elseif ixx==7 && ischar(input1)
            parsed_inputs.(poss_input{ixx}) = input1;
        elseif ixx==8
            parsed_inputs.(poss_input{ixx}) = input1;
        elseif ixx==9 && isgraphics(input1,'axes')
            parsed_inputs.(poss_input{ixx}) = input1;
        elseif ixx>9 && ixx<11 && isvector(input1)
            parsed_inputs.(poss_input{ixx}) = input1;
        elseif ixx==12 && ischar(input1)
            parsed_inputs.(poss_input{ixx}) = input1;
        elseif ixx==13 && isnumeric(input1)
            parsed_inputs.(poss_input{ixx}) = input1;
        else
            error(['Invalid input for ''',inputs{input_ind(ixx)},...
                '''. Please specify filepath or image structure.'])
        end
    end
end

% If missing background input, clear others:
if isempty(parsed_inputs.background)
    parsed_inputs.overlay1 = [];
    parsed_inputs.overlay2 = [];
    parsed_inputs.print = [];
    parsed_inputs.axes = [];
end

% Initialize Handles & Options:
handles = struct('figure',99.999,'background_axes',99.999,...
    'background_img',99.999,'drawing_axes',...
    99.999,'drawing_img',99.999,'overlay1_axes',99.999,'overlay1_img',...
    99.999,'overlay2_axes',99.999,'overlay2_img',99.999);

% Initialize Figure (auto-detect screen size):
set(0,'units','pixels'); screen_res = get(0,'ScreenSize');
figure_pos = [.25*screen_res(3), .065*screen_res(4), .5*screen_res(3), ...
    .86*screen_res(4)];
if any(figure_pos<=0) % avoid ScreenSize errors
    figure_pos = [342,50,683,660]; 
    screen_res = [1,1,(figure_pos(3:4)./[.5,.86])];
end
if isempty(parsed_inputs.axes) % if no axes supplied
    handles.figure = figure('Position',figure_pos,'MenuBar','none',...
        'Name','neuroimage_editor','NumberTitle','off','Color',figure_color,...
        'Visible','off','doublebuffer','on','Interruptible','off');
else
    handles.figure = get(parsed_inputs.axes,'Parent');
end
    
% Intialize GUI Data:
guidata1 = guidata(handles.figure);
guidata1.figure_color = figure_color; % [0,0,0]; % [.8,.88,.98];
guidata1.axes_color = axes_color; % ones(1,3)
last_draw = [];
drawing_idx = [];
guidata1.main_colormap_selection = parsed_inputs.colormap; % default colormap
guidata1.title_on = parsed_inputs.title_on;
guidata1.colorbar = parsed_inputs.colorbar; % default: On (option for any colorbar or no colorbars)
guidata1.axis_tick_on = parsed_inputs.axis_tick_on;
if isempty(parsed_inputs.background_caxis)
    guidata1.autoscale = true; % default setting
else
    guidata1.autoscale = false;
end
guidata1.nargout = nargout;
if ~isempty(parsed_inputs.axes)
    guidata1.ax_pos_input = get(parsed_inputs.axes,'Position');
end
guidata1.script_path = script_path;
guidata1.untouch_nii = untouch_nii;
guidata1.screen_res = screen_res;
guidata1.figure_pos = figure_pos;
guidata1.fig_height = guidata1.figure_pos(4);
guidata1.screen_mid_x = .5*guidata1.screen_res(3);
guidata1.parsed_inputs = parsed_inputs;
guidata1.poss_input = poss_input;
if isempty(parsed_inputs.axes)
    handles.background_axes = 16.48382;
else
    handles.background_axes = parsed_inputs.axes;
end
% Axes Positions:
guidata1.ax_pos_all_on =        [.04,  .036,  .875,  .93]; % All On
guidata1.ax_pos_all_off =       [0,     0,    1,    1]; % All Off
guidata1.ax_pos_no_title =      [.04,  .036,  .875,  1]; % No Title
guidata1.ax_pos_no_colorbar =   [.04,  .036,  .96,  .93]; % No Colorbar
guidata1.ax_pos_no_tick =       [0,     0,    .91, .97]; % No Tick
guidata1.ax_pos_title_only =    [0,     0,    1,    .97]; % Title only
guidata1.ax_pos_colorbar_only = [0,     0,    .91, 1];% Colorbar only
guidata1.ax_pos_tick_only =     [.04,  .036,  1,    1]; % Tick Only
guidata1.colorbar_pos = [.927,guidata1.ax_pos_all_on(2),0.0390,guidata1.ax_pos_all_on(4)];
% guidata1.x_ax_percent = .88; % for default view ('all on')
% guidata1.y_ax_percent = .93; % for default view ('all on')
guidata1.curr_axis_pos = guidata1.ax_pos_all_on;
guidata1.x_ax_percent = guidata1.curr_axis_pos(3);
guidata1.y_ax_percent = guidata1.curr_axis_pos(4);
guidata1.prev_state = 1; % various combo's of axes objects (see 'reposition_axes')
guidata1.curr_state = 1; % various combo's of axes objects
% Overlay Window Position:
guidata1.overlay_figure_pos = [guidata1.figure_pos(1),...
    guidata1.figure_pos(2)+5.5*guidata1.figure_pos(2),...
    guidata1.figure_pos(3),...
    guidata1.figure_pos(4)-.8*guidata1.figure_pos(4)];
% Colormap & Overlay Settings:
guidata1.update_colorbar = 0;
guidata1.overlay1_draw_color = 0;
guidata1.overlay2_draw_color = 0;
guidata1.prechange1_overlay1 = [];
guidata1.prechange2_overlay1 = [];
guidata1.redo_data1_overlay1 = [];
guidata1.redo_data2_overlay1 = [];
guidata1.prechange1_overlay2 = [];
guidata1.prechange2_overlay2 = [];
guidata1.redo_data1_overlay2 = [];
guidata1.redo_data2_overlay2 = [];
guidata1.overlay_colormaps = {'jet','hot','cool','spring','summer','winter'};
guidata1.overlay1_colormap_selection = 'jet'; % default: jet
guidata1.overlay2_colormap_selection = 'jet'; % default: jet
guidata1.overlay1_colorbar_on = false; % make sure to disable the other one while one is on
guidata1.overlay2_colorbar_on = false;
guidata1.m_string = '100';  % Number of colors in the current colormap
guidata1.m = str2double(guidata1.m_string);
% Colormaps & Colorbars:
guidata1.h_colorbar_main = 20.19347; % initialize handle to random
guidata1.h_colorbar_overlay1 = 21.38573; % initialize handle to random
guidata1.h_colorbar_overlay2 = 22.39553; % initialize handle to random
guidata1.main_colormap = colormap(eval([guidata1.main_colormap_selection,'(',guidata1.m_string,')']));
guidata1.drawing_colormap = [0,0,0; 0,1,0; 1,0,0]; % black (invisible),red,green,blue,cyan,magenta,white,black
guidata1.overlay1_colormap = eval([guidata1.overlay1_colormap_selection,'(',guidata1.m_string,')']);
guidata1.overlay2_colormap = eval([guidata1.overlay2_colormap_selection,'(',guidata1.m_string,')']); 
guidata1.overlay_alpha = 1; % overlay transparency (default: opaque)
guidata1.overlay1_alpha = parsed_inputs.overlay_alpha; % overlay transparency (default: 90%)
guidata1.overlay2_alpha = parsed_inputs.overlay_alpha; % overlay transparency (default: opaque)
guidata1.cmin_overlay1 = [];
guidata1.cmin_overlay2 = [];
guidata1.cmax_overlay1 = [];
guidata1.cmax_overlay2 = [];
guidata1.imgfull_idx = []; % actual colormap idx for image storage
guidata1.overlay1_idx = [];
guidata1.overlay2_idx = [];
guidata1.transparency_val1 = 1;
guidata1.transparency_val2 = 1;
guidata1.overlay1_colormap_val = 1;
guidata1.overlay2_colormap_val = 1;
guidata1.overlay_figure = 9.5869; % Initialize non-handle
% Drawing ROIs
guidata1.draw_path = pwd;
guidata1.curr_drawing = [];
guidata1.last_draw = []; % #vox x 3 matrix of indices [ix,iy,iz] of last drawing
guidata1.previously_editing = 'background';
guidata1.drawing_draw_color = 0;
guidata1.prechange1_draw = [];
guidata1.prechange2_draw = [];
guidata1.redo_data1_draw = [];
guidata1.redo_data2_draw = [];
guidata1.draw_filename = [];
guidata1.shape = 0;
guidata1.shape_ind = zeros(2,2);
guidata1.sphere_radius = 1;
guidata1.shape_coords = zeros(64,2);
guidata1.prev_shape_coords = [];
guidata1.shape_prev_colors = zeros(64,1);
% Other:
guidata1.customizable = {'guidata1.colorbar','guidata1.title_on',...
    'guidata1.axis_tick_on','guidata1.extensions','guidata1.last_nav_dir',...
    'guidata1.main_colormap_selection'};
guidata1.main_colorbar_on = true; % switches between main image colorbar and other axes colorbars
guidata1.first_dir = pwd; % remember initial Matlab path
guidata1.last_nav_dir = pwd;
guidata1.extensions = {'*.img';'*.nii';'*nii.gz'};
guidata1.slice_orientation = 3; % default = z-dim
guidata1.last_actions = [0,0]; % track whether recent changes were to background or overlay drawing
guidata1.h_title = 10.48487; % initialize non-handle
guidata1.num_voxels_crop = 0; % default until specified
guidata1.scroll_zoom_equiv = [1,.9,.8,.7,.6,.5,.4,.3,.2,.1];
% Set user interface callback functions:
set(handles.figure,'WindowKeyPressFcn',@keypress_callback);
% Don't turn on click function if embedded graphic in user spec axes:
if isempty(guidata1.parsed_inputs.axes)
    set(handles.figure,'WindowButtonDownFcn',@cursor_click_callback);
end
set(handles.figure,'WindowScrollWheelFcn',@scroll_zoom_callback);
if isempty(parsed_inputs.axes)
    set(handles.figure,'CloseRequestFcn',@closereq_callback)
else
    set(handles.figure,'CloseRequestFcn',@closereq_no_dlg)
end

% Save GUI Data:
guidata(handles.figure,guidata1);
get_settings; % extracts custom settings from .txt if available

guidata1.main_colormap = colormap(eval([guidata1.main_colormap_selection,...
    '(',guidata1.m_string,')']));

% Background: Get Input Filename & Load:
if isempty(guidata1.parsed_inputs.(guidata1.poss_input{1}))
    status = open_button_callback;
    if ~status
        return;
    end
elseif ischar(guidata1.parsed_inputs.(guidata1.poss_input{1}))
    [guidata1.path,guidata1.filename,ext] = fileparts(guidata1.parsed_inputs.(guidata1.poss_input{1}));
    guidata1.filename = [guidata1.filename,ext];
    guidata(handles.figure,guidata1);
    sort_exts(guidata1.path,guidata1.filename);
    load_img('char');
elseif isstruct(guidata1.parsed_inputs.(guidata1.poss_input{1}))
    load_img('struct');
end

%% Create Menu Items:
if isempty(parsed_inputs.axes)
    file_menu = uimenu(handles.figure,'Label','File');
        uimenu(file_menu,'Label','Open                                  ''o''','Callback',@open_button_callback);
        uimenu(file_menu,'Label','Save Image                        ''s''','Callback',@save_callback);
        uimenu(file_menu,'Label','Save Image As','Callback',@saveas_callback);
        uimenu(file_menu,'Label','Save Drawing','Callback',@save_drawing_callback); % add this callback
        uimenu(file_menu,'Label','Save Drawing As','Callback',@saveas_overlay_callback); % add this callback
        uimenu(file_menu,'Label','Exit                                      ''e''','Callback',@closereq_callback);
    edit_menu = uimenu(handles.figure,'Label','Edit');
        uimenu(edit_menu,'Label','Undo                         ''z''',...
            'Callback',@undo_callback);
        uimenu(edit_menu,'Label','Redo                          ''r''',...
            'Callback',@redo_callback);
        uimenu(edit_menu,'Label','Go to Slice...             ''g''',...
            'Callback',@go_to_slice_callback);
        uimenu(edit_menu,'Label','Revert to Defaults','Callback',@revert_defaults);
end
if isempty(parsed_inputs.axes)
    display_opts_menu = uimenu(handles.figure,'Label','Display');
else
    display_opts_menu = uimenu(handles.figure,'Label','NeuroImage Display');
end
    orient_menu = uimenu(display_opts_menu,'Label','Orientation');
        guidata1.h_x = uimenu(orient_menu,'Label','Coronal','Callback',@reorient_callback); 
        guidata1.h_y = uimenu(orient_menu,'Label','Sagittal','Callback',@reorient_callback);
        guidata1.h_z = uimenu(orient_menu,'Label','Axial','Callback',@reorient_callback);
    colormap_menu = uimenu(display_opts_menu,'Label','Colormap');
        uimenu(colormap_menu,'Label','gray','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','jet','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','hot','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','cool','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','hsv','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','bone','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','colorcube','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','copper','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','spring','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','summer','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','winter','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','pink','Callback',@colormap_callback);
    autoscale_menu = uimenu(display_opts_menu,'Label','Autoscale Slice Color');
        uimenu(autoscale_menu,'Label','On (default)','Callback',@autoscale_callback);
        autoscale_off_menu = uimenu(autoscale_menu,'Label','Off');
            uimenu(autoscale_off_menu,'Label','Use [min, max]','Callback',@autoscale_callback);
            uimenu(autoscale_off_menu,'Label','Specify...','Callback',@autoscale_callback);
    uimenu(display_opts_menu,'Label','Slice # (On/Off)','Callback',@title_toggle_callback);
if isempty(parsed_inputs.axes)
        uimenu(display_opts_menu,'Label','Colorbar (On/Off)','Callback',@colorbar_callback);
        uimenu(display_opts_menu,'Label','Axis Tick (On/Off)','Callback',@axis_tick_callback);
    draw_menu = uimenu(handles.figure,'Label','Draw');
        draw_color_menu = uimenu(draw_menu,'Label','Select Draw Color');
            uimenu(draw_color_menu,'Label','Erase','Callback',@draw_overlay_callback);
            uimenu(draw_color_menu,'Label','Signal ROI (Green)','Callback',@draw_overlay_callback);
            uimenu(draw_color_menu,'Label','Noise ROI (Red)','Callback',@draw_overlay_callback);
        draw_transparency_menu = uimenu(draw_menu,'Label','Select Draw Opacity');
            uimenu(draw_transparency_menu,'Label','Opaque (default)','Callback',@draw_transparency_callback);
            uimenu(draw_transparency_menu,'Label','80%','Callback',@draw_transparency_callback);
            uimenu(draw_transparency_menu,'Label','60%','Callback',@draw_transparency_callback);
            uimenu(draw_transparency_menu,'Label','40%','Callback',@draw_transparency_callback);
            uimenu(draw_transparency_menu,'Label','20%','Callback',@draw_transparency_callback);
            uimenu(draw_transparency_menu,'Label','Invisible','Callback',@draw_transparency_callback);
        draw_shapes_menu = uimenu(draw_menu,'Label','Shapes');
            uimenu(draw_shapes_menu,'Label','No Shape (Manual Trace)','Callback',{@draw_shapes_callback,1});
            uimenu(draw_shapes_menu,'Label','Circle','Callback',{@draw_shapes_callback,2});
            uimenu(draw_shapes_menu,'Label','Rectangle','Callback',{@draw_shapes_callback,3});
            uimenu(draw_shapes_menu,'Label','Sphere','Callback',{@draw_shapes_callback,4});
        uimenu(draw_menu,'Label','Propagate Through Slices','Callback',@propagate_draw_callback);    
        uimenu(draw_menu,'Label','Clear Drawing','Callback',@clear_drawing_callback);
    overlay_menu = uimenu(handles.figure,'Label','Overlay');
        uimenu(overlay_menu,'Label','Add/Edit Overlays','Callback',@add_overlay_callback);
    confirm_menu = uimenu(handles.figure,'Label','Confirm ROIs');
         uimenu(confirm_menu,'Label','Confirm','Callback',{@closereq_callback,1});
end

% Add customization settings if requested:
for ixx = 4:6
    input1 = guidata1.parsed_inputs.(guidata1.poss_input{ixx});
    if ~isempty(input1)
        if islogical(input1) || isnumeric(input1)
            eval(['guidata1.',poss_input{ixx},'=',num2str(input1),';']);
        end
    end
end

% Save GUI Data:
guidata1.refresh_img = true;
guidata(handles.figure,guidata1);

% Resize Figure Window Based on Aspect Ratio of Image:
update_image([],1);
if ~isempty(guidata1.parsed_inputs.colormap)
    guidata1.main_colormap_selection = guidata1.parsed_inputs.colormap;
    guidata(handles.figure,guidata1);
    update_colormap;
end
set(handles.figure,'Visible','on'); % Make Visible

% Add Overlays If Requested:
if ~isempty(guidata1.parsed_inputs.(guidata1.poss_input{2}))
    switch class(guidata1.parsed_inputs.(guidata1.poss_input{2}))
        case 'char'
            open_overlay_callback([], [], 1, 'char');
        case 'struct'
            open_overlay_callback([], [], 1, 'struct');
        otherwise
            if isnumeric(guidata1.parsed_inputs.(guidata1.poss_input{2}))
                open_overlay_callback([], [], 1, 'matrix');
            end
    end
end
if ~isempty(guidata1.parsed_inputs.(guidata1.poss_input{3}))
    switch class(guidata1.parsed_inputs.(guidata1.poss_input{3}))
        case 'char'
            open_overlay_callback([], [], 2, 'char');
        case 'struct'
            open_overlay_callback([], [], 2, 'struct');
        otherwise
            if isnumeric(guidata1.parsed_inputs.(guidata1.poss_input{3}))
                open_overlay_callback([], [], 2, 'matrix');
            end
    end
end

% Adjust CAXES if specified:
if ~isempty(parsed_inputs.background_caxis)
    if parsed_inputs.background_caxis(2)~=parsed_inputs.background_caxis(1) && parsed_inputs.background_caxis(2)>parsed_inputs.background_caxis(1)
        guidata1.cmin = parsed_inputs.background_caxis(1);
        guidata1.cmax = parsed_inputs.background_caxis(2);
        guidata1.autoscale = 0;
    else
        warning('Invalid background color axis specification. Ignoring inputs.')
        return;
    end
    guidata(handles.figure,guidata1);
    calc_full_ind
    guidata1.refresh_img = 0;
    guidata(handles.figure,guidata1);
    update_image
end

if ~isempty(parsed_inputs.overlay_caxis)
    if parsed_inputs.overlay_caxis(2)~=parsed_inputs.overlay_caxis(1) && parsed_inputs.overlay_caxis(2)>parsed_inputs.overlay_caxis(1)
        guidata1.cmin_overlay1 = parsed_inputs.overlay_caxis(1);
        guidata1.cmax_overlay1 =  parsed_inputs.overlay_caxis(2);
    else
        warning('Invalid background color axis specification. Ignoring inputs.')
        return;
    end
    guidata(handles.figure,guidata1);
    calc_full_ind
    guidata1.refresh_img = 0;
    guidata(handles.figure,guidata1);
    update_image
end

% If overlay is on, make it the active colorbar:
if guidata1.overlay1_on
    overlay1_colorbar_callback([],1,[]);
elseif guidata1.overlay2_on
    overlay2_colorbar_callback([],1,[]);
end

% Print if requested:
if ischar(guidata1.parsed_inputs.print)
    [requested_dir,~,~] = fileparts(guidata1.parsed_inputs.print);
    if ~isdir(requested_dir); return; end %#ok
    print_button_callback([],[],guidata1.parsed_inputs.print)
end

%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions:
%%%%%%%%%%%%%%%%%%%%%

% Load Image Function:
function load_img(img_type)
    guidata1 = guidata(handles.figure);
    switch img_type
        case 'struct'
            try
                guidata1.img = guidata1.parsed_inputs.(guidata1.poss_input{1});
                try
                    guidata1.path = guidata1.img.fileprefix;
                    [~,guidata1.filename,~] = fileparts(guidata1.path);
                catch
                    guidata1.filename = 'Image Structure';
                end
            catch
                disp('Error loading image structure. Check fields and retry or load through GUI.')
                return;
            end
        case 'char'
            guidata1.fname = fullfile(guidata1.path,guidata1.filename); 
            if guidata1.untouch_nii
                guidata1.img = load_untouch_nii(guidata1.fname);
            else
                try
                    guidata1.img = load_nii(guidata1.fname);
                catch
                    guidata1.img = load_untouch_nii(guidata1.fname);
                    guidata1.untouch_nii = true;
                    warning('Non-orthogonal shearing detected in affine matrix. Image loaded without applying affine.')
                end
            end
    end
    % Parse image input
    if isempty(parsed_inputs.title)
        guidata1.window_name = ['neuroimage_editor:    ',guidata1.filename];
    else
        guidata1.window_name = ['neuroimage_editor:                                  ',parsed_inputs.title];
    end 
    set(handles.figure,'Name',guidata1.window_name);  
    guidata1.dim = guidata1.img.hdr.dime.dim;
    guidata1.pixdim = guidata1.img.hdr.dime.pixdim;
    guidata1.xwidth = guidata1.dim(2)*guidata1.pixdim(2);
    guidata1.yheight = guidata1.dim(3)*guidata1.pixdim(3);             
    guidata1.imgfull = guidata1.img.img; 
    guidata1.imgfull = single(guidata1.imgfull); % conversion may enhance refresh speed
    guidata1.imgfull = permute(guidata1.imgfull,[2,1,3]); % Match SPM View
    [guidata1.xdim,guidata1.ydim,guidata1.z1] = size(guidata1.imgfull);
    guidata1.middle_slice = round(guidata1.z1/2);
    guidata1.curr_slice = guidata1.middle_slice;
    guidata1.drawing_active = false;
    guidata1.xind = guidata1.xdim:-1:1;
    guidata1.yind = 1:guidata1.ydim;
    guidata1.slice_orientation = 3; % default: Z-dim
    guidata1.cmin = min(guidata1.imgfull(:));
    guidata1.cmax = max(guidata1.imgfull(:));
    guidata1.xmax = max(guidata1.xind);
    guidata1.xmin = min(guidata1.xind);
    guidata1.ymax = max(guidata1.yind);
    guidata1.ymin = min(guidata1.yind);    
    guidata1.ax_xlim = [1,numel(guidata1.yind)];
    guidata1.ax_ylim = [1,numel(guidata1.xind)];
    if guidata1.cmin==guidata1.cmax
        guidata1.cmax = guidata1.cmin + 1;
    end
    guidata1.x_slice = []; guidata1.y_slice = []; guidata1.z_slice = []; 
    guidata1.curr_slice = guidata1.middle_slice;
    guidata1.in_motion = false;
    guidata1.background_draw_color = 0; % default
    guidata1.curr_drawing = [];
    guidata1.prechange1_main = [];
    guidata1.prechange2_main = [];
    guidata1.redo_data1_main = [];
    guidata1.redo_data2_main = [];
    % Clear Drawing and Overlays:
    guidata1.overlay1_on = false;
    guidata1.overlay2_on = false;
    if isgraphics(handles.overlay1_axes,'axes'); cla(handles.overlay1_axes); end
    if isgraphics(handles.overlay2_axes,'axes'); cla(handles.overlay2_axes); end
    if ishandle(handles.overlay1_img); cla(handles.overlay1_img); end
    if ishandle(handles.overlay2_img); cla(handles.overlay2_img); end
    % With these changes it now won't always load; if I Close overlay, it
    % gives colormap control back to main axes
    guidata1.drawing_active = false;
    guidata1.currently_editing = 'background';
    guidata1.shape = 0;
    guidata1.undo_num = 0;
    guidata1.scroll_count = 0;
    guidata1.drawing_idx = [];
    guidata1.drawing_overlay_alpha_data = []; 
    guidata1.overlay1_idx = [];
    guidata1.overlay2_idx = [];
    guidata1.draw_fname = [];
    guidata1.update_colorbar = true;
    guidata1.autoscale = true;
    if ishandle(handles.drawing_img)
        delete(handles.drawing_img);
    end
    if ishandle(handles.overlay1_img)
        delete(handles.overlay1_img);
    end
    if ishandle(handles.overlay2_img)
        delete(handles.overlay2_img);
    end
    if isgraphics(handles.background_img)
        delete(handles.background_img)
    end
    % Determine possible combinations of x,y indices:
    guidata1.poss_ind = zeros(guidata1.xdim*guidata1.ydim,2); count = 0;
    for ix = 1:guidata1.xdim
        for jx = 1:guidata1.ydim
            count = count + 1;
            guidata1.poss_ind(count,:) = [ix,jx];
        end
    end
    guidata(handles.figure,guidata1);
end

function sort_exts(path1,name1)
    guidata1 = guidata(handles.figure);
    guidata1.last_nav_dir = path1;
    [~,~,ext] = fileparts(name1);
    for ix = 1:3
        if ~isempty(strfind(guidata1.extensions{ix},ext))
            type = ix;
            others = setdiff(1:3,ix);
            break;
        end
    end
    guidata1.extensions = guidata1.extensions([type,others]);
    guidata(handles.figure,guidata1);
end

function resize_figure
    guidata1 = guidata(handles.figure);
    guidata1.aspect_ratio = (guidata1.xwidth*(1/guidata1.x_ax_percent))/(guidata1.yheight*(1/guidata1.y_ax_percent));
    guidata1.fig_width = guidata1.aspect_ratio*guidata1.fig_height;
    guidata1.figure_pos(1) = max(guidata1.screen_res(1)+8,...
        guidata1.screen_mid_x-(.5*guidata1.fig_width));
    guidata1.figure_pos(3) = min(guidata1.screen_res(3)-guidata1.figure_pos(1)-7,...
        guidata1.fig_width);
    handles.figure.Position = guidata1.figure_pos;
    guidata1.overlay_figure_pos = [guidata1.figure_pos(1),...
        guidata1.figure_pos(2)+5.5*guidata1.figure_pos(2),...
        guidata1.figure_pos(3),guidata1.figure_pos(4)-.8*guidata1.figure_pos(4)];
    guidata(handles.figure,guidata1);
end

function get_settings
    guidata1 = guidata(handles.figure);
    listing = dir(fullfile(guidata1.script_path,'neuroimage_editor_settings.txt'));
    if isempty(listing)
        disp('No customized settings file found: missing file ''neuroimage_editor_settings.txt''.')
        return;
    end
    gui_settings = importdata(fullfile(guidata1.script_path,listing(1).name));
    for ix1 = 1:length(gui_settings)
        ind = strfind(gui_settings{ix1},'=');
        eval([guidata1.customizable{ix1}, '=[', gui_settings{ix1}(ind+1:end),'];']);
    end
    if ~isdir(guidata1.last_nav_dir); guidata1.last_nav_dir = pwd; end %#ok 
    guidata(handles.figure,guidata1);
end

function write_settings
    guidata1 = guidata(handles.figure);
    % Change Settings Data:
    gui_settings = cell(6,1);
    for ix1 = 1:3
        setting1 = eval(guidata1.customizable{ix1});
        if length(setting1)==1
            gui_settings{ix1} = [guidata1.customizable{ix1},'=',num2str(setting1),';'];
        else
            gui_settings{ix1} = [guidata1.customizable{ix1},'=1;'];
        end     
    end
    exts = eval(guidata1.customizable{4});
    gui_settings{4} = [guidata1.customizable{4},'={''',exts{1},''';''',exts{2},''';''',exts{3},'''};'];
    for ix1 = 5:6
        gui_settings{ix1} = [guidata1.customizable{ix1},'=''',eval(guidata1.customizable{ix1}),''';'];
    end
    gui_settings=char(gui_settings);
    % Write Text File:
    fileID = fopen(fullfile(guidata1.script_path,'neuroimage_editor_settings.txt'),'w');
    for ix1 = 1:6
        fprintf(fileID,'%s\r\n',gui_settings(ix1,:));
    end
    fclose(fileID);
    guidata(handles.figure,guidata1);
end

function revert_defaults(~,~,~)
    guidata1 = guidata(handles.figure);
    guidata1.colorbar=1;                                                                                             
    guidata1.title_on=1;                                                                                             
    guidata1.axis_tick_on=1;                                                                                         
    guidata1.extensions={'*.img';'*.nii';'*nii.gz'};                                                                 
    guidata1.last_nav_dir=guidata1.first_dir;
    guidata1.main_colormap_selection='gray';    
    guidata1.refresh_img = true;
    guidata(handles.figure,guidata1);
    update_colormap;
    update_image;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Button Callback Functions: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% File Menu:

% Exit Request Callback
function closereq_callback(~,~,no_confirm) %#ok
%     if nargin<3 || ~no_confirm
%         selection = questdlg('Are you sure you want to exit?','Exit GUI',...
%           'Yes','No','Yes');
%     else
%         selection = 'Yes';
%     end
%     switch selection,
%         case 'Yes'
            guidata1 = guidata(handles.figure);
            write_settings % Save Custom Settings
            % Write 2nd & 3rd Outputs to Base, if Requested:
            if guidata1.nargout>1 && ~isempty(guidata1.last_draw)
                % Switch x-y:
                hold1 = guidata1.last_draw(:,1);
                guidata1.last_draw(:,1) = guidata1.last_draw(:,2);
                guidata1.last_draw(:,2) = hold1;
                assignin('base','last_draw',guidata1.last_draw)
                assignin('caller','last_draw',guidata1.last_draw)
            end
            if guidata1.nargout>2
                assignin('base','drawing_idx',guidata1.drawing_idx)
                assignin('caller','drawing_idx',guidata1.drawing_idx)
            end
            % Delete GUI & Figures:
            if ishandle(guidata1.overlay_figure); delete(guidata1.overlay_figure); end
            delete(handles.figure)
%         case 'No'
%           return
%     end   
end

% Exit Request No Dialogue:
function closereq_no_dlg(~, ~, ~)
    guidata1 = guidata(handles.figure);
%     write_settings % Save Custom Settings
    % Write 2nd & 3rd Outputs to Base, if Requested:
%     if guidata1.nargout>1
%         assignin('base','last_draw',guidata1.last_draw)
%         assignin('caller','last_draw',guidata1.last_draw)
%     end
%     if guidata1.nargout>2
%         assignin('base','drawing_idx',guidata1.drawing_idx)
%         assignin('caller','drawing_idx',guidata1.drawing_idx)
%     end
    % Delete GUI & Figures:
    if ishandle(guidata1.overlay_figure); delete(guidata1.overlay_figure); end
    delete(handles.figure)
end

% "Open" Callback:
function status = open_button_callback(~,~,~)
    guidata1 = guidata(handles.figure);
    % Get input filename:
    cd(guidata1.last_nav_dir);
    [guidata1.filename, guidata1.path] = uigetfile(guidata1.extensions,...
        'Select Image File:','MultiSelect','off');
    cd(guidata1.first_dir);
    guidata(handles.figure,guidata1);
    if guidata1.path ~= 0 % if not 'cancel'
        sort_exts(guidata1.path,guidata1.filename);
        load_img('char');
        % Save and Update:
        guidata1.refresh_img = true;
        guidata(handles.figure,guidata1);
        update_image([],1);
        update_colormap;
        status = 1;
    else
        status = 0;
        return;
    end
end

% "Save" Button Callback:
function save_callback(~,~,~)
    guidata1 = guidata(handles.figure);
    if guidata1.slice_orientation~=3
        reorient_callback(guidata1.h_z) % revert to original orientation
    end
    guidata1.img.img = permute(guidata1.imgfull,[2,1,3]);
    if guidata1.untouch_nii
        save_untouch_nii(guidata1.img,guidata1.fname);
    else 
        save_nii(guidata1.img,guidata1.fname);
    end
    disp(['Image successfully saved as ',guidata1.fname])
    guidata(handles.figure,guidata1);
end

% "Save As" Button Callback:
function saveas_callback(~,~,~)
    guidata1 = guidata(handles.figure);
    if guidata1.slice_orientation~=3
        reorient_callback(guidata1.h_z) % revert to original orientation
    end
    guidata1.img.img = permute(guidata1.imgfull,[2,1,3]);
    cd(guidata1.last_nav_dir);
    [guidata1.filename,PathName] = uiputfile(guidata1.extensions,'Specify filename:');
    cd(guidata1.first_dir);
    if isdir(PathName); guidata1.last_nav_dir = PathName; end; %#ok
    if ischar(guidata1.filename)
        guidata1.fname = fullfile(PathName,guidata1.filename);
        if guidata1.untouch_nii
            save_untouch_nii(guidata1.img,guidata1.fname);
        else
            save_nii(guidata1.img,guidata1.fname);
        end
        % Change Figure Window Title to Reflect New Filename:
        guidata1.window_name = ['neuroimage_editor:    ',guidata1.filename];
        set(handles.figure,'Name',guidata1.window_name)
    end
    disp(['Image successfully saved as ',guidata1.fname])
    guidata(handles.figure,guidata1);
end

% "Save Drawing":
function save_drawing_callback(~,~,~)
    guidata1 = guidata(handles.figure);
    if guidata1.slice_orientation~=3
        reorient_callback(guidata1.h_z) % revert to original orientation
    end
    guidata1.img.img = permute(guidata1.drawing_idx,[2,1,3]);
    nonzero = guidata1.img.img(:)~=0;
    guidata1.img.img(nonzero) = guidata1.img.img(nonzero)-1;
    if isempty(guidata1.draw_fname)
        cd(guidata1.last_nav_dir);
        [guidata1.draw_filename,guidata1.draw_path] = uiputfile(guidata1.extensions,'Specify filename:');
        cd(guidata1.first_dir);
        if ischar(guidata1.draw_filename)
            guidata1.draw_fname = fullfile(guidata1.draw_path,guidata1.draw_filename);
            if guidata1.untouch_nii
                save_untouch_nii(guidata1.img,guidata1.draw_fname);
            else
                save_nii(guidata1.img,guidata1.draw_fname);
            end
        end
    elseif ischar(guidata1.draw_fname)
        if guidata1.untouch_nii
            save_untouch_nii(guidata1.img,guidata1.draw_fname);
        else
            save_nii(guidata1.img,guidata1.draw_fname);
        end
    end
    disp(['Overlay drawing successfully saved as ',guidata1.draw_fname])
    guidata(handles.figure,guidata1);
end

% "Save Overlay Drawing As"
function saveas_overlay_callback(~,~,~)
    guidata1 = guidata(handles.figure);
    if guidata1.slice_orientation~=3
        reorient_callback(guidata1.h_z) % revert to original orientation
    end
    guidata1.img.img = permute(guidata1.drawing_idx,[2,1,3]);
    nonzero = guidata1.img.img(:)~=0;
    guidata1.img.img(nonzero) = guidata1.img.img(nonzero)-1;
    cd(guidata1.last_nav_dir);
    [guidata1.draw_filename,guidata1.draw_path] = uiputfile(guidata1.extensions,'Specify filename:');
    cd(guidata1.first_dir);
    if ischar(guidata1.draw_filename)
        guidata1.draw_fname = fullfile(guidata1.draw_path,guidata1.draw_filename);
        if guidata1.untouch_nii
            save_untouch_nii(guidata1.img,guidata1.draw_fname);
        else
            save_nii(guidata1.img,guidata1.draw_fname);
        end
    end
    disp(['Overlay drawing successfully saved as ',guidata1.draw_fname])
    guidata(handles.figure,guidata1);
end

% Print Button Callback:
function print_button_callback(~, ~, figure_title)
    guidata1 = guidata(handles.figure);
    if nargin==3 && ischar(figure_title)
        [~,~,file_ext] = fileparts(figure_title);
        % Identify File Extension:
        exts = {'-dpng','-dtiff','-dbmp'};
        for i = 1:3
            if ~isempty(strfind(exts{i},file_ext(2:end)))
                ext_ind = i; break;
            end
        end
        % Print Figure:
        set(handles.figure,'InvertHardcopy','off','PaperPositionMode','auto')
        print(handles.figure,exts{ext_ind},'-r300','-loose',figure_title)
        set(handles.figure,'CloseRequestFcn',@closereq_no_dlg)
        disp(['Printed figure: ',figure_title])
        return;
    end
    ext_opts = {'.png',' ';'.tiff',' ';'.bmp',' '};
    [title1,path1] = uiputfile(ext_opts,'Specify filename:');
    figure_title = fullfile(path1,title1);
    [~,~,file_ext] = fileparts(figure_title);
    % Identify File Extension:
    exts = {'-dpng','-dtiff','-dbmp'};
    for i = 1:3
        if ~isempty(strfind(exts{i},file_ext(2:end)))
            ext_ind = i; break;
        end
    end
%     figure_pos = [.2*guidata1.screen_res(3), .065*guidata1.screen_res(4), .605*guidata1.screen_res(3), .895*guidata1.screen_res(4)];
    print_figure = figure('Position',guidata1.figure_pos,'MenuBar','none','Name','neuroimage_editor','NumberTitle','off','Color',[1,1,1]);
    % Regenerate Figure:
    if guidata1.autoscale
        background_slice = guidata1.imgfull(guidata1.xind,guidata1.yind,guidata1.curr_slice);
        cmin1 = min(background_slice(:));  % Minimum color value
        cmax1 = max(background_slice(:));  % Maximum color value
        idx1 = min(guidata1.m,round((guidata1.m-1)*(background_slice-cmin1)/(cmax1-cmin1))+1);
        image(idx1);
    else
        image(guidata1.imgfull_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice));
    end
    hold on;
    h_axes = gca;
    colormap(print_figure,guidata1.main_colormap);
    set(h_axes,'Position',[0,0,1,1]);
    % Regenerate Drawing and Overlays:
    if ishandle(handles.drawing_img)
        h_drawing = axes('Visible','off','Position',[0,0,1,1]);
        image('Parent',h_drawing,'CData',guidata1.drawing_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
            'AlphaData',guidata1.drawing_overlay_alpha_data(guidata1.xind,guidata1.yind,guidata1.curr_slice));
    end
    if ishandle(handles.overlay1_img)
        overlay1_axes = axes('Visible','off','Position',[0,0,1,1]);
        image('Parent',overlay1_axes,'CData',guidata1.overlay1_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
            'AlphaData',guidata1.overlay1_alpha_data(guidata1.xind,guidata1.yind,guidata1.curr_slice));
        colormap(overlay1_axes,guidata1.overlay1_colormap)
    end
    if ishandle(handles.overlay2_img)
        overlay2_axes = axes('Visible','off','Position',[0,0,1,1]);
        image('Parent',overlay2_axes,'CData',guidata1.overlay2_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
            'AlphaData',guidata1.overlay1_alpha_data(guidata1.xind,guidata1.yind,guidata1.curr_slice));
        colormap(overlay2_axes,guidata1.overlay2_colormap)
    end
    % Print Figure:
    set(print_figure,'InvertHardcopy','off','PaperPositionMode','auto')
    print(print_figure,exts{ext_ind},'-r300','-loose',figure_title)
    close(print_figure);
    disp(['Printed figure: ',figure_title])
end

function title_toggle_callback(~,~,~)
    guidata1 = guidata(handles.figure);
    guidata1.title_on = ~guidata1.title_on;
    guidata1.refresh_img = true;
    guidata(handles.figure,guidata1);
    update_image
    reposition_axes
end

function axis_tick_callback(~,~,~)
    guidata1 = guidata(handles.figure);
    guidata1.axis_tick_on = ~guidata1.axis_tick_on;
    guidata1.refresh_img = true;
    guidata(handles.figure,guidata1);
    reposition_axes
end

%% Image Orientation

% Reorient Image (axial, coronal, sagittal):
function reorient_callback(hObject, ~, ~)
    guidata1 = guidata(handles.figure);
    reorient_type = get(hObject,'Label');
%     if nargin==1
%         reorient_type = 'Axial'; % If called from Load Image
%     end
    switch reorient_type
        case 'Coronal'
            if guidata1.slice_orientation~=1
                if guidata1.slice_orientation==2
                    guidata1.dimperm = [1,3,2];
                    guidata1.y_slice = guidata1.curr_slice;
                elseif guidata1.slice_orientation==3
                    guidata1.dimperm = [3,2,1];
                    guidata1.z_slice = guidata1.curr_slice;
                end       
                guidata1.slice_orientation = 1;
                if isempty(guidata1.x_slice)
                    guidata1.curr_slice = [];
                else
                    guidata1.curr_slice = guidata1.x_slice;
                end
            else
                return;
            end
            guidata1.xwidth = guidata1.dim(4)*guidata1.pixdim(4);
            guidata1.yheight = guidata1.dim(3)*guidata1.pixdim(3);
        case 'Sagittal'
            if guidata1.slice_orientation~=2
                if guidata1.slice_orientation==1
                    guidata1.dimperm = [1,3,2];
                    guidata1.x_slice = guidata1.curr_slice;
                elseif guidata1.slice_orientation==3
                    guidata1.dimperm = [3,1,2];
                    guidata1.z_slice = guidata1.curr_slice;
                end
                guidata1.slice_orientation = 2;
                if isempty(guidata1.y_slice)
                    guidata1.curr_slice = [];
                else
                    guidata1.curr_slice = guidata1.y_slice;
                end
            else
                return;
            end
            guidata1.xwidth = guidata1.dim(4)*guidata1.pixdim(4);
            guidata1.yheight = guidata1.dim(2)*guidata1.pixdim(2);
        case 'Axial'
            if guidata1.slice_orientation~=3
                if guidata1.slice_orientation==1
                    guidata1.dimperm = [3,2,1]; 
                    guidata1.x_slice = guidata1.curr_slice;
                elseif guidata1.slice_orientation==2
                    guidata1.dimperm = [2,3,1];
                    guidata1.y_slice = guidata1.curr_slice;
                end
                guidata1.slice_orientation = 3;  
                if isempty(guidata1.z_slice)
                    guidata1.curr_slice = [];
                else
                    guidata1.curr_slice = guidata1.z_slice;
                end
            else
                return;
            end
            guidata1.xwidth = guidata1.dim(2)*guidata1.pixdim(2);
            guidata1.yheight = guidata1.dim(3)*guidata1.pixdim(3);
    end
    guidata1.imgfull = permute(guidata1.imgfull,guidata1.dimperm);
    guidata1.imgfull_idx = permute(guidata1.imgfull_idx,guidata1.dimperm);
    if guidata1.drawing_active          
        guidata1.drawing_idx = permute(guidata1.drawing_idx,guidata1.dimperm);
        guidata1.drawing_overlay_alpha_data = permute(guidata1.drawing_overlay_alpha_data,guidata1.dimperm);
    end
    if guidata1.overlay1_on
        guidata1.overlay1_imgfull = permute(guidata1.overlay1_imgfull,guidata1.dimperm);
        guidata1.overlay1_idx = permute(guidata1.overlay1_idx,guidata1.dimperm);
        guidata1.overlay1_alpha_data = permute(guidata1.overlay1_alpha_data,guidata1.dimperm);
    end
    if guidata1.overlay2_on
        guidata1.overlay2_imgfull = permute(guidata1.overlay2_imgfull,guidata1.dimperm);
        guidata1.overlay2_idx = permute(guidata1.overlay2_idx,guidata1.dimperm);
        guidata1.overlay2_alpha_data = permute(guidata1.overlay2_alpha_data,guidata1.dimperm);
    end
    [guidata1.xdim,guidata1.ydim,guidata1.z1] = size(guidata1.imgfull);
%     guidata1.xind = 1:guidata1.xdim;
    guidata1.xind = guidata1.xdim:-1:1;
    if guidata1.slice_orientation==2
        guidata1.yind = guidata1.ydim:-1:1;
    else
        guidata1.yind = 1:guidata1.ydim;
    end
    guidata1.middle_slice = round(guidata1.z1/2);
    if isempty(guidata1.curr_slice)
        guidata1.curr_slice = guidata1.middle_slice;
    end
    % Permute Undo/Redo Data:
    guidata1.prechange1_main = permute(guidata1.prechange1_main,guidata1.dimperm);
    guidata1.prechange2_main = permute(guidata1.prechange2_main,guidata1.dimperm);
    guidata1.redo_data1_main = permute(guidata1.redo_data1_main,guidata1.dimperm);
    guidata1.redo_data2_main = permute(guidata1.redo_data2_main,guidata1.dimperm);
    guidata1.prechange1_draw = permute(guidata1.prechange1_draw,guidata1.dimperm);
    guidata1.prechange2_draw = permute(guidata1.prechange2_draw,guidata1.dimperm);
    guidata1.redo_data1_draw = permute(guidata1.redo_data1_draw,guidata1.dimperm);
    guidata1.redo_data2_draw = permute(guidata1.redo_data2_draw,guidata1.dimperm);
    guidata1.prechange1_overlay1 = permute(guidata1.prechange1_overlay1,guidata1.dimperm);
    guidata1.prechange2_overlay1 = permute(guidata1.prechange2_overlay1,guidata1.dimperm);
    guidata1.redo_data1_overlay1 = permute(guidata1.redo_data1_overlay1,guidata1.dimperm);
    guidata1.redo_data2_overlay1 = permute(guidata1.redo_data2_overlay1,guidata1.dimperm);
    guidata1.prechange1_overlay2 = permute(guidata1.prechange1_overlay2,guidata1.dimperm);
    guidata1.prechange2_overlay2 = permute(guidata1.prechange2_overlay2,guidata1.dimperm);
    guidata1.redo_data1_overlay2 = permute(guidata1.redo_data1_overlay2,guidata1.dimperm);
    guidata1.redo_data2_overlay2 = permute(guidata1.redo_data2_overlay2,guidata1.dimperm);    
    % Save Guidata and Update Image:
    guidata1.refresh_img = false;
    guidata(handles.figure,guidata1);
    update_image
    reposition_axes
    resize_figure
end

%% Edit Underlying Data (Background Figure) Colorspec, Colormap, Colorbar, Auto-Scale:

% Colormap Callback:
function colormap_callback(hObject, ~, ~)
    guidata1 = guidata(handles.figure);
    guidata1.main_colormap_selection = get(hObject,'Label');
    guidata(handles.figure,guidata1);
    update_colormap
    update_image
end

% Colorbar Callback:
function colorbar_callback(~, ~, ~)
    guidata1 = guidata(handles.figure);
    guidata1.colorbar = ~guidata1.colorbar;
    guidata1.refresh_img = true;
    guidata1.update_colorbar = true;
    guidata(handles.figure,guidata1);
    update_image
    reposition_axes
end

% Change Background Color Limits
function autoscale_callback(hObject, ~, ~) 
    guidata1 = guidata(handles.figure);
    opt = get(hObject,'Label');
    switch opt
        case 'On (default)'
            guidata1.autoscale = true;
            guidata1.cmin = [];
            guidata1.cmax = [];
            guidata1.imgfull_idx = [];
            guidata(handles.figure,guidata1);
            calc_full_ind
        case 'Use [min, max]'
            guidata1.autoscale = false;
            guidata1.cmin = min(guidata1.imgfull(:));
            guidata1.cmax = max(guidata1.imgfull(:));
            if guidata1.cmin==guidata1.cmax
                guidata1.cmax = guidata1.cmin + 1;
            end
            guidata(handles.figure,guidata1);
            calc_full_ind
        case 'Specify...'
            guidata1.autoscale = false;
            cmin1 = min(guidata1.imgfull(:));
            cmax1 = max(guidata1.imgfull(:));
            prompt = {['Specify range intensity values for color axis: ',...
                char(10),'Min:'],'Max:'}; %#ok
            dlg_title = 'Input CAxis Limits'; num_lines = [1,20;1,20]; 
            defaultans = {num2str(cmin1),num2str(cmax1)};
            answer1 = inputdlg(prompt,dlg_title,num_lines,defaultans);
            if ~isempty(answer1)
                cmin2 = str2double(answer1{1});
                cmax2 = str2double(answer1{2});
            else
                return
            end
            if ~isempty(cmin2) && ~isempty(cmax2) && cmin2~=cmax2
                guidata1.cmin = cmin2;
                guidata1.cmax = cmax2;
            else
                guidata1.cmin = cmin1;
                guidata1.cmax = cmax1;
                disp('Invalid color axis specification.')
            end
            guidata(handles.figure,guidata1);
            calc_full_ind
    end
    guidata1.refresh_img = 0;
    guidata(handles.figure,guidata1);
    update_image
end

%% Undo, Redo, Go to Slice Crop Border, Fill Image:

% Undo Callback:
function undo_callback(~, ~, ~) 
    guidata1 = guidata(handles.figure);
    draw_overlay_undo = false;
    if guidata1.undo_num==0
        switch guidata1.last_actions(2)
            case 1
                if ~isempty(guidata1.prechange1_main)
                    guidata1.redo_data1_main = guidata1.imgfull;
                    guidata1.imgfull = guidata1.prechange1_main;
                    guidata1.undo_num = -1;
                    draw_overlay_undo = false;
                end
            case 2
                if ~isempty(guidata1.prechange1_draw)
                    guidata1.redo_data1_draw = guidata1.drawing_idx;
                    guidata1.drawing_idx = guidata1.prechange1_draw;
                    guidata1.undo_num = -1;
                    draw_overlay_undo = true;
                end
            case 3
                if ~isempty(guidata1.prechange1_overlay1)
                    guidata1.redo_data1_overlay1 = guidata1.overlay1_imgfull;
                    guidata1.overlay1_imgfull = guidata1.prechange1_overlay1;
                    guidata1.undo_num = -1;
                    draw_overlay_undo = false;
                end
            case 4
                if ~isempty(guidata1.prechange1_overlay2)
                    guidata1.redo_data1_overlay2 = guidata1.overlay2_imgfull;
                    guidata1.overlay2_imgfull = guidata1.prechange1_overlay2;
                    guidata1.undo_num = -1;
                    draw_overlay_undo = false;
                end
        end
    elseif guidata1.undo_num==-1
        switch guidata1.last_actions(1)
            case 1
                if ~isempty(guidata1.prechange2_main)
                    guidata1.redo_data2_main = guidata1.imgfull;
                    guidata1.imgfull = guidata1.prechange2_main;
                    guidata1.undo_num = -2;
                    draw_overlay_undo = false;
                end
            case 2
                if ~isempty(guidata1.prechange2_draw)
                    guidata1.redo_data2_draw = guidata1.drawing_idx;
                    guidata1.drawing_idx = guidata1.prechange2_draw;
                    guidata1.undo_num = -2;
                    draw_overlay_undo = true;
                end
            case 3
                if ~isempty(guidata1.prechange2_overlay1)
                    guidata1.redo_data2_overlay1 = guidata1.overlay1_imgfull;
                    guidata1.overlay1_imgfull = guidata1.prechange2_overlay1;
                    guidata1.undo_num = -2;
                    draw_overlay_undo = false;
                end
            case 4
                if ~isempty(guidata1.prechange2_overlay2)
                    guidata1.redo_data2_overlay2 = guidata1.overlay2_imgfull;
                    guidata1.overlay2_imgfull = guidata1.prechange2_overlay2;
                    guidata1.undo_num = -2;
                    draw_overlay_undo = false;
                end
        end
    else
        return;
    end
    guidata(handles.figure,guidata1);
    if draw_overlay_undo
        calc_full_ind(1)
        update_image(1)
    else
        calc_full_ind
        update_image
    end
end

% Redo Callback:
function redo_callback(~,~,~)
    guidata1 = guidata(handles.figure);
    if guidata1.undo_num==-1
        switch guidata1.last_actions(2)
            case 1
                if ~isempty(guidata1.redo_data1_main)
                    guidata1.imgfull = guidata1.redo_data1_main;
                    guidata1.undo_num = 0;
                    draw_overlay_undo = false;
                end
            case 2
                if ~isempty(guidata1.redo_data1_draw)
                    guidata1.drawing_idx = guidata1.redo_data1_draw;
                    guidata1.undo_num = 0;
                    draw_overlay_undo = true;
                end
            case 3
                if isempty(guidata1.redo_data1_overlay1)
                    guidata1.overlay1_imgfull = guidata1.redo_data1_overlay1;
                    guidata1.undo_num = 0;
                    draw_overlay_undo = false;
                end               
            case 4
                if isempty(guidata1.redo_data1_overlay2)
                    guidata1.overlay2_imgfull = guidata1.redo_data1_overlay2;
                    guidata1.undo_num = 0;
                    draw_overlay_undo = false;
                end                      
        end
    elseif guidata1.undo_num==-2
        switch guidata1.last_actions(1)
            case 1
                if ~isempty(guidata1.redo_data2_main)
                    guidata1.imgfull = guidata1.redo_data2_main;
                    guidata1.undo_num = -1;
                    draw_overlay_undo = false;
                end
            case 2
                if ~isempty(guidata1.redo_data2_draw)
                    guidata1.drawing_idx = guidata1.redo_data2_draw;
                    guidata1.undo_num = -1;
                    draw_overlay_undo = true;
                end
            case 3
                if ~isempty(guidata1.redo_data2_overlay1)
                    guidata1.overlay1_imgfull = guidata1.redo_data2_overlay1;
                    guidata1.undo_num = -1;
                    draw_overlay_undo = false;
                end                
            case 4
                if ~isempty(guidata1.redo_data2_overlay2)
                    guidata1.overlay2_imgfull = guidata1.redo_data2_overlay2;
                    guidata1.undo_num = -1;
                    draw_overlay_undo = false;
                end                           
        end
    else
        return;
    end
    guidata(handles.figure,guidata1);
    if draw_overlay_undo
        calc_full_ind(1)
        update_image(1)
    else
        calc_full_ind
        update_image
    end
end

% Go to Slice Callback:
function go_to_slice_callback(~,~,~)
    guidata1 = guidata(handles.figure);
    prompt = {'Specify slice: '};
    dlg_title = 'Go to Slice'; num_lines = [1,15;]; defaultans = {num2str(guidata1.z1)};
    answer1 = inputdlg(prompt,dlg_title,num_lines,defaultans);
    if isempty(answer1); return; end
    slice_spec = str2double(answer1(1)); 
    if slice_spec <= guidata1.z1 && slice_spec > 0
        guidata1.curr_slice = slice_spec;
    else
        disp(['Slice specification out of range. Slice # must be between 1 and ',...
            num2str(guidata1.z1),'.']);
        return;
    end
    guidata(handles.figure,guidata1);
    update_image
end

% Crop Border Button Callback:
function crop_border_callback(~,~,~)
    prompt = {'# voxels from border: '};
    dlg_title = 'Crop'; num_lines = [1,25];
    defaultans = {'0'};
    answer1 = inputdlg(prompt,dlg_title,num_lines,defaultans,'on');
    if ~isempty(answer1)
        guidata1 = guidata(handles.figure);
        guidata1.num_voxels_crop = str2double(answer1);
        if guidata1.num_voxels_crop ~= 0
            guidata1.prechange2_main = guidata1.prechange1_main;
            guidata1.prechange1_main = guidata1.imgfull;
            guidata1.last_actions(1) = guidata1.last_actions(2);
            guidata1.last_actions(2) = 1; % background edit
            guidata1.undo_num = 0;
            slice_data = guidata1.imgfull(guidata1.xind,guidata1.yind,guidata1.curr_slice);
            slice_data(1:guidata1.num_voxels_crop,:) = guidata1.background_draw_color;
            slice_data(end-guidata1.num_voxels_crop:end,:) = guidata1.background_draw_color;
            slice_data(:,1:guidata1.num_voxels_crop) = guidata1.background_draw_color;
            slice_data(:,end-guidata1.num_voxels_crop:end) = guidata1.background_draw_color;
            guidata1.imgfull(guidata1.xind,guidata1.yind,guidata1.curr_slice) = slice_data;
            guidata1.currently_editing = 'background';
            guidata(handles.figure,guidata1);
            calc_full_ind
            update_image
        end
    end
end    

%% Create User Interface Callbacks (clicks, arrow keys, mouse scrolls)

% Arrow Key Callback Function:
function keypress_callback(varargin)
    keypress = varargin{2}.Key;
    switch keypress
        % Up Arrow Key
        case 'uparrow' 
            guidata1 = guidata(handles.figure);
            if guidata1.curr_slice < guidata1.z1
                guidata1.curr_slice = guidata1.curr_slice + 1;
                guidata(handles.figure,guidata1);
                update_image
            end
        % Down Arrow Key
        case 'downarrow'
            guidata1 = guidata(handles.figure);
            if guidata1.curr_slice > 1
                guidata1.curr_slice = guidata1.curr_slice - 1;
                guidata(handles.figure,guidata1);
                update_image
            end
        case 'o'
            open_button_callback
        case 's'
            save_callback
        case 'e'
            closereq_callback
        case 'z'
            undo_callback
        case 'r'
            redo_callback
        case 'g'
            go_to_slice_callback
        case 'c'
            crop_border_callback
    end     
end

% Cursor Click Callback:
function cursor_click_callback(~,~,~)
    if strcmpi(get(gco,'type'),'image') 
        % Retrieve Data
        guidata1 = guidata(handles.figure);
        % Get mouse location
        pt = get(handles.background_axes, 'CurrentPoint');
        pt = round(pt(1,1:2));
        if pt(1) > 0 && pt(2) > 0
            % Set Callbacks
            set(handles.figure, 'WindowButtonUpFcn', @unclick_callback);
            set(handles.figure, 'WindowButtonMotionFcn', @motion_callback);
            % Add color:
            guidata1.last_actions(1) = guidata1.last_actions(2);
            switch guidata1.currently_editing
                case 'background'
                    guidata1.last_actions(2) = 1; % background edit
                    guidata1.prechange2_main = guidata1.prechange1_main;
                    guidata1.prechange1_main = guidata1.imgfull; % update cache
                    guidata1.imgfull(guidata1.xind(pt(2)),guidata1.yind(pt(1)),guidata1.curr_slice) = guidata1.background_draw_color; % update data
                case 'drawing'
                    guidata1.last_actions(2) = 2; % drawing edit
                    guidata1.prechange2_draw = guidata1.prechange1_draw;
                    guidata1.prechange1_draw = guidata1.drawing_idx; % update cache
                    if guidata1.shape == 0
                        guidata1.drawing_idx(guidata1.xind(pt(2)),guidata1.yind(pt(1)),guidata1.curr_slice) = guidata1.drawing_draw_color; % update data
                    else
                        guidata1.shape_ind(1,:) = [guidata1.xind(pt(2)),guidata1.yind(pt(1))];
                    end
                    if guidata1.shape==3 % sphere
                        pixdim = guidata1.pixdim(2:4);
                        dim = guidata1.dim(2:4);
                        s_center = [guidata1.xind(pt(2)),guidata1.yind(pt(1)),guidata1.curr_slice];
                        s_center = s_center.*pixdim; % convert vox to physical units
                        Min = s_center - [guidata1.sphere_radius,guidata1.sphere_radius,guidata1.sphere_radius];
                        Min = fix(max(pixdim,Min));
                        Max = s_center + [guidata1.sphere_radius,guidata1.sphere_radius,guidata1.sphere_radius];
                        Max = ceil(min(dim.*pixdim,Max));
                        guidata1.curr_drawing = zeros(size(guidata1.drawing_idx)); % used for propagating last drawing
                        for ix = Min(1):pixdim(1):Max(1) %Min(1):Max(1)
                            for jx = Min(2):pixdim(2):Max(2)  %Min(2):Max(2)
                                for kx = Min(3):pixdim(3):Max(3)  %Min(3):Max(3)
                                    if sqrt(sum(([ix,jx,kx] - s_center).^2 )) <= guidata1.sphere_radius
                                        x_coord = ceil(ix/pixdim(1));
                                        y_coord = ceil(jx/pixdim(2));
                                        z_coord = ceil(kx/pixdim(3));
                                        guidata1.shape_prev_colors(x_coord,y_coord,z_coord) = guidata1.drawing_idx(x_coord,y_coord,z_coord);
                                        guidata1.drawing_idx(x_coord,y_coord,z_coord) = guidata1.drawing_draw_color;
                                        guidata1.curr_drawing(x_coord,y_coord,z_coord) = guidata1.drawing_draw_color;
                                    end
                                end
                            end
                        end
                        idx_out = find(guidata1.curr_drawing(:)>0);
                        [ix,iy,iz] = ind2sub(size(guidata1.curr_drawing),idx_out);
                        guidata1.last_draw = [ix,iy,iz];
                        guidata(handles.figure,guidata1);
                        calc_full_ind
                        update_image;
                    end
                case 'overlay 1'
                    guidata1.last_actions(2) = 3; % overlay 1 edit
                    guidata1.prechange2_overlay1 = guidata1.prechange1_overlay1;
                    guidata1.prechange1_overlay1 = guidata1.overlay1_imgfull; % update cache
                    guidata1.overlay1_imgfull(guidata1.xind(pt(2)),guidata1.yind(pt(1)),guidata1.curr_slice) = guidata1.overlay1_draw_color; % update data
                case 'overlay 2'
                    guidata1.last_actions(2) = 4; % overlay 2 edit
                    guidata1.prechange2_overlay2 = guidata1.prechange1_overlay2;
                    guidata1.prechange1_overlay2 = guidata1.overlay2_imgfull; % update cache
                    guidata1.overlay2_imgfull(guidata1.xind(pt(2)),guidata1.yind(pt(1)),guidata1.curr_slice) = guidata1.overlay2_draw_color; % update data                    
            end
            guidata1.undo_num = 0;
            guidata1.save_points = [pt(2),pt(1)];
            guidata1.in_motion = false;
            guidata(handles.figure,guidata1);
        end
    end
end

% Mouse Motion after Click
function motion_callback(~,~,~)
    % Load Data:
    guidata1 = guidata(handles.figure);
    set(handles.figure, 'WindowButtonUpFcn', @unclick_callback);
    % Get mouse location
    pt = get(handles.background_axes,'CurrentPoint');
    pt = round(pt(1,1:2));
    % Add color:
    pt_check = sum([(pt>0),(pt<=[length(guidata1.yind),length(guidata1.xind)])]);
    if pt_check~=4; return; end
    guidata1 = guidata(handles.figure);
    guidata1.in_motion = true;
    switch guidata1.currently_editing
        case 'background'
            guidata1.imgfull(guidata1.xind(pt(2)),guidata1.yind(pt(1)),guidata1.curr_slice) = guidata1.background_draw_color;
        case 'drawing'
            if guidata1.shape==0
                guidata1.drawing_idx(guidata1.xind(pt(2)),guidata1.yind(pt(1)),guidata1.curr_slice) = guidata1.drawing_draw_color;
            elseif guidata1.shape<3
                % Erase Previous Shape if Still Dragging
                if ~isempty(guidata1.prev_shape_coords)
                    for ix1 = 1:size(guidata1.prev_shape_coords,1)
                    guidata1.drawing_idx(guidata1.prev_shape_coords(ix1,1),...
                        guidata1.prev_shape_coords(ix1,2),guidata1.curr_slice) = guidata1.shape_prev_colors(ix1);
                    end
                end
                guidata1.shape_ind(2,:) = [guidata1.xind(pt(2)),guidata1.yind(pt(1))];
                if guidata1.shape==1 % circle
                    center = sum(guidata1.shape_ind,1).*.5;
                    radius = .5*sqrt(sum(diff(guidata1.shape_ind).^2));
                  % % for circle outline
%                         angles = linspace(0.0175,2*pi,64);
%                         guidata1.shape_coords = zeros(64,2);
%                         guidata1.shape_coords(:,1) = max(round(center(1) + radius.*sin(angles)),1); 
%                         guidata1.shape_coords(:,2) = max(round(center(2) + radius.*cos(angles)),1);
                    [rr,cc] = meshgrid(1:guidata1.ydim,1:guidata1.xdim);
                    C = sqrt((rr-center(2)).^2+(cc-center(1)).^2)<=radius;
                    [ix, iy] = ind2sub(guidata1.dim([3,2]),find(C(:)));
                    guidata1.shape_coords = [ix,iy];
                elseif guidata1.shape==2 % rectangle
                    min_x = min(guidata1.shape_ind(:,1)); min_y = min(guidata1.shape_ind(:,2)); 
                    max_x = max(guidata1.shape_ind(:,1)); max_y = max(guidata1.shape_ind(:,2)); 
                    list_x_n = max_x-min_x+1; list_y_n = max_y-min_y+1;
                    guidata1.shape_coords = zeros(list_x_n*list_y_n,2);
                    count = 0;
                    for ix = min_x:max_x
                        for iy = min_y:max_y
                            count = count+1;
                            guidata1.shape_coords(count,:) = [ix,iy];
                        end
                    end
                    % % for outline:
%                         list_x = (min_x:max_x)'; list_y = (min_y:max_y)';
%                         guidata1.shape_coords = [list_x,repmat(guidata1.shape_ind(1,2),list_x_n,1);...
%                             list_x,repmat(guidata1.shape_ind(2,2),list_x_n,1);...
%                             repmat(guidata1.shape_ind(1,1),list_y_n,1),list_y;...
%                             repmat(guidata1.shape_ind(2,1),list_y_n,1),list_y];
                end
                guidata1.shape_coords = unique(guidata1.shape_coords,'rows');
                guidata1.shape_coords(((guidata1.shape_coords(:,1)>max(guidata1.xind))+(guidata1.shape_coords(:,2)>max(guidata1.yind))>0),:) = [];
                guidata1.prev_shape_coords = guidata1.shape_coords;
                % Draw New Shape
                guidata1.curr_drawing = zeros(size(guidata1.drawing_idx)); % used for propagating last drawing
                for ix1 = 1:size(guidata1.shape_coords,1)
                    guidata1.shape_prev_colors(ix1) = guidata1.drawing_idx(guidata1.shape_coords(ix1,1),...
                        guidata1.shape_coords(ix1,2),guidata1.curr_slice);
                    guidata1.drawing_idx(guidata1.shape_coords(ix1,1),...
                        guidata1.shape_coords(ix1,2),guidata1.curr_slice) = guidata1.drawing_draw_color;
                    guidata1.curr_drawing(guidata1.shape_coords(ix1,1),...
                        guidata1.shape_coords(ix1,2),guidata1.curr_slice) = guidata1.drawing_draw_color;
                end
                idx_out = find(guidata1.curr_drawing(:)>0);
                [ix,iy,iz] = ind2sub(size(guidata1.curr_drawing),idx_out);
                guidata1.last_draw = [ix,iy,iz];
                guidata(handles.figure,guidata1);
                calc_full_ind
                update_image;
                return;
            end
        case 'overlay 1'
            guidata1.overlay1_imgfull(guidata1.xind(pt(2)),guidata1.yind(pt(1)),guidata1.curr_slice) = guidata1.overlay1_draw_color;
        case 'overlay 2'
            guidata1.overlay2_imgfull(guidata1.xind(pt(2)),guidata1.yind(pt(1)),guidata1.curr_slice) = guidata1.overlay2_draw_color;
    end
    guidata1.save_points = [guidata1.save_points; pt(2),pt(1)];
    guidata(handles.figure,guidata1);
%         update_image % this turns on tracing of drawings, but decreases performance
end

% Cursor Unclick Callback
function unclick_callback(~,~,~)
    guidata1 = guidata(handles.figure);
    set(handles.figure, 'WindowButtonMotionFcn', []);
    set(handles.figure, 'WindowButtonUpFcn', []);
    if guidata1.shape==0
        if guidata1.in_motion
            guidata1.save_points = unique(guidata1.save_points,'rows');
            if guidata1.shape>0
                guidata1.save_points = guidata1.shape_coords;
            end
            guidata1.in_motion = false;
            guidata1.prev_shape_coords = [];
            if size(guidata1.save_points,1) > 2
                drag_interpolate;
            end
        else
            calc_full_ind
            update_image
        end
    else
        guidata1.in_motion = false;
        guidata1.prev_shape_coords = [];
    end
    guidata(handles.figure,guidata1);
end

% Apply Mouse Motion Data after Unclick:
function drag_interpolate
    guidata1 = guidata(handles.figure);
    % Fill-In Area:
    inside = inpolygon(guidata1.poss_ind(:,1),guidata1.poss_ind(:,2),guidata1.save_points(:,1),guidata1.save_points(:,2));
    guidata1.save_points = unique([guidata1.save_points(:,[1,2]); guidata1.poss_ind(inside,:)],'rows');
    guidata1.save_points(((guidata1.save_points(:,1)>max(guidata1.xind))+(guidata1.save_points(:,2)>max(guidata1.yind)))>0,:) = [];
    switch guidata1.currently_editing
        case 'background'
            for i = 1:size(guidata1.save_points,1)
                guidata1.imgfull(guidata1.xind(guidata1.save_points(i,1)),...
                    guidata1.yind(guidata1.save_points(i,2)),guidata1.curr_slice) = guidata1.background_draw_color;
            end
        case 'drawing'
            guidata1.curr_drawing = zeros(size(guidata1.drawing_idx)); % used for propagating last drawing
            for i = 1:size(guidata1.save_points,1)
                guidata1.drawing_idx(guidata1.xind(guidata1.save_points(i,1)),...
                    guidata1.yind(guidata1.save_points(i,2)),guidata1.curr_slice) = guidata1.drawing_draw_color;
                guidata1.curr_drawing(guidata1.xind(guidata1.save_points(i,1)),...
                    guidata1.yind(guidata1.save_points(i,2)),guidata1.curr_slice) = guidata1.drawing_draw_color;
            end
            % Save Coordinates of Last Drawing
            idx_out = find(guidata1.curr_drawing(:)>0);
            [ix,iy,iz] = ind2sub(size(guidata1.curr_drawing),idx_out);
            guidata1.last_draw = [ix,iy,iz];
        case 'overlay 1'
            for i = 1:size(guidata1.save_points,1)
                guidata1.overlay1_imgfull(guidata1.xind(guidata1.save_points(i,1)),...
                    guidata1.yind(guidata1.save_points(i,2)),guidata1.curr_slice) = guidata1.overlay1_draw_color;
            end
        case 'overlay 2'
            for i = 1:size(guidata1.save_points,1)
                guidata1.overlay2_imgfull(guidata1.xind(guidata1.save_points(i,1)),...
                    guidata1.yind(guidata1.save_points(i,2)),guidata1.curr_slice) = guidata1.overlay2_draw_color;
            end            
    end
    guidata(handles.figure,guidata1);
    calc_full_ind
    update_image
end

function scroll_zoom_callback(~, eventdata, ~)
    % Load Data:
    guidata1 = guidata(handles.figure);
    % Get mouse location
    pt = get(handles.background_axes,'CurrentPoint');
    pt = round(pt(1,1:2));
    % Check if Image Click
    if pt(2) <= guidata1.xdim && pt(2) > 0 && pt(1) <= guidata1.ydim && pt(1) > 0
        guidata1.scroll_count = min(max(guidata1.scroll_count + -eventdata.VerticalScrollCount,0),9);
        if guidata1.scroll_count ~= 0
            zoom_factor = guidata1.scroll_zoom_equiv(guidata1.scroll_count+1);
            xlength = zoom_factor*guidata1.xdim; ylength = zoom_factor*guidata1.ydim;
            guidata1.xmax = min(round(guidata1.xind(pt(2))+.5*xlength),guidata1.xdim);
            guidata1.xmin = max(round(guidata1.xind(pt(2))-.5*xlength),1);
            guidata1.ymax = min(round(guidata1.yind(pt(1))+.5*ylength),guidata1.ydim);
            guidata1.ymin = max(round(guidata1.yind(pt(1))-.5*ylength),1);
            guidata1.xind = guidata1.xmax:-1:guidata1.xmin;
            if guidata1.slice_orientation==2
                guidata1.yind = guidata1.ymax:-1:guidata1.ymin;
            else
                guidata1.yind = guidata1.ymin:guidata1.ymax;
            end
%             guidata1.xind = guidata1.xmin:guidata1.xmax;
%             guidata1.yind = guidata1.ymin:guidata1.ymax;
        else
            guidata1.xind = guidata1.xdim:-1:1;
            if guidata1.slice_orientation==2
                guidata1.yind = guidata1.ydim:-1:1;
            else
                guidata1.yind = 1:guidata1.ydim;
            end
%             guidata1.xind = 1:guidata1.xdim;
%             guidata1.yind = 1:guidata1.ydim;
        end
        guidata1.ax_xlim = [1,numel(guidata1.yind)];
        guidata1.ax_ylim = [1,numel(guidata1.xind)];
        guidata1.refresh_img = true;
        guidata(handles.figure,guidata1);
        update_image
        reposition_axes
    end
end

%% Overlay Callbacks:

% Add/Edit Overlay Callback (Open Overlay GUI):
function add_overlay_callback(~,~,~)
    guidata1 = guidata(handles.figure);
    % Open Overlay Window (auto-detect screen size):
    if ~ishandle(guidata1.overlay_figure)
        guidata1.overlay_figure = figure('Position',guidata1.overlay_figure_pos,...
            'MenuBar','none','Name','Add/Edit Overlays','NumberTitle','off',...
            'Color',[1,1,1],'Units','Normalized','doublebuffer','on');
    else
        figure(guidata1.overlay_figure);
    end
        % Create Buttons:
    % Create "Open Overlay" Button:
    guidata1.open_overlay_button1 = uicontrol(guidata1.overlay_figure,'style',...
        'pushbutton','Units','normalized','Position',[.008,.66,.08,.2],...
        'BackgroundColor', [0,0,1],'FontWeight','bold','String','Open',...
        'callback',{@open_overlay_callback,1}); 
    % Create "Open Overlay" Button 2: 
    guidata1.open_overlay_button1 = uicontrol(guidata1.overlay_figure,'style',...
        'pushbutton','Units','normalized','Position',[.008,.23,.08,.2],...
        'BackgroundColor', [0,0,1],'FontWeight','bold','String','Open',...
        'callback',{@open_overlay_callback,2}); 
    % Create "Close Overlay" button
    guidata1.close_overlay_button1 = uicontrol(guidata1.overlay_figure,'style',...
        'pushbutton','Units','normalized','Position',[.088,.66,.08,.2],...
        'BackgroundColor', [1,0,0],'FontWeight','bold','String','Close',...
        'callback',{@close_overlay_callback,1}); 
    % Create "Close Overlay" button 2
    guidata1.close_overlay_button2 = uicontrol(guidata1.overlay_figure,'style',...
        'pushbutton','Units','normalized','Position',[.088,.23,.08,.2],...
        'BackgroundColor', [1,0,0],'FontWeight','bold','String','Close',...
        'callback',{@close_overlay_callback,2}); 
        % Create "Save Overlay" button
    guidata1.save_overlay_button1 = uicontrol(guidata1.overlay_figure,'style',...
        'pushbutton','Units','normalized','Position',[.17,.66,.08,.2],...
        'BackgroundColor', [0,1,0],'FontWeight','bold','String','Save',...
        'callback',{@save_overlay_callback,1}); 
    % Create "Save Overlay" button 2
    guidata1.save_overlay_button2 = uicontrol(guidata1.overlay_figure,'style',...
        'pushbutton','Units','normalized','Position',[.17,.23,.08,.2],...
        'BackgroundColor', [0,1,0],'FontWeight','bold','String','Save',...
        'callback',{@save_overlay_callback,2}); 
    % Create "Colormap:" Textbox
    uicontrol('Parent',guidata1.overlay_figure,'Style','text','String','Colormap:',...
        'Units', 'normalized','FontSize',9,'FontName','Helvetica','Position',...
        [.26, .61, .1, .2],'BackgroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','center');
    % Create "Colormap:" Textbox 2
    uicontrol('Parent',guidata1.overlay_figure,'Style','text','String','Colormap:',...
        'Units', 'normalized','FontSize',9,'FontName','Helvetica','Position',...
        [.26, .185, .1, .2],'BackgroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','center');
    % Create Colormap Selection Drop-Down:
    guidata1.overlay1_colormap_dropdown = uicontrol('Parent',guidata1.overlay_figure,'Style','popupmenu',...
        'String','jet|hot|cool|spring|summer|winter',...
        'Units', 'normalized','FontSize',9,'FontName','Helvetica','Position',...
        [.37, .64, .09, .2],'BackgroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','center','callback',{@overlay_colormap_callback,1},...
        'Value',guidata1.overlay1_colormap_val);
    % Create Colormap Selection Drop-Down 2:
    guidata1.overlay2_colormap_dropdown = uicontrol('Parent',guidata1.overlay_figure,'Style','popupmenu',...
        'String','jet|hot|cool|spring|summer|winter',...
        'Units', 'normalized','FontSize',9,'FontName','Helvetica','Position',...
        [.37, .22, .09, .2],'BackgroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','center','callback',{@overlay_colormap_callback,2},...
        'Value',guidata1.overlay2_colormap_val);
    % Create Opacity Selection Drop-Down:
    guidata1.opacity1 = uicontrol('Parent',guidata1.overlay_figure,'Style','popupmenu',...
        'String','Opaque|90%|80%|70%|60%|50%|40%|30%|20%|10%|Invisible',...
        'Units','normalized','FontSize',9,'FontName','Helvetica','Position',...
        [.475, .64, .11, .2],'BackgroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','center','callback',{@overlay_opacity_callback,1},...
        'Value',guidata1.transparency_val1);
    % Create Opacity Selection Drop-Down 2:
    guidata1.opacity2 = uicontrol('Parent',guidata1.overlay_figure,'Style','popupmenu',...
        'String','Opaque|90%|80%|70%|60%|50%|40%|30%|20%|10%|Invisible',...
        'Units', 'normalized','FontSize',9,'FontName','Helvetica','Position',...
        [.475, .22, .11, .2],'BackgroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','center','callback',{@overlay_opacity_callback,2},...
        'Value',guidata1.transparency_val2);
    % Create "Color Axis:" Textbox
    uicontrol('Parent',guidata1.overlay_figure,'Style','text','String','Color Axis:',...
        'Units', 'normalized','FontSize',9,'FontName','Helvetica','Position',...
        [.595, .61, .1, .2],'BackgroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','center');
    % Create "Color Axis:" Textbox 2
    uicontrol('Parent',guidata1.overlay_figure,'Style','text','String','Color Axis:',...
        'Units', 'normalized','FontSize',9,'FontName','Helvetica','Position',...
        [.595, .185, .1, .2],'BackgroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','center');
    % Create CMIN Entry Box:
    guidata1.overlay1_cmin = uicontrol(guidata1.overlay_figure,'style','edit',...
        'Units','normalized','Position',[.7,.67,.04,.16],'String',guidata1.cmin_overlay1,...
        'BackgroundColor',[1,1,1],'callback',@overlay1_cmin_callback); 
    % Create CMIN Entry Box 2:
    guidata1.overlay2_cmin = uicontrol(guidata1.overlay_figure,'style','edit',...
        'Units','normalized','Position',[.7,.25,.04,.16],'String',guidata1.cmin_overlay2,...
        'BackgroundColor',[1,1,1],'callback',@overlay2_cmin_callback); 
    % Create " - " Textbox
    uicontrol('Parent',guidata1.overlay_figure,'Style','text','String',' - ',...
        'Units', 'normalized','FontSize',12,'FontName','Helvetica','Position',...
        [.741,.64,.02,.2],'BackgroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','center');
    % Create " - " Textbox 2
    uicontrol('Parent',guidata1.overlay_figure,'Style','text','String',' - ',...
        'Units', 'normalized','FontSize',12,'FontName','Helvetica','Position',...
        [.741,.215,.02,.2],'BackgroundColor',[1,1,1],'FontWeight','Bold',...
        'HorizontalAlignment','center');
    % Create CMAX Entry Box:
    guidata1.overlay1_cmax = uicontrol(guidata1.overlay_figure,'style','edit',...
        'Units','normalized','Position',[.765,.67,.04,.16],'String',guidata1.cmax_overlay1,...
        'BackgroundColor',[1,1,1],'callback',@overlay1_cmax_callback); 
    % Create CMAX Entry Box 2:
    guidata1.overlay2_cmax = uicontrol(guidata1.overlay_figure,'style','edit',...
        'Units','normalized','Position',[.765,.25,.04,.16],'String',guidata1.cmax_overlay2,...
        'BackgroundColor',[1,1,1],'callback',@overlay2_cmax_callback); 
    % Overlay Colorbar 1:
    guidata1.overlay1_colorbar = uicontrol(guidata1.overlay_figure,'style',...
        'togglebutton','Units','normalized','Position',[.823,.66,.08,.2],...
        'BackgroundColor', [1,0,1],'FontWeight','bold','String','Colorbar',...
        'callback',@overlay1_colorbar_callback); 
    % Overlay Colorbar 2:
    guidata1.overlay2_colorbar = uicontrol(guidata1.overlay_figure,'style',...
        'togglebutton','Units','normalized','Position',[.823,.227,.08,.2],...
        'BackgroundColor', [1,0,1],'FontWeight','bold','String','Colorbar',...
        'callback',@overlay2_colorbar_callback); 
    % Edit Overlay 1:
    guidata1.overlay1_edit_button = uicontrol(guidata1.overlay_figure,'style',...
        'togglebutton','Units','normalized','Position',[.92,.66,.065,.2],...
        'BackgroundColor', [0,1,1],'FontWeight','bold','String','Edit',...
        'callback',{@overlay_edit_callback,1}); 
    % Edit Overlay 2:
    guidata1.overlay2_edit_button = uicontrol(guidata1.overlay_figure,'style',...
        'togglebutton','Units','normalized','Position',[.92,.227,.065,.2],...
        'BackgroundColor', [0,1,1],'FontWeight','bold','String','Edit',...
        'callback',{@overlay_edit_callback,2}); 
    guidata(handles.figure,guidata1);
end

function open_overlay_callback(hObject, ~, overlay_num, load_type)
    guidata1 = guidata(handles.figure);
    % Before update, go to default orientation:
    switch guidata1.slice_orientation
        case 1
            save_orientation = guidata1.h_x;
            reorient_callback(guidata1.h_z) % revert to original orientation
            change_orientation = 1;
        case 2
            save_orientation = guidata1.h_y;
            reorient_callback(guidata1.h_z) % revert to original orientation
            change_orientation = 1;
        case 3 
            change_orientation = 0;
    end     
    if nargin<4; load_type = 'char'; end
    switch load_type
        case 'char'
            if nargin < 4
                % Get input filename:
                cd(guidata1.last_nav_dir);
                [overlay_name, overlay_path] = uigetfile(guidata1.extensions,'Select Overlay Image:','MultiSelect','off');
                cd(guidata1.first_dir);
            elseif nargin==4
                if overlay_num==1
                    [overlay_path,overlay_name,ext] = fileparts(guidata1.parsed_inputs.(guidata1.poss_input{2}));
                elseif overlay_num==2
                    [overlay_path,overlay_name,ext] = fileparts(guidata1.parsed_inputs.(guidata1.poss_input{3}));
                end
                overlay_name = [overlay_name,ext];  
            end
            if ischar(overlay_name)
                fname1 = fullfile(overlay_path,overlay_name); 
                % Remember Extension Chosen:
                if ~isempty(overlay_path)
                    sort_exts(overlay_path,overlay_name)
                end                    
                if guidata1.untouch_nii
                    if overlay_num==1
                        guidata1.overlay1_img = load_untouch_nii(fname1);
                        guidata1.overlay1_on = true;
                    elseif overlay_num==2
                        guidata1.overlay2_img = load_untouch_nii(fname1);
                        guidata1.overlay2_on = true;
                    end
                else
                    if overlay_num==1
                        guidata1.overlay1_img = load_nii(fname1);
                        guidata1.overlay1_on = true;
                    elseif overlay_num==2
                        guidata1.overlay2_img = load_nii(fname1);
                        guidata1.overlay2_on = true;
                    end
                end
            else
                errordlg('ERROR: Could not find overlay file.','ERROR')
                return;   
            end
        case 'struct'
            if overlay_num==1
                try
                    guidata1.overlay1_img = guidata1.parsed_inputs.(guidata1.poss_input{2});
                    guidata1.overlay1_on = true;
                catch
                    disp('Error loading overlay 1 image structure. Check fields and retry or load through GUI.')
                    return;
                end
            elseif overlay_num==2
                try
                    guidata1.overlay2_img = guidata1.parsed_inputs.(guidata1.poss_input{3});
                    guidata1.overlay2_on = true;
                catch
                    disp('Error loading overlay 2 image structure. Check fields and retry or load through GUI.')
                    return;
                end
            end
        case 'matrix'
            if overlay_num==1
                guidata1.overlay1_img = guidata1.img;
                guidata1.overlay1_img.img = guidata1.parsed_inputs.(guidata1.poss_input{2});
                guidata1.overlay1_on = true;
            elseif overlay_num==2
                guidata1.overlay2_img = guidata1.img;
                guidata1.overlay2_img.img = guidata1.parsed_inputs.(guidata1.poss_input{3});
                guidata1.overlay2_on = true;
            end
    end
    % Parse image input
    if overlay_num==1
        guidata1.overlay1_imgfull = single(guidata1.overlay1_img.img); % conversion may enhance refresh speed
        if ~all(size(guidata1.overlay1_imgfull)==guidata1.dim(2:4))
            disp('Error: Overlay dimensions do not match original image.');
            return;
        end
        guidata1.overlay1_imgfull = permute(guidata1.overlay1_imgfull,[2,1,3]);
        % Calculate Image Color Indices (Intialize to caxis([min,max]))
        guidata1.cmin_overlay1 = min(guidata1.overlay1_imgfull(:));
        guidata1.cmax_overlay1 = max(guidata1.overlay1_imgfull(:));
    elseif overlay_num==2
        % Parse image input
        guidata1.overlay2_imgfull = single(guidata1.overlay2_img.img); % conversion may enhance refresh speed
        if ~all(size(guidata1.overlay2_imgfull)==guidata1.dim(2:4))
            disp('Error: Overlay dimensions do not match original image.');
            return;
        end
        guidata1.overlay2_imgfull = permute(guidata1.overlay2_imgfull,[2,1,3]);
        % Calculate Image Color Indices (Intialize to caxis([min,max]))
        guidata1.cmin_overlay2 = min(guidata1.overlay2_imgfull(:));
        guidata1.cmax_overlay2 = max(guidata1.overlay2_imgfull(:));                
    end
    % Update Figure:
    guidata(handles.figure,guidata1);
    calc_full_ind
    update_image % overlays are initialized here
    % Return to default colormap:
    if nargin<4
        add_overlay_callback % re-calling function updates cmin,cmax in add/edit overlays window
        if overlay_num==1
            overlay_colormap_callback(hObject,[],1) 
        elseif overlay_num==2
            overlay_colormap_callback(hObject,[],2)
        end
    end
    % Return to Original Orientation if Changed:
    if change_orientation
        reorient_callback(save_orientation,[]) % have to have two inputs, see line 1047
    end
end

function close_overlay_callback(~, ~, overlay_num) 
    guidata1 = guidata(handles.figure);
    if overlay_num==1
        guidata1.overlay1_on = false;
        guidata1.cmin_overlay1 = [];
        guidata1.cmax_overlay1 = [];
        if ishandle(handles.overlay1_img)
            delete(handles.overlay1_img)
        end
        guidata1.overlay1_colorbar_on = false;
        if ~guidata1.overlay2_on
            guidata1.overlay2_colorbar_on = false;
            guidata1.main_colorbar_on = true;
        end
    elseif overlay_num==2
        guidata1.overlay2_on = false;
        guidata1.cmin_overlay2 = [];
        guidata1.cmax_overlay2 = [];
        if ishandle(handles.overlay2_img)
            delete(handles.overlay2_img)
        end
        guidata1.overlay2_colorbar_on = false;
        if ~guidata1.overlay1_on
            guidata1.overlay1_colorbar_on = false;
            guidata1.main_colorbar_on = true;
        end        
    end
    guidata1.update_colorbar = true;
    guidata(handles.figure,guidata1);
    update_image
    add_overlay_callback
end

function save_overlay_callback(~, ~, overlay_num)
    guidata1 = guidata(handles.figure);
    if guidata1.overlay1_on && overlay_num==1
        save_num = 1;
    elseif guidata1.overlay2_on && overlay_num==2
        save_num = 2;
    else
        return;
    end
    % Get saveas filename:
    cd(guidata1.last_nav_dir);
    [overlay_name, overlay_path] = uiputfile(guidata1.extensions,'Specify overlay filename:');
    cd(guidata1.first_dir);
    if any(overlay_name == 0) % if cancel, or X
        return;
    end
    if ~ischar(overlay_name) 
        disp('Please enter a valid filename.')
        return
    end
    sort_exts(overlay_path,overlay_name)
    if guidata1.slice_orientation~=3
        reorient_callback(guidata1.h_z) % revert to original orientation
    end
    % Save Image File:
    fname = fullfile(overlay_path,overlay_name); 
    if save_num == 1
        guidata1.overlay1_img.img = permute(guidata1.overlay1_imgfull,[2,1,3]);
        if guidata1.untouch_nii
            save_untouch_nii(guidata1.overlay1_img, fname);
        else
            save_nii(guidata1.overlay1_img, fname);
        end
    elseif save_num == 2
        guidata1.overlay2_img.img = permute(guidata1.overlay2_imgfull,[2,1,3]);
        if guidata1.untouch_nii
            save_untouch_nii(guidata1.overlay2_img, fname);
        else
            save_nii(guidata1.overlay2_img, fname);
        end
    end
    disp(['Overlay successfully saved as ',fname])
    guidata(handles.figure,guidata1);
end

function overlay_colormap_callback(hObject, ~, overlay_num) 
    guidata1 = guidata(handles.figure);
    map_num = get(hObject,'Value');
    if overlay_num==1 && guidata1.overlay1_on
        guidata1.overlay1_colormap_val = map_num;
        guidata1.overlay1_colormap_selection = guidata1.overlay_colormaps{map_num};
    elseif overlay_num==2 && guidata1.overlay2_on
        guidata1.overlay2_colormap_val = map_num;
        guidata1.overlay2_colormap_selection = guidata1.overlay_colormaps{map_num};
    else
        return;
    end
    guidata(handles.figure,guidata1);
    update_colormap
    calc_full_ind
    update_image
end

function overlay_opacity_callback(hObject, ~, overlay_num) 
    guidata1 = guidata(handles.figure);
    transparency_val = get(hObject,'Value');
    switch transparency_val
        case 1; overlay_alpha = 1;
        case 2; overlay_alpha = .9;
        case 3; overlay_alpha = .8;
        case 4; overlay_alpha = .7;    
        case 5; overlay_alpha = .6;
        case 6; overlay_alpha = .5;
        case 7; overlay_alpha = .4;
        case 8; overlay_alpha = .3;
        case 9; overlay_alpha = .2;
        case 10; overlay_alpha = .1;
        case 11; overlay_alpha = 0;
    end
    if overlay_num==1 && guidata1.overlay1_on
        guidata1.transparency_val1 = transparency_val;
        guidata1.overlay1_alpha = overlay_alpha;
    elseif overlay_num==2 && guidata1.overlay2_on
        guidata1.transparency_val2 = transparency_val;
        guidata1.overlay2_alpha = overlay_alpha;
    else
        return;
    end
    guidata(handles.figure,guidata1);
    calc_full_ind
    update_image
end

function overlay1_cmin_callback(hObject, ~, ~)
    guidata1 = guidata(handles.figure);
    if guidata1.overlay1_on    
        val = str2double(get(hObject,'String'));
        if isnumeric(val) && val < guidata1.cmax_overlay1
            guidata1.cmin_overlay1 = val;
        else
            errordlg('Enter a numeric value less than cmax.');
        end
        guidata(handles.figure,guidata1);
        calc_full_ind
        update_image
    end
end

function overlay2_cmin_callback(hObject, ~, ~)
    guidata1 = guidata(handles.figure);
    if guidata1.overlay2_on
        val = str2double(get(hObject,'String'));
        if isnumeric(val) && val < guidata1.cmax_overlay2
            guidata1.cmin_overlay2 = val;
        else
            errordlg('Enter a numeric value less than cmax.');
        end
        guidata(handles.figure,guidata1);
        calc_full_ind
        update_image
    end
end

function overlay1_cmax_callback(hObject, ~, ~)
    guidata1 = guidata(handles.figure);
    if guidata1.overlay1_on
        val = str2double(get(hObject,'String'));
        if isnumeric(val) && val > guidata1.cmin_overlay1
            guidata1.cmax_overlay1 = val;
        else
            errordlg('Enter a numeric value greater than cmin.');
        end
        guidata(handles.figure,guidata1);
        calc_full_ind
        update_image
    end
end

function overlay2_cmax_callback(hObject, ~, ~)
    guidata1 = guidata(handles.figure);
    if guidata1.overlay2_on    
        val = str2double(get(hObject,'String'));
        if isnumeric(val) && val > guidata1.cmin_overlay2
            guidata1.cmax_overlay2 = val;
        else
            errordlg('Enter a numeric value greater than cmin.');
        end
        guidata(handles.figure,guidata1);
        calc_full_ind
        update_image
    end
end

% Add Overlay 1 Colorbar:
function overlay1_colorbar_callback(hObject, calledFromInitial, ~)
    guidata1 = guidata(handles.figure);
    if ~isempty(calledFromInitial)
        value = 1;
    else
        value = get(hObject,'Value');
    end
    if value
        if guidata1.overlay1_on
            guidata1.overlay1_colorbar_on = true;
            guidata1.main_colorbar_on = false;
        else
            guidata1.main_colorbar_on = true;
        end
        guidata1.refresh_img = true;
        guidata1.overlay2_colorbar_on = false;
        try
            set(guidata1.overlay2_colorbar,'Value',false);
        catch
        end
    else
        if guidata1.overlay1_on
            guidata1.refresh_img = true;
        end
        guidata1.overlay1_colorbar_on = false;
        if ~guidata1.overlay2_colorbar_on
            guidata1.main_colorbar_on = true;
            guidata1.refresh_img = true;
        end
    end
    guidata1.update_colorbar = true;
    guidata(handles.figure,guidata1);
    update_image
end

% Add Overlay 2 Colorbar:
function overlay2_colorbar_callback(hObject, calledFromInitial, ~)
    guidata1 = guidata(handles.figure);
    if ~isempty(calledFromInitial)
        value = 1;
    else
        value = get(hObject,'Value');
    end
    if value
        if guidata1.overlay2_on
            guidata1.overlay2_colorbar_on = true;
            guidata1.main_colorbar_on = false;
        else
            guidata1.main_colorbar_on = true;
        end
        guidata1.refresh_img = true;
        guidata1.overlay1_colorbar_on = false;
        try
            set(guidata1.overlay1_colorbar,'Value',false);
        catch
        end
    else
        if guidata1.overlay2_on
            guidata1.refresh_img = true;
        end
        guidata1.overlay2_colorbar_on = false;
        if ~guidata1.overlay1_colorbar_on
            guidata1.main_colorbar_on = true;
            guidata1.refresh_img = true;
        end
    end
    guidata1.update_colorbar = true;
    guidata(handles.figure,guidata1);
    update_image
end

function overlay_edit_callback(hObject, ~, overlay_num)
    guidata1 = guidata(handles.figure);
    if overlay_num==1 && guidata1.overlay1_on
        if get(hObject,'Value')
            guidata1.previously_editing = guidata1.currently_editing;
            guidata1.currently_editing = 'overlay 1';
            set(guidata1.overlay2_edit_button,'Value',false);
        else % turn off
            guidata1.currently_editing = guidata1.previously_editing;            
        end
    elseif overlay_num==2 && guidata1.overlay2_on
        if get(hObject,'Value')
            guidata1.previously_editing = guidata1.currently_editing;
            guidata1.currently_editing = 'overlay 2';
            set(guidata1.overlay1_edit_button,'Value',false);
        else % turn off
            guidata1.currently_editing = guidata1.previously_editing;
        end
    else % invalid call
        if overlay_num==1
            set(guidata1.overlay1_edit_button,'Value',false);
        elseif overlay_num==2
            set(guidata1.overlay2_edit_button,'Value',false);
        end
        return;
    end
    % Get user input for draw color:
    if get(hObject,'Value')
        prompt = {'Enter Overlay Draw Color: '};
        dlg_title = 'Color Spec'; num_lines = [1,30];
        defaultans = {'0'};
        answer1 = inputdlg(prompt,dlg_title,num_lines,defaultans,'on');
        if ~isempty(answer1)
            draw_color = str2double(answer1);
            if isnan(draw_color)
              errordlg('You must enter a numeric value',...
                  'Invalid Input','modal')
              uicontrol(hObject)
              return 
            else  
                switch guidata1.currently_editing
                    case 'overlay 1'
                        guidata1.overlay1_draw_color = draw_color;
                    case 'overlay 2'
                        guidata1.overlay2_draw_color = draw_color;
                end
            end
        end
    end
    guidata(handles.figure,guidata1)  ;
end
    
%% Drawing Overlay Functions

% Select Color for Overlay Drawing:
function draw_overlay_callback(hObject, ~, ~)
    guidata1 = guidata(handles.figure);
    guidata1.drawing_active = true;
    guidata1.currently_editing = 'drawing';
    guidata1.refresh_img = true;
    color_select = get(hObject,'Label');
    switch color_select
        case 'Erase'; guidata1.drawing_draw_color = 1;
        case 'Signal ROI (Green)'; guidata1.drawing_draw_color = 2;
        case 'Noise ROI (Red)'; guidata1.drawing_draw_color = 3;
    end
    if isempty(guidata1.drawing_idx) && guidata1.drawing_active
        guidata1.drawing_idx = ones(size(guidata1.imgfull));
    end
    guidata(handles.figure,guidata1);
end

% Determine transparency (alpha) for drawing:
function draw_transparency_callback(hObject, ~, ~)
    guidata1 = guidata(handles.figure);
    transparency_val = get(hObject,'Label');
    switch transparency_val
        case 'Opaque (default)'; guidata1.overlay_alpha = 1; 
        case '80%'; guidata1.overlay_alpha = .8;
        case '60%'; guidata1.overlay_alpha = .6;
        case '40%'; guidata1.overlay_alpha = .4;
        case '20%'; guidata1.overlay_alpha = .2;
        case 'Invisible'; guidata1.overlay_alpha = 0;
    end
    guidata(handles.figure,guidata1);
    if guidata1.drawing_active || guidata1.overlay1_on || guidata1.overlay2_on
        calc_full_ind
        update_image
    end
end

function draw_shapes_callback(~, ~, type)
    guidata1 = guidata(handles.figure);
    switch type
        case 1 % no shape
            guidata1.shape = 0;
        case 2 % circle
            guidata1.shape = 1;
        case 3 % rectangle
            guidata1.shape = 2;
        case 4 % sphere
            prompt = {'Specify sphere radius in physical units of image:'};
            dlg_title = ''; num_lines = [1,25]; 
            answer = inputdlg(prompt,dlg_title,num_lines);
            if isempty(answer); return; end
            guidata1.sphere_radius = str2double(answer(1));
            guidata1.shape = 3;
    end
    guidata(handles.figure,guidata1);
end

function propagate_draw_callback(~, ~, ~)
    guidata1 = guidata(handles.figure);
    if guidata1.drawing_active && ~isempty(guidata1.curr_drawing)
        prompt = {sprintf('Slice Range: \n\nFirst:'),'Last:'};
        dlg_title = ''; num_lines = [1,15;1,15]; 
        answer = inputdlg(prompt,dlg_title,num_lines);
        if isempty(answer); return; end
        slice_first = str2double(answer(1)); slice_last = str2double(answer(2));
        ind = find(guidata1.curr_drawing(:)>0);
        [ix,iy,~] = ind2sub(size(guidata1.curr_drawing),ind);
        guidata1.curr_drawing(:) = 0;
        for ix1 = 1:length(ix)
            guidata1.drawing_idx(ix(ix1),iy(ix1),slice_first:slice_last) = guidata1.drawing_draw_color;
            guidata1.curr_drawing(ix(ix1),iy(ix1),slice_first:slice_last) = 1;
        end
        % Save Coordinates of Last Drawing
        idx_out = find(guidata1.curr_drawing(:)==1);
        [ix,iy,iz] = ind2sub(size(guidata1.curr_drawing),idx_out);
        guidata1.last_draw = [ix,iy,iz];
        guidata(handles.figure,guidata1);
        calc_full_ind
        update_image
    end
end

function clear_drawing_callback(~,~,~)
    answer1 = questdlg('Are you sure you want to clear the drawing?','Confirmation');
    if strcmp(answer1,'Yes')
        guidata1 = guidata(handles.figure);
        guidata1.drawing_idx = [];
        guidata1.drawing_overlay_alpha_data = []; 
        guidata1.draw_fname = [];
        if ishandle(handles.drawing_img)
            delete(handles.drawing_img)
        end
        guidata1.drawing_active = false;
        if strcmp(guidata1.currently_editing,'drawing')
            guidata1.currently_editing = 'background';
        end
        guidata(handles.figure,guidata1);
    end
end

%% Colormap, Colormap Indices, and Update_Image Functions:

% Update Colormap
function update_colormap
    guidata1 = guidata(handles.figure);
    guidata1.main_colormap = eval([guidata1.main_colormap_selection,'(',guidata1.m_string,')']);
    colormap(handles.background_axes,guidata1.main_colormap);
    if guidata1.overlay1_on
        guidata1.overlay1_colormap = eval([guidata1.overlay1_colormap_selection,'(',guidata1.m_string,')']);
        colormap(handles.overlay1_axes,guidata1.overlay1_colormap);
    end
    if guidata1.overlay2_on
        guidata1.overlay2_colormap = eval([guidata1.overlay2_colormap_selection,'(',guidata1.m_string,')']);
        colormap(handles.overlay2_axes,guidata1.overlay2_colormap);
    end
    guidata(handles.figure,guidata1);
end

% Calculate Colormap Indices
function calc_full_ind(~)
    guidata1 = guidata(handles.figure);
    % Find ColorMap Indices from Main Image:
    if ~guidata1.autoscale
        guidata1.imgfull_idx = min(guidata1.m,round((guidata1.m-1).*(guidata1.imgfull-guidata1.cmin)./(guidata1.cmax-guidata1.cmin))+1);
        guidata1.imgfull_idx(guidata1.imgfull_idx<=0) = 1;
        guidata1.main_colorbar_lim = [min(guidata1.imgfull_idx(:)),max(guidata1.imgfull_idx(:))+1];
        if guidata1.main_colorbar_lim(1)>=guidata1.main_colorbar_lim(2)
            guidata1.main_colorbar_lim(2) = guidata1.main_colorbar_lim(2)+1;
        end
    end
    % Edit Alpha Data for Drawing Overlay:
    if guidata1.drawing_active || nargin==1
        guidata1.drawing_overlay_alpha_data = zeros(size(guidata1.drawing_idx)); 
        guidata1.drawing_overlay_alpha_data(guidata1.drawing_idx(:)>1) = guidata1.overlay_alpha;
    end
    % Find ColorMap Indices for Overlay 1:
    if guidata1.overlay1_on
        guidata1.overlay1_idx = min(guidata1.m,round((guidata1.m-1).*(guidata1.overlay1_imgfull-guidata1.cmin_overlay1)./(guidata1.cmax_overlay1-guidata1.cmin_overlay1))+1); % add m to these
        guidata1.overlay1_idx(guidata1.overlay1_idx<=0) = 1;
        guidata1.overlay1_colorbar_lim = [min(guidata1.overlay1_idx(:)),max(guidata1.overlay1_idx(:))+1];
        guidata1.overlay1_alpha_data = zeros(size(guidata1.overlay1_idx)); 
        guidata1.overlay1_alpha_data(guidata1.overlay1_imgfull(:)~=0) = guidata1.overlay1_alpha;
    end
    % Find Colormap Indices for Overlay 2:
    if guidata1.overlay2_on
        guidata1.overlay2_idx = min(guidata1.m,round((guidata1.m-1).*(guidata1.overlay2_imgfull-guidata1.cmin_overlay2)./(guidata1.cmax_overlay2-guidata1.cmin_overlay2))+1); % add m to these
        guidata1.overlay2_idx(guidata1.overlay2_idx<=0) = 1;
        guidata1.overlay2_colorbar_lim = [min(guidata1.overlay2_idx(:)),max(guidata1.overlay2_idx(:))+1];        
        guidata1.overlay2_alpha_data = zeros(size(guidata1.overlay2_idx)); 
        guidata1.overlay2_alpha_data(guidata1.overlay2_imgfull(:)~=0) = guidata1.overlay2_alpha;
    end
    guidata(handles.figure,guidata1);
end

function update_image(~,~)
    guidata1 = guidata(handles.figure);
    % Set main figure as current if overlay window is open:
    if ishandle(guidata1.overlay_figure) || nargin==2
        set(0, 'CurrentFigure', handles.figure)
    end
    % Update Main Image (Background):
    if isgraphics(handles.background_axes,'axes')
        guidata1.xticks = get(handles.background_axes,'XTick');
        guidata1.yticks = get(handles.background_axes,'YTick');
        guidata1.ax_xlim = [1,numel(guidata1.yind)];
        guidata1.ax_ylim = [1,numel(guidata1.xind)];
    end
    guidata1.refresh_img = false;
    if guidata1.autoscale 
        guidata1.update_colorbar = 1;
        background_slice = guidata1.imgfull(guidata1.xind,guidata1.yind,guidata1.curr_slice);
        guidata1.cmin = min(background_slice(:));  % Minimum color value
        guidata1.cmax = max(background_slice(:));  % Maximum color value
        if guidata1.cmin==guidata1.cmax; guidata1.cmax = guidata1.cmin + 1; end
        idx1 = min(guidata1.m,round((guidata1.m-1)*(background_slice-guidata1.cmin)/(guidata1.cmax-guidata1.cmin))+1);
        guidata1.main_colorbar_lim = [min(idx1(:)),max(idx1(:))+1];
        if guidata1.main_colorbar_lim(1)>=guidata1.main_colorbar_lim(2)
            guidata1.main_colorbar_lim(2) = guidata1.main_colorbar_lim(2)+1;
        end
%         if ~guidata1.refresh_img && nargin==0
        if nargin==0
            handles.background_img.CData = idx1;
        else
%             if isgraphics(handles.background_img); delete(handles.background_img); end
            if nargin==2 && isempty(guidata1.parsed_inputs.axes)
%                 if isgraphics(handles.background_axes); delete(handles.background_axes); end
                handles.background_img = image(idx1);
                handles.background_axes = gca;
                set(handles.background_axes,'Color',guidata1.axes_color,...
                    'XColor',guidata1.axes_color,'YColor',...
                    guidata1.axes_color,'ZColor',guidata1.axes_color,...
                    'GridColor',guidata1.axes_color)
            else 
                handles.background_img = image('CData',idx1,...
                    'Parent',handles.background_axes);
                if ~isempty(guidata1.parsed_inputs.axes) && nargin==2
                    set(handles.background_axes,'YDir','reverse');
                    handles.background_img.Parent = guidata1.parsed_inputs.axes;
                end
            end
            guidata1.xticks = get(handles.background_axes,'XTick');
            guidata1.yticks = get(handles.background_axes,'YTick');
            guidata1.ax_xlim = [1,numel(guidata1.yind)];
            guidata1.ax_ylim = [1,numel(guidata1.xind)];
            hold 'on';
        end
    else
        background_slice = guidata1.imgfull(guidata1.xind,guidata1.yind,guidata1.curr_slice);
        idx1 = min(guidata1.m,round((guidata1.m-1)*(background_slice-guidata1.cmin)/(guidata1.cmax-guidata1.cmin))+1);
        guidata1.main_colorbar_lim = [min(idx1(:)),max(idx1(:))+1];
        if guidata1.main_colorbar_lim(1)>=guidata1.main_colorbar_lim(2)
            guidata1.main_colorbar_lim(2) = guidata1.main_colorbar_lim(2)+1;
        end
%         if ~guidata1.refresh_img
            set(handles.background_img,'CData',guidata1.imgfull_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice));
%         else
%             hold off; delete(handles.background_img)
%             if nargin==2 && isempty(guidata1.parsed_inputs.axes)
%                 if isgraphics(handles.background_axes); delete(handles.background_axes); end
%                 handles.background_img = image('CData',...
%                     guidata1.imgfull_idx(guidata1.xind,guidata1.yind,...
%                     guidata1.curr_slice),'Parent',handles.background_axes);
%                 handles.background_axes = gca;
%                 hold(handles.background_axes,'on')
%             else 
%                 handles.background_img = image('Parent',handles.background_axes,...
%                     'CData',guidata1.imgfull_idx(guidata1.xind,guidata1.yind,...
%                     guidata1.curr_slice));
%             end
%             hold(handles.background_axes,'on');
%         end
    end
    % Update Drawing:
    if guidata1.drawing_active || nargin==1
        if ishandle(handles.drawing_img)
            if ~guidata1.refresh_img
                set(handles.drawing_img,'Parent',handles.drawing_axes,...
                    'CData',guidata1.drawing_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
                    'AlphaDataMapping','none','AlphaData',guidata1.drawing_overlay_alpha_data(guidata1.xind,guidata1.yind,guidata1.curr_slice));
%             else
%                 delete(handles.drawing_img)
%                 handles.drawing_img = image('CData',guidata1.drawing_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
%                     'AlphaDataMapping','none','AlphaData',guidata1.drawing_overlay_alpha_data(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
%                     'Parent',handles.drawing_axes); % parent axes should be already present);
            end      
        elseif ~isempty(guidata1.drawing_idx)
            if ~isgraphics(handles.drawing_axes,'axes')
                calc_full_ind
                handles.drawing_axes = axes('Visible','off'); % must remain invisible
                handles.drawing_axes.CLim = [1,8];
                handles.drawing_axes.YDir = 'reverse';
                colormap(handles.drawing_axes,guidata1.drawing_colormap);
            end
            handles.drawing_img = image('Parent',handles.drawing_axes,'AlphaDataMapping','none',...
                'CData',guidata1.drawing_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
                'AlphaData',guidata1.drawing_overlay_alpha_data(guidata1.xind,guidata1.yind,guidata1.curr_slice)); 
            guidata1.refresh_img = 1;
        end
    end
    % Update Overlay 1:
    if guidata1.overlay1_on
        if isgraphics(handles.overlay1_img,'image')
            if ~guidata1.refresh_img
                set(handles.overlay1_img,'CData',guidata1.overlay1_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
                    'AlphaData',guidata1.overlay1_alpha_data(guidata1.xind,guidata1.yind,guidata1.curr_slice));
%             else
%                 delete(handles.overlay1_img)
%                 handles.overlay1_img = image(guidata1.overlay1_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
%                     'AlphaData',guidata1.overlay1_alpha_data(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
%                     'Parent',handles.overlay1_axes); % parent axes should be already present);
            end
        else
            if ~isgraphics(handles.overlay1_axes,'axes')
                handles.overlay1_axes = axes('Visible','off'); % must remain invisible
                handles.overlay1_axes.YDir = 'reverse';
                if ~isempty(guidata1.parsed_inputs.axes)
                    set(handles.overlay1_axes,'Position',guidata1.ax_pos_input)
                end
            end
            colormap(handles.overlay1_axes,guidata1.overlay1_colormap);
            handles.overlay1_img = image('Parent',handles.overlay1_axes,'CData',guidata1.overlay1_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
                'AlphaDataMapping','none','AlphaData',guidata1.overlay1_alpha_data(guidata1.xind,guidata1.yind,guidata1.curr_slice)); 
            guidata1.refresh_img = 1;
        end
    end
    % Update Overlay 2:
    if guidata1.overlay2_on
        if isgraphics(handles.overlay2_img,'image')
            if ~guidata1.refresh_img
                set(handles.overlay2_img,'CData',guidata1.overlay2_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
                    'AlphaData',guidata1.overlay2_alpha_data(guidata1.xind,guidata1.yind,guidata1.curr_slice));
%             else
%                 delete(handles.overlay2_img)
%                 handles.overlay2_img = image('Parent',handles.overlay2_axes,'CData',guidata1.overlay2_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
%                     'AlphaData',guidata1.overlay2_alpha_data(guidata1.xind,guidata1.yind,guidata1.curr_slice)); % parent axes should be already present
            end
        else
            if ~isgraphics(handles.overlay2_axes,'axes')
                handles.overlay2_axes = axes('Visible','off'); % must remain invisible
                handles.overlay2_axes.YDir = 'reverse';
                if ~isempty(guidata1.parsed_inputs.axes)
                    set(handles.overlay2_axes,'Position',guidata1.ax_pos_input)
                end
            end
            colormap(handles.overlay2_axes,guidata1.overlay2_colormap);
            handles.overlay2_img = image('Parent',handles.overlay2_axes,'CData',guidata1.overlay2_idx(guidata1.xind,guidata1.yind,guidata1.curr_slice),...
                'AlphaDataMapping','none','AlphaData',guidata1.overlay2_alpha_data(guidata1.xind,guidata1.yind,guidata1.curr_slice)); 
            guidata1.refresh_img = 1;
        end
    end
    % Add/Remove Colorbar(s)
    if guidata1.colorbar && guidata1.update_colorbar % why is this char?
        guidata1.refresh_img = 1;
        if guidata1.autoscale || guidata1.refresh_img
            guidata1.refresh_img = 1;
            % Determine which axes to associate colorbar with   
            if guidata1.main_colorbar_on
                if isgraphics(guidata1.h_colorbar_main,'colorbar')
                    set(guidata1.h_colorbar_main,'Visible','on','Color',guidata1.axes_color)
                else
                    guidata1.h_colorbar_main = colorbar('Peer',handles.background_axes,...
                        'Position',guidata1.colorbar_pos,'Color',guidata1.axes_color);
                end
                % Set Colorbar Ticks
                guidata1.h_colorbar_main.Limits = guidata1.main_colorbar_lim;
                guidata1.h_colorbar_main.LimitsMode = 'auto';
                guidata1.colorbar_main_cvec = linspace(guidata1.cmin,guidata1.cmax,guidata1.m);
                if (guidata1.cmax-guidata1.cmin)>(.15*guidata1.m)
                    guidata1.colorbar_main_cvec = round(guidata1.colorbar_main_cvec);
                    % TODO
%                     guidata1.h_colorbar_main.TickLabels = cellstr(sprintf('%1g\n',guidata1.colorbar_main_cvec(guidata1.h_colorbar_main.Ticks)));                    
                else
                    % TODO
%                     guidata1.h_colorbar_main.TickLabels = cellstr(sprintf('%4.2g\n',guidata1.colorbar_main_cvec(guidata1.h_colorbar_main.Ticks)));
                end
                if isgraphics(guidata1.h_colorbar_overlay1,'colorbar')
                    set(guidata1.h_colorbar_overlay1,'Visible','off');
                end
                if isgraphics(guidata1.h_colorbar_overlay2,'colorbar')
                    set(guidata1.h_colorbar_overlay2,'Visible','off');
                end
            elseif guidata1.overlay1_colorbar_on && guidata1.overlay1_on
                if isgraphics(guidata1.h_colorbar_main,'colorbar')
                    set(guidata1.h_colorbar_main,'Visible','off');
                end
                if isgraphics(guidata1.h_colorbar_overlay2,'colorbar')
                    set(guidata1.h_colorbar_overlay2,'Visible','off');
                end                
                if isgraphics(guidata1.h_colorbar_overlay1,'colorbar')
                    set(guidata1.h_colorbar_overlay1,'Visible','on','Color',guidata1.axes_color)
                else
                    guidata1.h_colorbar_overlay1 = colorbar('Peer',handles.overlay1_axes,...
                    'Position',guidata1.colorbar_pos,'Color',guidata1.axes_color);
                end
                % Set Colorbar Ticks
                guidata1.h_colorbar_overlay1.Limits = guidata1.overlay1_colorbar_lim;
                guidata1.h_colorbar_overlay1.LimitsMode = 'manual';
                guidata1.overlay1_cvec = linspace(guidata1.cmin_overlay1,guidata1.cmax_overlay1,guidata1.m);
                if (guidata1.cmax_overlay1 - guidata1.cmin_overlay1)>(.15*guidata1.m)
                    guidata1.overlay1_cvec = round(guidata1.overlay1_cvec);
                    % TODO
%                     guidata1.h_colorbar_overlay1.TickLabels = cellstr(sprintf('%1g\n',guidata1.overlay1_cvec(guidata1.h_colorbar_overlay1.Ticks)));
                else
                    % TODO
%                     guidata1.h_colorbar_overlay1.TickLabels = cellstr(sprintf('%4.2g\n',guidata1.overlay1_cvec(guidata1.h_colorbar_overlay1.Ticks)));
                end                
            elseif guidata1.overlay2_colorbar_on && guidata1.overlay2_on
                if isgraphics(guidata1.h_colorbar_main,'colorbar')
                    set(guidata1.h_colorbar_main,'Visible','off');
                end
                if isgraphics(guidata1.h_colorbar_overlay1,'colorbar')
                    set(guidata1.h_colorbar_overlay1,'Visible','off');
                end                
                if isgraphics(guidata1.h_colorbar_overlay2,'colorbar')
                    set(guidata1.h_colorbar_overlay2,'Visible','on','Color',guidata1.axes_color)
                else
                    guidata1.h_colorbar_overlay2 = colorbar(handles.overlay2_axes,...
                    'Position',guidata1.colorbar_pos,'Color',guidata1.axes_color);
                end
                % Set Colorbar Ticks
                guidata1.h_colorbar_overlay2.Limits = guidata1.overlay2_colorbar_lim;
                guidata1.h_colorbar_overlay2.LimitsMode = 'manual';
                guidata1.overlay2_cvec = linspace(guidata1.cmin_overlay2,guidata1.cmax_overlay2,guidata1.m);
                if (guidata1.cmax_overlay2 - guidata1.cmin_overlay2)>(.15*guidata1.m)
                    guidata1.overlay2_cvec = round(guidata1.overlay2_cvec);
                    guidata1.h_colorbar_overlay2.TickLabels = cellstr(sprintf('%1g\n',guidata1.overlay2_cvec(guidata1.h_colorbar_overlay2.Ticks)));
                else
                    guidata1.h_colorbar_overlay2.TickLabels = cellstr(sprintf('%4.2g\n',guidata1.overlay2_cvec(guidata1.h_colorbar_overlay2.Ticks)));
                end  
            end
        end
    else
        if isgraphics(guidata1.h_colorbar_main,'colorbar')
            set(guidata1.h_colorbar_main,'Visible','off');
        end
        if isgraphics(guidata1.h_colorbar_overlay1,'colorbar')
            set(guidata1.h_colorbar_overlay1,'Visible','off');
        end                
        if isgraphics(guidata1.h_colorbar_overlay2,'colorbar')
            set(guidata1.h_colorbar_overlay2,'Visible','off')
        end
    end
    % Update Title:
    if guidata1.title_on
        guidata1.refresh_img = 1;
        if ~isgraphics(guidata1.h_title)
            guidata1.h_title = title(['Slice ',num2str(guidata1.curr_slice)],...
                'FontName','Helvetica','FontSize',10,'FontWeight','Bold',...
                'Units','normalized','Visible','on','Color',guidata1.axes_color); 
            guidata1.title_pos = get(guidata1.h_title,'Position');
            guidata1.title_pos(2) = 1.003*guidata1.title_pos(2);
            set(guidata1.h_title,'Position',guidata1.title_pos)
        else
            guidata1.h_title.String = ['Slice ',num2str(guidata1.curr_slice)];
            guidata1.h_title.Visible = 'on';
        end
    elseif isgraphics(guidata1.h_title)
        guidata1.h_title.Visible = 'off';
    end
    % Save Data:
    guidata(handles.figure,guidata1);
    if guidata1.refresh_img
        if nargin==2 % initial call
            reposition_axes(true);
        else
            reposition_axes;
        end
    end
end

function reposition_axes(~)
    guidata1 = guidata(handles.figure);
    % Determine Positioning
    if guidata1.title_on && guidata1.colorbar && guidata1.axis_tick_on % 1
        guidata1.curr_axis_pos = guidata1.ax_pos_all_on;
        guidata1.prev_state = guidata1.curr_state; guidata1.curr_state = 1;
    elseif ~guidata1.title_on && ~guidata1.colorbar && ~guidata1.axis_tick_on % 2
        guidata1.curr_axis_pos = guidata1.ax_pos_all_off;
        guidata1.prev_state = guidata1.curr_state; guidata1.curr_state = 2;
    elseif guidata1.title_on && ~guidata1.colorbar && guidata1.axis_tick_on % 3
        guidata1.curr_axis_pos = guidata1.ax_pos_no_colorbar;
        guidata1.prev_state = guidata1.curr_state; guidata1.curr_state = 3;
    elseif guidata1.title_on && guidata1.colorbar && ~guidata1.axis_tick_on % 4
        guidata1.curr_axis_pos = guidata1.ax_pos_no_tick;
        guidata1.prev_state = guidata1.curr_state; guidata1.curr_state = 4;
    elseif ~guidata1.title_on && guidata1.colorbar && guidata1.axis_tick_on % 5
        guidata1.curr_axis_pos = guidata1.ax_pos_no_title;
        guidata1.prev_state = guidata1.curr_state; guidata1.curr_state = 5;
    elseif guidata1.title_on && ~guidata1.colorbar && ~guidata1.axis_tick_on % 6
        guidata1.curr_axis_pos = guidata1.ax_pos_title_only;
        guidata1.prev_state = guidata1.curr_state; guidata1.curr_state = 6;
    elseif ~guidata1.title_on && guidata1.colorbar && ~guidata1.axis_tick_on % 7
        guidata1.curr_axis_pos = guidata1.ax_pos_colorbar_only;
        guidata1.prev_state = guidata1.curr_state; guidata1.curr_state = 7;
    elseif ~guidata1.title_on && ~guidata1.colorbar && guidata1.axis_tick_on % 8
        guidata1.curr_axis_pos = guidata1.ax_pos_tick_only;
        guidata1.prev_state = guidata1.curr_state; guidata1.curr_state = 8;
    end
    if guidata1.prev_state~=guidata1.curr_state || nargin==1
        guidata1.x_ax_percent = guidata1.curr_axis_pos(3);
        guidata1.y_ax_percent = guidata1.curr_axis_pos(4);
        resize1 = true;
    else
        resize1 = false;
    end
    guidata1.colorbar_pos = [guidata1.curr_axis_pos(1)+guidata1.curr_axis_pos(3),...
        guidata1.curr_axis_pos(2),guidata1.colorbar_pos(3),guidata1.curr_axis_pos(4)];
    % Check if axes were specified by user (in which case, don't reposition)
    if ~isempty(parsed_inputs.axes) 
        guidata1.curr_axis_pos = guidata1.ax_pos_input; resize1 = false;
    end
    % Apply Positioning:
    set(handles.background_axes,'XLim',guidata1.ax_xlim,'YLim',guidata1.ax_ylim,...
        'XTick',guidata1.xticks,'YTick',guidata1.yticks,'Position',guidata1.curr_axis_pos)
    if guidata1.main_colorbar_on && ishandle(guidata1.h_colorbar_main)
        set(guidata1.h_colorbar_main,'Position',[guidata1.colorbar_pos(1),...
            guidata1.curr_axis_pos(2),guidata1.colorbar_pos(3),guidata1.curr_axis_pos(4)])
    end
    if guidata1.drawing_active
        set(handles.drawing_axes,'XLim',guidata1.ax_xlim,'YLim',guidata1.ax_ylim,...
            'XTick',guidata1.xticks,'YTick',guidata1.yticks,'Position',guidata1.curr_axis_pos)
    end
    if guidata1.overlay1_on
        set(handles.overlay1_axes,'XLim',guidata1.ax_xlim,'YLim',guidata1.ax_ylim,...
            'XTick',guidata1.xticks,'YTick',guidata1.yticks,'Position',guidata1.curr_axis_pos)
        if guidata1.overlay1_colorbar_on && ishandle(guidata1.h_colorbar_overlay1)
            set(guidata1.h_colorbar_overlay1,'Position',[guidata1.colorbar_pos(1),...
                guidata1.curr_axis_pos(2),guidata1.colorbar_pos(3),guidata1.curr_axis_pos(4)])
        end
    end
    if guidata1.overlay2_on
        set(handles.overlay2_axes,'XLim',guidata1.ax_xlim,'YLim',guidata1.ax_ylim,...
            'XTick',guidata1.xticks,'YTick',guidata1.yticks,'Position',guidata1.curr_axis_pos)
        if guidata1.overlay2_colorbar_on && ishandle(guidata1.h_colorbar_overlay2) 
            set(guidata1.h_colorbar_overlay2,'Position',[guidata1.colorbar_pos(1),...
                guidata1.curr_axis_pos(2),guidata1.colorbar_pos(3),guidata1.curr_axis_pos(4)])
        end
    end
    if resize1; resize_figure; end
    guidata(handles.figure,guidata1);
end
            
end % End Function