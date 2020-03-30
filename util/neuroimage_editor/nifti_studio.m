function [handles] = nifti_studio(varargin)
% NIfTI Studio:  
%   A GUI for navigating, visualizing, and editing 3D NIfTI images 
%   (file types: .nii, .nii.gz, .img/.hdr)
% 
% Author:
%   Elliot A. Layden, The University of Chicago, 2016-19
% 
% Usage: 
% To begin, simply type "nifti_studio" into the command line, adding 
% any desired name-value pair arguments. Alternatively, simply right-click 
% "nifti_studio.m" -> Run. Next, a file selection menu will be 
% displayed; select your desired 3D image file, and begin viewing/editing.
% 
% Optional Output:
%   'handles',          Handles structure for GUI figure, axes, etc.
% 
% Optional Inputs (Name-Value Pair Arguments):
%   (Note: these can also be loaded within the GUI)
%   'background',       filename or path-filename to crossa background image to
%                       be loaded; can also be loaded image structure 
%   'overlay',          filename or path-filename to an image to be loaded
%                       as an overlay; can also be loaded image structure
%   'colorbar',         customization: true/false
%   'title_on',         customization: true/false
%   'axis_tick_on',     customization: true/false
%   'colormap',         colormap for main figure (background), specified as
%                       string, e.g., 'colormap','jet'
%   'axes',             an axes graphics object (handle); used to embed a 
%                       NIfTI Studio window within an already present
%                       figure/axes
%   'apply_header',     1 or 0, denoting whether to apply affine matrix in
%                       header (Default: 1)
%   'background_caxis', 1x2 vector [cmin, cmax] to scale colors of
%                       background image
%   'overlay_caxis',    1x2 vector [cmin, cmax] to scale colors of
%                       overlay image
% 
% [Note: if you encounter error messages when loading a file that includes 
% "non-orthogonal shearing", this means that applying the affine to the 
% image is introducing distortion beyond the tolerance range of NIFTI Tools. 
% In this case, try loading the image without applying the affine/header 
% info (nifti_studio('apply_header',0); equivalent to using 
% load_untouch_nii.m in NIFTI Tools)]
%
% General Features:
% -use up and down arrow keys to navigate slices
% -use the mouse scroll wheel to zoom in or out
% -right click and drag OR Ctrl + left click:  pans left/right, up/down
% -click and drag the mouse: 
%   in Crosshair mode: to obtain coordinates and pixel intensity value at
%     the clicked location
%   in Drawing mode:  draw shapes; for instance, drawing an arc will 
%     cause the interior of the arc to become filled with color; any closed 
%     figure will also be filled with color (particularly useful for erasing 
%     undesirable parts of images like artifacts, or for creating ROIs)
%     -can draw to edit background image data, to draw new ROIs which can
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
% 's' hotkey saves whichever file is currently being edited
% 
% GUI Menus:
% FILE
%   -Open:  select new file to open as a background or overlay image
%   -Close Overlay...:  close an open overlay image
%   -Save:  save a background or overlay image along with any edits made 
%           using the same filename as was loaded
%   -Save As:  save a background or overlay image with a new file name
%   -Exit (hotkey: 'esc'):  exits NIfTI Studio, prompting user to verify
% SELECT
%   -select either the background image or a loaded overlay to edit
%   -New Overlay...:  create a new overlay image
% TOOLS
%   -Crosshair (hotkey: 'c'):  enables the crosshair tool, which outputs
%   voxel location and value when the image is clicked
%   -Draw (hotkey: 'd'):  enables the drawing tool, which allows edits to
%   the selected image (e.g., create new spherical ROIs)
%   -Pan (hotkey:  'p'):  enables the pan tool, which allows the user to
%   click and drag to move the image (adjust axes limits). Note that the
%   pan tool can also be accessed while drawing via ctrl + click or right
%   click
% EDIT
%   -Undo (hotkey: 'u'):  undo last action (drawing or orientation change)
%   -Redo (hotkey: 'r'):  redo actions following calls to undo
%   -Go to Slice... (hotkey: 'g'):  navigate to a specified slice
%   -Revert to Defaults:  change display settings back to defaults
% DISPLAY
%   -Orientation:   change orientation (Coronal, Sagittal, Axial, 
%                   3D Display, Mosaic of slices)
%   -Colorbar (On/Off):   turn on/off colorbar
%   -Axis Tick (On/Off):  turn on/off axis ticks
%   -Slice # (On/Off):    turn on/off title which displays slice #
% DRAW
%   -Select Draw Color: specify an image intensity to use for drawing 
%   -Shapes:  select a shape to draw (No Shape (Manual Trace), Circle,
%             Rectangle, Sphere)
%   -Propagate through Slices:  propagate the most recent drawing through
%                               specified slices
%   -Add Border:  specify a # of voxels to create a border within the 
%                 current slice
%   -Fill Slice:  fill an entire slice with current draw color
% COLORS
%   -Colormap:      change colormap of the selected image
%   -Color Scale: 	change the color scale (climit) for the selected image
% 
% Dependencies: 
% NIFTI_tools (Shen, 2005) must be on Matlab's path
%   Link: http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
%   Note that the necessary functions have been included in the 
%   NIfTI Studio download folder, so no further action is required.

% Note: colorbars and overlays will not function properly for Matlab
% versions prior to 2014b, due to the major graphics update which came in
% 2014b

% Author:  Elliot Layden, University of Chicago, 2016-2019
% Contact:  elliot.layden@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Image and Initialize to Display a Middle Slice:
figure_color = [.2,.2,.2];
axes_color = ones(1,3);

% Identify Function Path and Add Helper Scripts:
script_fullpath = mfilename('fullpath');
[script_path,~,~] = fileparts(script_fullpath);
addpath(genpath(script_path))

% Get Inputs:
inputs = varargin;
parsed_inputs = struct('background',[],'overlay',[],...
    'colorbar_on',[],'title_on',[],'axis_tick_on',[],...
    'colormap',[],'axes',[],'apply_header',1,...
    'background_caxis',[],'overlay_caxis',[]);
poss_input = {'background','overlay','colorbar_on','title_on',...
    'axis_tick_on','colormap','axes','apply_header',...
    'background_caxis','overlay_caxis'};

% column of internal cells == OR, row of internal cells == AND
input_types = {{'char','file';'struct',''},{'char','file';'struct',''},...
    {'logical';'numeric'},{'logical';'numeric'},{'logical';'numeric'},...
    {'char'},{'axes'},{'logical';'numeric'},...
    {'vector'},{'vector'}}; 
parsed_inputs = getInputs(inputs, parsed_inputs, input_types);

% If missing background input, clear others:
if isempty(parsed_inputs.background)
    parsed_inputs.overlay = [];
    parsed_inputs.colormap = [];
    parsed_inputs.axes = [];
end

% Determine whether to apply header:
apply_header = parsed_inputs.apply_header;
if apply_header
    untouch_nii = false;
elseif ~apply_header
    untouch_nii = true;
end

% Initialize Handles & Options:
handles = struct('figure',99.999,'axes',99.999,'images',99.999);
handles.axes = {99.999}; handles.images = {99.999};

% Initialize Figure (auto-detect screen size):
set(0,'units','pixels'); 
screen_res = get(0,'ScreenSize');
figure_pos = [.25*screen_res(3), .065*screen_res(4), ...
    .5*screen_res(3), .86*screen_res(4)];
if any(figure_pos<=0) % avoid ScreenSize errors
    figure_pos = [342,50,683,660]; 
    screen_res = [1,1,(figure_pos(3:4)./[.5,.86])];
end

if isempty(parsed_inputs.axes) % if no axes supplied
    handles.figure = figure('Position',figure_pos,'MenuBar','none',...
        'Name','NIfTI Studio','NumberTitle','off','Color',figure_color,...
        'Visible','off','doublebuffer','on','Interruptible','off');
else
    handles.figure = get(parsed_inputs.axes,'Parent');
end
figure(handles.figure)
 
% Intialize GUI Data:
selectedImage = 1;
nOverlays = 0; % global counter that never decrements (even when delete previous overlay, numbering continues based on this)
if ~isempty(parsed_inputs.axes)
    ax_pos_input = get(parsed_inputs.axes,'Position');
end
fig_height = figure_pos(4); 
screen_mid_x = .5*screen_res(3);
if isempty(parsed_inputs.axes)
    handles.axes{selectedImage} = 16.48382;
else
    handles.axes{selectedImage} = parsed_inputs.axes;
end

% Main Data Storage:
imageData = {0}; alphaData = {0};
cmin = {0}; cmax = {0};
img = []; 
fullPaths = {[]}; % stores full paths of background and overlay images

origin = [];
physical_units = [];
units = 'physical';
filename = []; fpath = []; aspect_ratio = []; fig_width = []; dimperm = [];
draw_on = false; pan_on = false;
window_name = []; dim = []; pixdim = []; voxSize = []; xwidth = []; yheight = [];
xdim = []; ydim = []; zdim = []; middle_slice = [];
curr_slice = []; curr_drawing = []; 
xind = []; yind = []; 
xmax = []; xmin = []; ymax = []; ymin = [];
ax_xlim = []; ax_ylim = [];
x_slice = []; y_slice = []; z_slice = []; 
in_motion = false;
draw_color = 1; 
erase_draw = false;
idx_draw = [];
scroll_count = 0; 

launch_mosaic_string = 'Mosaic (Beta)';
launch_3d_string = '3D Display';

% Initialize undo/redo:
undoNum = 0; undoLimit = 10;
undoCache = struct('selectedImage', [], 'action', [], ...
    'orientation', [], 'idx', [], 'color', [], 'alpha', []);
redoCache = struct('selectedImage', [], 'action', [], ...
    'orientation', [], 'idx', [], 'color', [], 'alpha', []);

poss_ind = []; save_points = []; 
n_ticks_x = 10; n_ticks_y = 10;
xticks = []; yticks = [];

% Menu Handles:
orientation_labels = {'Coronal','Sagittal','Axial',launch_3d_string,launch_mosaic_string}; % 'Mosaic'
menu_orientations = zeros(1,length(orientation_labels)); 
menu_slice_number = []; menu_colorbar = []; menu_axis_tick = [];

% Axes Positions:
ax_pos_all_on =        [.09,  .04,  .8,  .92]; % All On
ax_pos_all_off =       [0,     0,    1,    1]; % All Off
ax_pos_no_title =      [.09,  .04,  .8,  .95]; % No Title
ax_pos_no_colorbar =   [.09,  .04,  .88,  .92]; % No Colorbar
ax_pos_no_tick =       [0,     0,    .88,  .95]; % No Tick
ax_pos_title_only =    [0,     0,    1,    .95]; % Title only
ax_pos_colorbar_only = [0,     0,    .88,  1];% Colorbar only
ax_pos_tick_only =     [.09,  .04, .88, .95]; % Tick Only
colorbar_pos = [.927, ax_pos_all_on(2), 0.0390, ax_pos_all_on(4)];
curr_axis_pos = ax_pos_all_on; 
x_ax_percent = curr_axis_pos(3);
y_ax_percent = curr_axis_pos(4); 
prev_state = 1; curr_state = 1; % various combo's of axes objects

% Colormap & Overlay Settings:
colormap_n = 200; % Number of distinct colors in the current colormap
n_colorbar_ticks = 10; % number of ticks on colorbar

% Colormaps & Colorbars:
h_colorbar = 20.19347;
colorMapStr = {'gray','jet'};
colormap_opts = {'gray','jet','hot','cool','hsv','bone','colorcube','copper',...
    'spring','summer','winter','pink'};
colormaps{selectedImage} = colormap(eval([colorMapStr{selectedImage},'(',num2str(colormap_n),')']));
colorscaleType = {1, 2}; % (1) Slice min/max, (2) Global min/max, (3) Custom

% Opacity / alpha
alpha_opts = {'Opaque','90%','80%','70%','60%','50%','40%','30%',...
    '20%','10%','Invisible'};
alpha_values = 1:-.1:0;
alphaValue = {1, .6}; % overlay transparency (default: 60%)

% Drawing ROIs
curr_drawing = []; 
shape = 0; shape_ind = zeros(2,2); sphere_radius = 1;
shape_coords = zeros(64,2); 
prev_shape_colors = zeros(64,1); 
prev_alpha = [];

% Other:
customizable = {'colorbar_on','title_on','axis_tick_on',...
    'extensions','last_nav_dir','colorMapStr'};
first_dir = pwd; last_nav_dir = pwd; 
extensions = {'*.img';'*.nii';'*.nii.gz'};
slice_orientation = 3; % default = z-dim
h_title = 10.48487; num_voxels_border = 0; 
scroll_zoom_equiv = [1,.9,.8,.7,.6,.5,.4,.3,.2,.1];

% Set user interface callback functions:
set(handles.figure,'WindowKeyPressFcn',@keypress_callback);
% Don't turn on click function if embedded graphic in user spec axes:
if isempty(parsed_inputs.axes)
    set(handles.figure,'WindowButtonDownFcn',@cursor_click_callback);
end
set(handles.figure,'WindowScrollWheelFcn',@scroll_zoom_callback);
if isempty(parsed_inputs.axes)
    set(handles.figure,'CloseRequestFcn',@closereq_callback)
else
    set(handles.figure,'CloseRequestFcn',@closereq_no_dlg)
end

% Get Settings:
succeeded = getSettings; % extracts custom settings from .txt if available

if ~isempty(parsed_inputs.title_on)
    title_on = parsed_inputs.title_on; 
elseif ~succeeded
    title_on = true;
end
if ~isempty(parsed_inputs.axis_tick_on)
    axis_tick_on = parsed_inputs.axis_tick_on;
elseif ~succeeded
    axis_tick_on = true;
end
if ~isempty(parsed_inputs.colorbar_on)
    colorbar_on = parsed_inputs.colorbar_on;
elseif ~succeeded
    colorbar_on = true;
end

% Background: Get Input Filename & Load
if isempty(parsed_inputs.(poss_input{1}))
    status = openNewBackground;
    if ~status; return; end
elseif ischar(parsed_inputs.(poss_input{1}))
    [fpath,filename,ext] = fileparts(parsed_inputs.(poss_input{1}));
    filename = [filename,ext];
    sort_exts(fpath,filename);
    load_img('char');
elseif isstruct(parsed_inputs.(poss_input{1}))
    load_img('struct');
end

%% Create Menu Items:
if isempty(parsed_inputs.axes)
    
    % FILE
    file_menu = uimenu(handles.figure,'Label','File');
    
    % OPEN
    open_menu = uimenu(file_menu,'Label','Open');
    uimenu(open_menu,'Label','Background Image','Callback',{@openNewBackground,1});
    uimenu(open_menu,'Label','Overlay Image','Callback',@openNewOverlay);
    uimenu(open_menu,'Label','New Overlay','Callback',@createOverlay);
        
    % CLOSE
    close_overlays_menu = uimenu(file_menu,'Label','Close Overlay...');
    h_close = [];
    
    % SAVE  
    save_menu = uimenu(file_menu,'Label','Save');
    uimenu(save_menu,'Label','Save Background Image','Callback',{@save_callback,1});
    uimenu(save_menu,'Label','Save Current Overlay','Callback',{@save_callback,2});
    
    % SAVE AS
    save_as_menu = uimenu(file_menu,'Label','Save As');
    uimenu(save_as_menu,'Label','Save Background Image As','Callback',{@saveas_callback,1});
    uimenu(save_as_menu,'Label','Save Current Overlay As','Callback',{@saveas_callback,2});         
    uimenu(file_menu,'Label','Exit                                          ''esc''','Callback',@closereq_callback);
        
    % SELECT
    select_menu = uimenu(handles.figure,'Label','Select');
    h_image(1) = uimenu(select_menu,'Label','Background Image',...
        'Checked','on','Callback',{@changeSelection,1});
    h_new_image = uimenu(select_menu,'Label','New Overlay...',...
        'Callback',@createOverlay);

    % TOOLS
    tools_menu = uimenu(handles.figure,'Label','Tools');
    tool_crosshair = uimenu(tools_menu,'Label','Crosshair                         ''c''',...
        'Callback',@crosshair_callback, 'Checked','on');
    tool_draw = uimenu(tools_menu,'Label','Draw                                ''d''',...
        'Callback',@draw_callback);
%     tool_zoom = uimenu(tools_menu,'Label','Zoom                               ''z''',...
%         'Callback',@zoom_callback);
    tool_pan = uimenu(tools_menu,'Label','Pan                                   ''p''',...
        'Callback',@pan_callback);

    % EDIT    
    edit_menu = uimenu(handles.figure,'Label','Edit');
    uimenu(edit_menu,'Label','Undo                         ''u''',...
        'Callback',@undoCallback);
    uimenu(edit_menu,'Label','Redo                          ''r''',...
        'Callback',@redoCallback);
    uimenu(edit_menu,'Label','Go to Slice...             ''g''',...
        'Callback',@goToSlice);
    uimenu(edit_menu,'Label','Go to Origin...           ''o''',...
        'Callback',@goToOrigin);
    uimenu(edit_menu,'Label','Revert to Defaults','Callback',@revert_defaults);    
end

% DISPLAY
display_opts_menu = uimenu(handles.figure,'Label','Display');
orient_menu = uimenu(display_opts_menu,'Label','Orientation');
for ixx = 1:length(orientation_labels)
    menu_orientations(ixx) = uimenu(orient_menu,'Label',...
        orientation_labels{ixx},'Checked','off',...
        'Callback',@reorient_callback); 
end
set(menu_orientations(3),'Checked','on');

% Units:
menu_units = uimenu(display_opts_menu,'Label','Units'); 
h_units(1) = uimenu(menu_units,'Label','Physical','Checked','on',...
    'Callback', @changeUnits); 
h_units(2) = uimenu(menu_units,'Label','Voxel','Checked','off',...
    'Callback', @changeUnits); 

if isempty(parsed_inputs.axes)
    if colorbar_on
        menu_colorbar = uimenu(display_opts_menu,'Label','Colorbar',...
            'Checked','on','Callback',@colorbar_callback);
    else
        menu_colorbar = uimenu(display_opts_menu,'Label','Colorbar',...
            'Checked','off','Callback',@colorbar_callback);
    end
    if axis_tick_on
        menu_axis_tick = uimenu(display_opts_menu,'Label','Axis Tick',...
            'Checked','on','Callback',@axis_tick_callback);
    else
        menu_axis_tick = uimenu(display_opts_menu,'Label','Axis Tick',...
            'Checked','off','Callback',@axis_tick_callback);
    end
    if title_on
        menu_slice_number = uimenu(display_opts_menu,'Label','Slice #',...
            'Checked','on','Callback',@title_toggle_callback);
    else
        menu_slice_number = uimenu(display_opts_menu,'Label','Slice #',...
            'Checked','off','Callback',@title_toggle_callback);
    end
    
    % DRAW    
    draw_menu = uimenu(handles.figure,'Label','Draw');
    h_color_menu = uimenu(draw_menu,'Label','Select Draw Color');
    h_color = zeros(1,11);
    h_color(1) = uimenu(h_color_menu,'Label','Erase',...
        'Callback',{@change_color_callback, 0});
    for j = 2:10
        h_color(j) = uimenu(h_color_menu,'Label',num2str(j-1),...
            'Callback',{@change_color_callback, j-1});
    end
    set(h_color(2), 'Checked', 'on')
    h_color(11) = uimenu(h_color_menu,'Label','Custom',...
        'Callback',{@change_color_callback, 10});
    h_shapes = zeros(1,4);
    h_shapes_menu = uimenu(draw_menu,'Label','Shapes');
    h_shapes(1) = uimenu(h_shapes_menu,'Label','No Shape (Manual Trace)',...
        'Checked','on','Callback',{@draw_shapes_callback, 0});
    h_shapes(2) = uimenu(h_shapes_menu,'Label','Circle','Callback',{@draw_shapes_callback,1});
    h_shapes(3) = uimenu(h_shapes_menu,'Label','Rectangle','Callback',{@draw_shapes_callback,2});
    h_shapes(4) = uimenu(h_shapes_menu,'Label','Sphere','Callback',{@draw_shapes_callback,3});
    uimenu(draw_menu,'Label','Propagate last draw through slices','Callback',@propagate_draw_callback);    
    uimenu(draw_menu,'Label','Add Border','Callback',@add_border_callback);
    uimenu(draw_menu,'Label','Fill Slice','Callback',@apply2whole_slice_callback);
    
    % COLOR
    color_menu = uimenu(handles.figure, 'Label', 'Colors');
    
    % Colormap / cmap
    colormap_menu = uimenu(color_menu,'Label','Colormap'); 
    h_colormap_menu = zeros(1, length(colormap_opts));
    for j = 1:length(colormap_opts)
        h_colormap_menu(j) = uimenu(colormap_menu, ...
            'Label', colormap_opts{j}, ...
            'Callback', {@colormap_callback, j});
    end
        
    % Colorscale / caxis
    h_colorscale_menu = zeros(1,3);   
    colorscale_menu = uimenu(color_menu,'Label','Color Scale');
    h_colorscale_menu(1) = uimenu(colorscale_menu,...
        'Label','Slice Min/Max','Callback',{@colorscale_callback, 1},...
        'Checked','on');
    h_colorscale_menu(2) = uimenu(colorscale_menu,...
        'Label','Global Min/Max','Callback',{@colorscale_callback, 2});
    h_colorscale_menu(3) = uimenu(colorscale_menu,...
        'Label','Custom...','Callback',{@colorscale_callback, 3});
        
    % Opacity / alpha
    opacity_menu = uimenu(color_menu,'Label','Opacity');
    h_opacity_menu = zeros(1, length(alpha_opts));
    for j = 1:length(alpha_opts)
        h_opacity_menu(j) = uimenu(opacity_menu,...
            'Label', alpha_opts{j},'Callback', {@opacity_callback, j});
    end
%     if overlay_present
%         set(h_opacity_menu(strcmp(alpha_opts, alphaValue{2})), 'Checked','on')
%     else
        set(h_opacity_menu(1), 'Checked','on')
        set(opacity_menu, 'Visible','off')
%     end
end

%% Customizations

% Must run to initialize axes etc., prior to the subsequent lines
updateImage

% Add Overlays If Requested:
if ~isempty(parsed_inputs.overlay)
    switch class(parsed_inputs.overlay)
        case 'char'
            openNewOverlay([], [], 'char');
        case 'struct'
            openNewOverlay([], [], 'struct');
        otherwise
            if isnumeric(parsed_inputs.overlay)
                openNewOverlay([], [], 'matrix');
            end
    end
end

% Adjust colormap of background image based on input
selectedImageHolder = selectedImage;
selectedImage = 1;
if ~isempty(parsed_inputs.colormap)
    if ismember(parsed_inputs.colormap, colormap_opts) 
        colormap_callback([], [], find(strcmp(parsed_inputs.colormap, colormap_opts)))
    else
        warning('Invalid ''colormap'' input.')
    end
else
    colormap_callback([], [], find(strcmp(colorMapStr{selectedImage}, colormap_opts)))
end
selectedImage = selectedImageHolder;
set(handles.figure,'Visible','on'); % Make Visible

% Adjust CAXES if specified:
if ~isempty(parsed_inputs.background_caxis)
    if parsed_inputs.background_caxis(2)~=parsed_inputs.background_caxis(1) && parsed_inputs.background_caxis(2)>parsed_inputs.background_caxis(1)
        cmin{1} = parsed_inputs.background_caxis(1);
        cmax{1} = parsed_inputs.background_caxis(2);
        colorscaleType{1} = 3;
        % Adjust colorscale menu checkmarks:
        set(h_colorscale_menu(2), 'Checked','off') % default
        set(h_colorscale_menu(3), 'Checked','on') % Custom
    else
        warning('Invalid background color axis specification. Ignoring inputs.')
        return;
    end
    refresh_img = 0; 
    updateImage;
end

if ~isempty(parsed_inputs.overlay_caxis)
    if parsed_inputs.overlay_caxis(2)~=parsed_inputs.overlay_caxis(1) && parsed_inputs.overlay_caxis(2)>parsed_inputs.overlay_caxis(1)
        cmin{2} = parsed_inputs.overlay_caxis(1);
        cmax{2} = parsed_inputs.overlay_caxis(2);
        colorscaleType{2} = 3;
        % Adjust colorscale menu checkmarks:
        set(h_colorscale_menu(2), 'Checked','off') % default
        set(h_colorscale_menu(3), 'Checked','on') % Custom
    else
        warning('Invalid background color axis specification. Ignoring inputs.')
        return;
    end
    refresh_img = 0; 
    updateImage;
end

repositionAxes(1);

% Turn on default tool 'crosshair':
crosshair_callback;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Background Image:
function load_img(img_type)
    switch img_type
        case 'struct'
            try
                img = parsed_inputs.(poss_input{1});
                try
                    fpath = img.fileprefix;
                    [~,filename,~] = fileparts(fpath);
                catch
                    filename = 'Image Structure';
                end
            catch
                error('Error loading image structure. Check fields and retry or load through GUI.')
            end
        case 'char'
            fullPaths{1} = fullfile(fpath, filename); 
            if untouch_nii
                img = load_untouch_nii(fullPaths{1});
            else
                try
                    img = load_nii(fullPaths{1});
                catch
                    img = load_untouch_nii(fullPaths{1});
                    untouch_nii = true;
                    warning('Non-orthogonal shearing detected in affine matrix. Image loaded without applying affine.')
                end
            end
    end
    
    % Determine voxel units & origin:
    try
        origin = img.hdr.hist.originator(1:3);
        switch bitand(img.hdr.dime.xyzt_units, 7) % see xform_nii.m, extra_nii_hdr.m
            case 1, physical_units = 'm';
            case 2, physical_units = 'mm';
            case 3, physical_units = 'microns';
            otherwise, physical_units = '';
        end
    catch
        physical_units = '';
        origin = img.hdr.dime.dim(2:4)/2;
        warning('Failed to retrieve voxel units.')
    end
    if any(origin==0)
       origin = img.hdr.dime.dim(2:4)/2;
    end
            
    % Change selected image back to background:
    if selectedImage ~= 1
        changeSelection([], [], 1)
    end
    
    % Close any open overlays:
    selectedImage = 1;
    nOverlays = 0; % reset
    if length(imageData)>1
       for i = 2:length(imageData)
          closeOverlay([], [], i)
       end
    end
    
    % Clear overlay menu items if present
    if exist('h_image', 'var') && numel(h_image) > 1
        h_image = h_image(1);
    end
    
    if exist('h_close', 'var') && numel(h_close) > 1
        h_close = h_close(1);
    end
        
    % Parse image input
    window_name = ['NIfTI Studio:    ',filename];
    set(handles.figure,'Name',window_name);  
    dim = img.hdr.dime.dim; 
    pixdim = img.hdr.dime.pixdim; 
    voxSize = pixdim(2:4);
    xwidth = dim(2)*pixdim(2); 
    yheight = dim(3)*pixdim(3); 
    imageData = {single(img.img)}; % conversion may enhance refresh speed
    imageData{1} = permute(imageData{1},[2,1,3]); % Match SPM View
    voxSize = voxSize([2,1,3]); % d=0
    origin = origin([2,1,3]);
    [xdim, ydim, zdim] = size(imageData{1});
    middle_slice = round(zdim/2); 
    curr_slice = middle_slice;
    xind = xdim:-1:1; 
    yind = 1:ydim;
    slice_orientation = 3; % default: Z-dim
    cmin = {min(imageData{1}(:))}; 
    cmax = {max(imageData{1}(:))}; 
    if cmin{1}==cmax{1}; cmax{1} = cmin{1} + 1; end
    xmax = max(xind); xmin = min(xind); 
    ymax = max(yind); ymin = min(yind);    
    ax_xlim = [1,numel(yind)]; ax_ylim = [1,numel(xind)]; 
    
    % Clear previous axes & images:
    for i = 1:numel(handles.axes)
        if isgraphics(handles.axes{1},'axes'); cla(handles.axes{1}); end
        if ishandle(handles.images{1}); delete(handles.images{1}); end
    end

    % Reset
    scroll_count = 0;
    
    % Reset undo's:
    undoNum = 0; 
    undoCache = struct('selectedImage', [], 'action', [], ...
        'orientation', [], 'idx', [], 'color', [], 'alpha', []);
    redoCache = [];
    
    % Determine possible combinations of x,y indices:
    poss_ind = zeros(xdim*ydim,2); count = 0;
    for ix = 1:xdim
        for jx = 1:ydim
            count = count + 1;
            poss_ind(count,:) = [ix,jx];
        end
    end
    
end

% Load Overlay Image:
function openNewOverlay(~, ~, load_type)
    
    if nargin<3
        load_type = 'char';
    end
    
    fullpath = [];
    switch load_type
        case 'char'
            if nargin < 3
                % Get input filename:
                cd(last_nav_dir);
                [overlayName, overlay_path] = uigetfile(extensions,...
                    'Select Overlay Image:','MultiSelect','off');
                cd(first_dir);
            elseif nargin==3
                [overlay_path, overlayName, ext] = fileparts(parsed_inputs.(poss_input{2}));
                overlayName = [overlayName, ext];  
            end
            if ischar(overlayName)
                fullpath = fullfile(overlay_path, overlayName);
                % Remember Extension Chosen:
                if ~isempty(overlay_path)
                    sort_exts(overlay_path,overlayName)
                end                    
                if untouch_nii
                    overlay1_img = load_untouch_nii(fullpath);
                else
                    overlay1_img = load_nii(fullpath);
                end
            else
                return
            end
        case 'struct'
            overlayName = 'Overlay 1';
            try
                overlay1_img = parsed_inputs.(poss_input{2});
            catch
                error('Error loading overlay image. Check fields and retry or load through GUI.')
            end
        case 'matrix'
            overlayName = 'Overlay 1';
            overlay1_img = img;
            overlay1_img.img = parsed_inputs.(poss_input{2});
    end
    
    % Before update, go to default orientation:
    save_orientation = [];
    if slice_orientation==1
            save_orientation = 1;
            reorient_callback(menu_orientations(3)) % revert to original orientation
    elseif slice_orientation==2
            save_orientation = 2;
            reorient_callback(menu_orientations(3)) % revert to original orientation
    end     
    
    % Extract & store image
    if ~all(size(overlay1_img.img)==dim(2:4))
        error('Error: Overlay dimensions do not match original image.');
    end
    imageData{end + 1} = single(overlay1_img.img); % conversion may enhance refresh speed
    selectedImage = numel(imageData);
    nOverlays = nOverlays + 1;
    fullPaths{selectedImage} = fullpath; % add filepath for image saving
    
    % Permute dimensions to match:
    imageData{selectedImage} = permute(imageData{selectedImage},[2,1,3]);
    alphaValue{selectedImage} = .6; 
    alphaData{selectedImage} = single(zeros(size(imageData{selectedImage})));
    use_idx = ((imageData{selectedImage}(:)~=0) + (~isnan(imageData{selectedImage}(:))))==2;
    alphaData{selectedImage}(use_idx) = alphaValue{selectedImage};
    
    % Add new menu item:
    incrementOverlayMenus(overlayName)
    
    % Calculate Image Color Indices
    colorscaleType{selectedImage} = 2; % Global min/max (default for overlays)
    cmin{selectedImage} = min(imageData{selectedImage}(:));
    cmax{selectedImage} = max(imageData{selectedImage}(:));
    if cmin{selectedImage} == cmax{selectedImage}
        cmax{selectedImage} = cmax{selectedImage} + 1;
    end
    
    % Create new axes:
    handles.axes{selectedImage} = axes('parent', handles.figure, 'Visible','off','YDir','reverse'); % must remain invisible
    handles.axes{selectedImage}.CLim = [cmin{selectedImage}, cmax{selectedImage}];

    % Add new colormap:
    colormap_callback([], [], 2)
    
    % Handle several uimenu and colobar changes:
    changeSelection([], [], selectedImage)
    
    % Return to Original Orientation if Changed:
    if ~isempty(save_orientation)
        reorient_callback(menu_orientations(save_orientation)) % have to have two inputs, see line 1047
    end
    
    % Update Figure:
    updateImage 
    updateColormap
    
end

function createOverlay(~,~,~)
    % Create new entries in imageData, alphaData, etc.
    selectedImage = numel(imageData) + 1;
    nOverlays = nOverlays + 1;
    
    % Image & Opacity data
    imageData{selectedImage} = single(nan(size(imageData{1}))); 
    alphaData{selectedImage} = single(zeros(size(imageData{1})));
    alphaValue{selectedImage} = .6; 
    fullPaths{selectedImage} = [];
    
    % Color mapping
    colorscaleType{selectedImage} = 2; % Global min/max (default for overlays)
    cmin{selectedImage} = 0; 
    cmax{selectedImage} = 1;
    
    % Create new axes:
    handles.axes{selectedImage} = axes('parent', handles.figure, 'Visible','off','YDir','reverse'); % must remain invisible
    handles.axes{selectedImage}.CLim = [cmin{selectedImage}, cmax{selectedImage}];
    
    % Add new colormap:
    colormap_callback([], [], 2)
    
    % Add additional uimenu's
    incrementOverlayMenus
    
    % Handle several uimenu and colobar changes:
    changeSelection([], [], selectedImage)
    
    updateImage
end

function incrementOverlayMenus(overlayName)
    % Remove New Overlay menu item, then regenerate it at the end of list:
    if exist('h_new_image','var') && ishandle(h_new_image)
        delete(h_new_image)
    end
    
    % Add new SELECT menu:
    if nargin==0 || isempty(overlayName)
        h_image(end + 1) = uimenu(select_menu,'Label',...
            sprintf('Overlay %1g', nOverlays),'Callback',...
            {@changeSelection, nOverlays + 1});
    else
        h_image(end + 1) = uimenu(select_menu,'Label',...
            overlayName,'Callback',...
            {@changeSelection, nOverlays + 1});
    end
    for i = 1:length(h_image)
        if ishandle(h_image(i))
            set(h_image(i),'Checked','off') % Uncheck all
        end
    end
    set(h_image(end),'Checked','on') % Check new
    
    % Re-make New Overlay button in SELECT menu:
    h_new_image = uimenu(select_menu,'Label','New Overlay...',...
        'Callback',@createOverlay);
    
    % Add new CLOSE overlay menu:
    if nargin==0 || isempty(overlayName)
        if ~isempty(h_close)
            h_close(end + 1) = uimenu(close_overlays_menu,'Label',...
                sprintf('Overlay %g', nOverlays),'Callback',{@closeOverlay, nOverlays + 1});
        else
            h_close(2) = uimenu(close_overlays_menu,'Label',...
                sprintf('Overlay %g', 1),'Callback',{@closeOverlay, 2}); % start at index 2 to match h_image
        end 
    else
        if ~isempty(h_close)
            h_close(end + 1) = uimenu(close_overlays_menu,'Label',...
                overlayName,'Callback',{@closeOverlay, nOverlays + 1});
        else
            h_close(2) = uimenu(close_overlays_menu,'Label',...
                overlayName,'Callback',{@closeOverlay, 2}); % start at index 2 to match h_image
        end 
    end
end

function sort_exts(path1,name1)
    last_nav_dir = path1;
    [~,~,ext] = fileparts(name1);
    for ix = 1:3
        if ~isempty(strfind(extensions{ix},ext))
            type = ix;
            others = setdiff(1:3,ix);
            break;
        end
    end
    extensions = extensions([type,others]);
end

function resizeFigure
    aspect_ratio = (xwidth/x_ax_percent)/(yheight/y_ax_percent);    
    fig_width = aspect_ratio*fig_height;
    figure_pos(1) = max(screen_res(1)+8, screen_mid_x-(.5*fig_width));
    figure_pos(3) = min(screen_res(3)-figure_pos(1)-7, fig_width);
    handles.figure.Position = figure_pos;
end

function success = getSettings
    listing = dir(fullfile(script_path, 'helpers', 'nifti_studio_settings.txt'));
    if isempty(listing)
        disp('No customized settings file found:  missing file ''nifti_studio_settings.txt''.')
        success = false;
        return
    end
    gui_settings = importdata(fullfile(script_path, 'helpers', listing(1).name));
    for ix1 = 1:length(gui_settings)
        ind = strfind(gui_settings{ix1},'=');
        eval([customizable{ix1}, '=[', gui_settings{ix1}(ind+1:end),'];']);
    end
    if ~isdir(last_nav_dir); last_nav_dir = pwd; end %#ok
    success = true;
end

function write_settings
    
    % Change Settings Data / Convert to char array:
    gui_settings = cell(6,1);
    for ix1 = 1:3
        setting1 = eval(customizable{ix1});
        if length(setting1)==1
            gui_settings{ix1} = [customizable{ix1},'=',num2str(setting1),';'];
        else
            gui_settings{ix1} = [customizable{ix1},'=1;'];
        end     
    end
    exts = eval(customizable{4});
    gui_settings{4} = [customizable{4},'={''',exts{1},''';''',exts{2},''';''',exts{3},'''};'];
    gui_settings{5} = [customizable{5},'=''',eval(customizable{5}),''';'];
    gui_settings{6} = [customizable{6},'={''',colorMapStr{1},'''};'];
    gui_settings = char(gui_settings);
    
    % Write as Text File:
    try 
        fileID = fopen(fullfile(script_path, 'helpers', 'nifti_studio_settings.txt'),'w');
        for ix1 = 1:6
            fprintf(fileID,'%s\r\n',gui_settings(ix1,:));
        end
        fclose(fileID);
    catch
        warning('Failed to cache settings.')
    end
end

function revert_defaults(~,~,~)
    % Default display options
    colorbar_on = 1;                                                                                             
    title_on = 1;                                                                                             
    axis_tick_on = 1;       
    set(menu_colorbar,'Checked','on')
    set(menu_axis_tick,'Checked','on')
    set(menu_slice_number,'Checked','on')

    extensions = {'*.img';'*.nii';'*.nii.gz'};                                                                 
    last_nav_dir = first_dir;
    
    % Change colormaps to defaults
    save_selected = selectedImage;
    selectedImage = 1;
    colorMapStr{1} = 'gray'; 
    colormap_callback([], [], 1)
    if length(colorMapStr)>1
        for k = 2:length(colorMapStr)
            selectedImage = k;
            colorMapStr{k} = 'jet';
            colormap_callback([], [], 2)
        end
    end
    selectedImage = save_selected;
    
    refresh_img = true;
    updateImage
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback Functions: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% File Menu:

% Exit Request Callback
function closereq_callback(~,~,~)
    selection = questdlg('Are you sure you want to exit NIfTI Studio?','Exit NIfTI Studio',...
      'Yes','No','Yes');
    switch selection
        case 'Yes'
            write_settings % Save Custom Settings
            % Delete GUI & Figures:
            delete(handles.figure)
        case 'No'
          return
    end
end

% Exit Request No Dialogue:
function closereq_no_dlg(~, ~, ~)
    write_settings % Save Custom Settings
    % Delete GUI & Figures:
    delete(handles.figure)
end

% "Open" Callback:
function status = openNewBackground(~,~,isCallback)
    
    % Get input filename:
    cd(last_nav_dir); 
    [filename, fpath] = uigetfile(extensions, 'Select Image File:',...
        'MultiSelect', 'off'); 
    cd(first_dir);
    
    if fpath ~= 0 
        sort_exts(fpath, filename);
        load_img('char');
        % Save and Update:
        refresh_img = true; 
        updateImage
        updateColormap
        repositionAxes(1)
        status = 1;
    else % "Cancel" selected
        status = 0; 
    end
    
    % Clear Undo/Redo cache:
    undoCache = struct('selectedImage', [], 'action', [], ...
        'orientation', [], 'idx', [], 'color', [], 'alpha', []);
    redoCache = struct('selectedImage', [], 'action', [], ...
        'orientation', [], 'idx', [], 'color', [], 'alpha', []);
    
    if nargin == 3 && isCallback
        % Assure orientation menu reverts to axial:
        for i = 1:length(menu_orientations)
            if ishandle(menu_orientations(i))
                set(menu_orientations(i), 'Checked','off'); 
            end
        end
        if ishandle(menu_orientations(3))
            set(menu_orientations(3),'Checked','on');
        end

        % Clear 
        colorMapStr = {'gray'}; 
    end

end

% "Save" Button Callback:
function save_callback(~, ~, save_type)

    % Check that overlay is selected
    if save_type==2 
  
        if selectedImage==1
            warning('First, select an overlay image to save.')
            return
        end
       
        if isempty(fullPaths{selectedImage})
            % Get save path
            cd(last_nav_dir);
            [filename, PathName] = uiputfile(extensions, 'Specify filename:');
            cd(first_dir);

            % Check whether canceled
            if ~ischar(filename); return; end
            if isdir(PathName); last_nav_dir = PathName; end %#ok
            
            fullPaths{selectedImage} = fullfile(PathName, filename);
        end
    end
            
    % Revert to original orientation
    saveOrientation = [];
    if slice_orientation ~= 3
        saveOrientation = slice_orientation;
        reorient_callback(menu_orientations(3)); 
    end
    
    img.img = permute(imageData{selectedImage},[2,1,3]);
    if untouch_nii
        save_untouch_nii(img, fullPaths{selectedImage})
    else
        save_nii(img, fullPaths{selectedImage})
    end
    
    if save_type == 1
        disp(['Background image successfully saved as ', fullPaths{1}])
    elseif save_type == 2
        disp(['Overlay image successfully saved as ', fullPaths{selectedImage}])
    end
    
    if ~isempty(saveOrientation)
        reorient_callback(menu_orientations(saveOrientation));  
    end
end

% "Save As" Button Callback:
function saveas_callback(~, ~, save_type)
    
    % Check that overlay is selected
    if selectedImage==1 && save_type==2
        warning('First, select an overlay image to save.')
        return
    end
    
    % Get save path
    cd(last_nav_dir);
    [filename, PathName] = uiputfile(extensions, 'Specify filename:');
    cd(first_dir);
   
    % Check whether canceled
    if ~ischar(filename); return; end
    if isdir(PathName); last_nav_dir = PathName; end %#ok
            
    % Revert to original orientation
    saveOrientation = [];
    if slice_orientation ~= 3 
        saveOrientation = slice_orientation;
        reorient_callback(menu_orientations(3)); 
    end

    % Save:
    img.img = permute(imageData{selectedImage},[2,1,3]);
    fullPaths{selectedImage} = fullfile(PathName, filename);
    if untouch_nii
        save_untouch_nii(img, fullPaths{selectedImage})
    else
        save_nii(img, fullPaths{selectedImage})
    end

    % Save types:
    if save_type == 1 % save background
        
        disp(['Background image successfully saved as ',fullPaths{selectedImage}])

        % Change Figure Window Title to Reflect New Filename:
        window_name = ['NIfTI Studio:    ',filename];
        set(handles.figure,'Name',window_name)
            
    elseif save_type == 2 % save current overlay
        
        disp(['Overlay image successfully saved as ', fullPaths{selectedImage}])
        
    end
    
    if ~isempty(saveOrientation)
        reorient_callback(menu_orientations(saveOrientation));  
    end
   
end

function title_toggle_callback(~,~,~)
    if title_on
        set(menu_slice_number,'Checked','off');
    else
        set(menu_slice_number,'Checked','on');
    end
    title_on = ~title_on;
    refresh_img = true; 
    updateImage 
    repositionAxes
end

function axis_tick_callback(~,~,~)
    if axis_tick_on
        set(menu_axis_tick,'Checked','off');
    else
        set(menu_axis_tick,'Checked','on');
    end
    axis_tick_on = ~axis_tick_on;
    refresh_img = true; 
    repositionAxes
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Image Orientation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reorient Image (axial, coronal, sagittal):
function reorient_callback(hObject, ~, ~)
    reorient_type = get(hObject,'Label');
    
    % Uncheck all menu items
    if ~strcmp(reorient_type, launch_3d_string)
        for i = 1:length(menu_orientations)
            set(menu_orientations(i),'Checked','off'); 
        end
    end
    
    % Add to undoCache
    if ~strcmp(reorient_type, launch_3d_string) && nargin < 3
        cacheState(selectedImage, 'reorient', slice_orientation, [], [],[])
    end
    
    switch reorient_type
%         case 'Multiplanar'
%             set(menu_orientations(1),'Checked','on')
        case 'Coronal' % (1)
            set(menu_orientations(1),'Checked','on')
            if slice_orientation~=1
                if slice_orientation==2
                    dimperm = [1,3,2];
                    y_slice = curr_slice;
                elseif slice_orientation==3
                    dimperm = [3,2,1];
                    z_slice = curr_slice;
                end       
                slice_orientation = 1;
                if isempty(x_slice)
                    curr_slice = [];
                else
                    curr_slice = x_slice;
                end
            else
                return;
            end
            xwidth = dim(2)*pixdim(2);
            yheight = dim(4)*pixdim(4);
        case 'Sagittal' % (2)
            set(menu_orientations(2),'Checked','on')
            if slice_orientation~=2
                if slice_orientation==1
                    dimperm = [1,3,2];
                    x_slice = curr_slice;
                elseif slice_orientation==3
                    dimperm = [3,1,2];
                    z_slice = curr_slice;
                end
                slice_orientation = 2;
                if isempty(y_slice)
                    curr_slice = [];
                else
                    curr_slice = y_slice;
                end
            else
                return;
            end
            xwidth = dim(3)*pixdim(3);
            yheight = dim(4)*pixdim(4);
        case 'Axial' % (3)
            set(menu_orientations(3),'Checked','on')
            if slice_orientation~=3
                if slice_orientation==1
                    dimperm = [3,2,1]; 
                    x_slice = curr_slice;
                elseif slice_orientation==2
                    dimperm = [2,3,1];
                    y_slice = curr_slice;
                end
                slice_orientation = 3;  
                if isempty(z_slice)
                    curr_slice = [];
                else
                    curr_slice = z_slice;
                end
            else
                return;
            end
            xwidth = dim(2)*pixdim(2);
            yheight = dim(3)*pixdim(3);
            
        case {launch_3d_string, launch_mosaic_string}
            
            % Change mouse pointer for busy
            pointer = get(handles.figure, 'Pointer');
            set(handles.figure, 'Pointer', 'watch')
            drawnow nocallbacks
            
            % Change orientations if need be
            switch slice_orientation
                case 1
                    save_orientation = 1;
                    reorient_callback(menu_orientations(3)) % revert to original orientation
                    change_orientation = 1;
                case 2
                    save_orientation = 2;
                    reorient_callback(menu_orientations(3)) % revert to original orientation
                    change_orientation = 1;
                case 3
                    save_orientation = 3;
                    change_orientation = 0;
            end
            
            % Collect background image
            background_img = img; 
            background_img.img = permute(imageData{1},[2,1,3]);
            
            % Threshold background image
            background_thresh_3d = .6; % filter out bottom 60%

            % Combine overlays
            if numel(imageData) > 1
                combined_overlay = zeros(size(imageData{1}));
                empty_idx = true(size(imageData{1}));
                for i = 2:length(imageData)
                    if ~isempty(imageData{i})
                        combined_overlay(empty_idx) = imageData{i}(empty_idx);
                        empty_idx(((combined_overlay(:)~=0) + (~isnan(combined_overlay(:))))==2) = false;
                    end
                end
                overlay1_img = img;
                overlay1_img.img = permute(combined_overlay,[2,1,3]);
            end
            
            % Generate new plot
            if strcmp(reorient_type, launch_3d_string)
                
                if numel(imageData) > 1
                    [~] = nifti_studio_3D('background',background_img,...
                        'ROI',overlay1_img,'roi_type',1, ... %'titles','Overlay 1',...
                        'vox_thresh', quantile(background_img.img(:), background_thresh_3d));
                else
                    [~] = nifti_studio_3D('background',background_img,... % 'titles','Background', ...
                        'vox_thresh', quantile(background_img.img(:), background_thresh_3d));
                end
                
            elseif strcmp(reorient_type, launch_mosaic_string)
                                
                if numel(imageData) > 1
                    [~] = nifti_studio_mosaic('background',background_img,...
                        'overlay',overlay1_img',...
                        'title','NIfTI Studio Mosaic',...
                        'dimension',save_orientation);
                else
                    [~] = nifti_studio_mosaic('background',background_img, ...
                        'title','NIfTI Studio Mosaic',...
                        'dimension',save_orientation); 
                end
                
            end
            
            % Return to Original Orientation if Changed:
            if change_orientation
                reorient_callback(menu_orientations(save_orientation),[]) % have to have two inputs, see line 1047
            end
            
            % Change mouse pointer back
            set(handles.figure, 'Pointer', pointer)
            drawnow
            return
            
    end
    
    for i = 1:numel(imageData)
        imageData{i} = permute(imageData{i}, dimperm);
        alphaData{i} = permute(alphaData{i}, dimperm);
    end
    [xdim, ydim, zdim] = size(imageData{1}); 
    voxSize = voxSize(dimperm);
    origin = origin(dimperm);
    
    % Determine possible combinations of x,y indices:
    poss_ind = zeros(xdim*ydim,2); count = 0;
    for ix = 1:xdim
        for jx = 1:ydim
            count = count + 1;
            poss_ind(count,:) = [ix,jx];
        end
    end

    xind = xdim:-1:1;
    if slice_orientation==2
        yind = ydim:-1:1; 
    else
        yind = 1:ydim;
    end
    middle_slice = round(zdim/2);
    if isempty(curr_slice); curr_slice = middle_slice; end
    
    % Update Image:
    refresh_img = false; 
    updateImage
    repositionAxes
    resizeFigure;
end

function changeUnits(~,~,~)
    if strcmp(units, 'physical')
        units = 'voxel';
        set(h_units(1), 'Checked','off'); set(h_units(2), 'Checked','on')
    else
        units = 'physical';
        set(h_units(1), 'Checked','on'); set(h_units(2), 'Checked','off')
    end
    % Update Image:
    refresh_img = false; 
    updateImage
    repositionAxes
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Edit Underlying Data (Background Figure) Colorspec, Colormap, Colorbar, Auto-Scale:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Colorbar Callback:
function colorbar_callback(~, ~, ~)
    if colorbar_on
        set(menu_colorbar,'Checked','off');
    else
        set(menu_colorbar,'Checked','on');
    end
    colorbar_on = ~colorbar_on;
    refresh_img = true; 
    updateImage
    repositionAxes
end

%% Color menu callbacks

% Colormap Callback:
function colormap_callback(~, ~, cmap_num)
    colorMapStr{selectedImage} = colormap_opts{cmap_num};
    % Update checkmarks
    for k = 1:length(colormap_opts)
       if k ~= cmap_num
          set(h_colormap_menu(k), 'Checked', 'off')
       else
          set(h_colormap_menu(k), 'Checked', 'on') 
       end
    end
    % Update graphics
    updateColormap 
end

% Update Colormap
function updateColormap
    colormaps{selectedImage} = eval([colorMapStr{selectedImage},'(',num2str(colormap_n),')']);
    colormap(handles.axes{selectedImage}, colormaps{selectedImage});
end

% Change Color Limits
function colorscale_callback(~, ~, cscale_num) 
    colorscaleType{selectedImage} = cscale_num;
    switch colorscaleType{selectedImage}
        case 1 % 'Slice Min/Max' 
%             cmin{selectedImage} = NaN;
%             cmax{selectedImage} = NaN;
        case 2 % 'Global Min/Max'
            cmin{selectedImage} = min(imageData{selectedImage}(:));
            cmax{selectedImage} = max(imageData{selectedImage}(:));
            if cmin{selectedImage}==cmax{selectedImage}
                cmax{selectedImage} = cmin{selectedImage} + 1;
            end
        case 3 % 'Custom'
            cmin1 = min(imageData{selectedImage}(:));
            cmax1 = max(imageData{selectedImage}(:));
            prompt = {['Specify range intensity values for color axis: ',...
                char(10),'Min:'],'Max:'}; %#ok
            dlg_title = 'Input CAxis Limits'; num_lines = [1,20;1,20]; 
            defaultans = {num2str(cmin1),num2str(cmax1)};
            answer1 = inputdlg(prompt,dlg_title,num_lines,defaultans);
            if ~isempty(answer1)
                cmin2 = str2double(answer1{1});
                cmax2 = str2double(answer1{2});
            else
                return;
            end
            if ~isempty(cmin2) && ~isempty(cmax2) && cmin2~=cmax2
                cmin{selectedImage} = cmin2;
                cmax{selectedImage} = cmax2;
            else
                cmin{selectedImage} = cmin1;
                cmax{selectedImage} = cmax1;
                disp('Invalid color axis specification.')
            end
    end
    
    % Update checkmarks
    for k = 1:3
        if k ~= cscale_num
            set(h_colorscale_menu(k), 'Checked','off')
        else
            set(h_colorscale_menu(k), 'Checked','on')
        end
    end
    
    % Update graphics
    refresh_img = 0; 
    updateImage
end

% Set transparency (alpha):
function opacity_callback(~, ~, alpha_num)
    % Set alpha value
    set(opacity_menu, 'Visible','on')
    
    alphaValue{selectedImage} = alpha_values(alpha_num);
    alphaData{selectedImage} = zeros(size(imageData{selectedImage})); 
    use_idx = ((imageData{selectedImage}(:)~=0) + (~isnan(imageData{selectedImage}(:))))==2;
    alphaData{selectedImage}(use_idx) = alphaValue{selectedImage};

    % Update checkmarks
    for k = 1:length(alpha_opts)
        if k ~= alpha_num
            set(h_opacity_menu(k), 'Checked','off')
        else
            set(h_opacity_menu(k), 'Checked','on')
        end
    end
    
    % Update graphics
    updateImage
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Images Callback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeSelection(~,~,whichImage)
    
    % Adjust Select menu checkmarks
    selectedImage = whichImage;
    for ix = 1:length(h_image)
        if ishandle(h_image(ix))
            set(h_image(ix),'Checked','off')
        end
    end
    set(h_image(selectedImage),'Checked','on')
    
    % Update colormap checkmarks
    cmap_num = find(strcmp(colorMapStr{selectedImage}, colormap_opts));
    for k = 1:length(colormap_opts)
       if k ~= cmap_num
          set(h_colormap_menu(k), 'Checked', 'off')
       else
          set(h_colormap_menu(k), 'Checked', 'on') 
       end
    end
    
    % Colorscale Type Checkmarks:
    for k = 1:3
        if k ~= colorscaleType{selectedImage}
            set(h_colorscale_menu(k), 'Checked','off')
        else
            set(h_colorscale_menu(k), 'Checked','on')
        end
    end
    
    % Adjust opacity checkmarks
    alpha_num = find(alphaValue{selectedImage}==alpha_values);
    for k = 1:length(alpha_opts)
        if k ~= alpha_num
            set(h_opacity_menu(k), 'Checked','off')
        else
            set(h_opacity_menu(k), 'Checked','on')
        end
    end

    % Disable opacity menu if main image
    if whichImage == 1
       set(opacity_menu, 'Visible','off') 
    else
       set(opacity_menu, 'Visible','on') 
    end
    
    % Switch colorbar to new axes  
    updateColorbar(1)
end

%% Tools Callbacks:

function adjustToolsCheckmarks(whichTool)
    set(tool_crosshair,'Checked','off')
    set(tool_draw,'Checked','off')
%     set(tool_zoom,'Checked','off')
    set(tool_pan,'Checked','off')
    switch whichTool
        case 1
            set(tool_crosshair,'Checked','on')
        case 2
            set(tool_draw,'Checked','on')
%         case 3
%             set(tool_zoom,'Checked','on')
        case 4
            set(tool_pan,'Checked','on')
    end 
end

function disableInteractiveModeHijack(whichTool)
    % Disable Matlab hijacking keyboard shortcuts: 
    %   https://undocumentedmatlab.com/blog/enabling-user-callbacks-during-zoom-pan
    
    hManager = uigetmodemanager(handles.figure);
    try 
        set(hManager.WindowListenerHandles, 'Enable', 'off');  % HG1
    catch
        [hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2
    end
    
    % Clear the new interactive mode keypress callbacks
    set(handles.figure, 'WindowKeyPressFcn', []);
    if whichTool
        set(handles.figure, 'KeyPressFcn', []);
    end
    
    % Add back callbacks
    set(handles.figure,'WindowKeyPressFcn',@keypress_callback);
    set(handles.figure,'WindowScrollWheelFcn',@scroll_zoom_callback);
    
end

% Tools menu callbacks:
function dcm_obj = crosshair_callback(~,~,~)
%     zoom off
    dcm_obj = datacursormode(handles.figure);
    set(dcm_obj, 'Enable', 'on', 'UpdateFcn', @dataCursorCallback)

    draw_on = false;
    pan_on = false;
    
    % Disable callback hijack
    disableInteractiveModeHijack(1)
    
    % Adjust checkmarks
    adjustToolsCheckmarks(1)
end

% Callback for formatting Data Cursor text box
function output_txt = dataCursorCallback(~, event_obj)
    % event_obj    Object containing event data structure
    % output_txt   Data cursor text
    
    % Get position
    pt = get(event_obj, 'Position');
    
    % Format text:
    switch slice_orientation
        case 1 % coronal
            output_txt = ['[',...
                num2str(round((yind(pt(1))-origin(2)) * voxSize(2),2)), ', ',...
                num2str(round((curr_slice - origin(3)) * voxSize(3),2)),', ',...
                num2str(round((xind(pt(2))-origin(1)) * voxSize(1),2)), '] ',...
                physical_units, '\n[',...
                num2str(yind(pt(1))), ', ',...
                num2str(curr_slice), ', ', ...
                num2str(xind(pt(2))),'] vox',...
                '\nBackground:  ',num2str(round(imageData{1}(xind(pt(2)),yind(pt(1)),curr_slice),2))];
        case 2 % sagittal
            output_txt = ['[',...
                num2str(round((curr_slice - origin(3)) * voxSize(3),2)),', ',...
                num2str(round((yind(pt(1))-origin(2)) * voxSize(2),2)), ', ',...
                num2str(round((xind(pt(2))-origin(1)) * voxSize(1),2)), '] ',...
                physical_units, '\n[',...
                num2str(curr_slice), ', ',...
                num2str(yind(pt(1))), ', ',...
                num2str(xind(pt(2))), ']  vox',...
                '\nBackground:  ',num2str(round(imageData{1}(xind(pt(2)),yind(pt(1)),curr_slice),2))];
        case 3 % axial
            output_txt = ['[',...
                num2str(round((yind(pt(1))-origin(2)) * voxSize(2),2)), ', ',...
                num2str(round((xind(pt(2))-origin(1)) * voxSize(1),2)), ', ',...
                num2str(round((curr_slice - origin(3)) * voxSize(3),2)),'] ',...
                physical_units, '\n[',...
                num2str(yind(pt(1))), ', ',...
                num2str(xind(pt(2))), ', ',...
                num2str(curr_slice),'] vox',...
                '\nBackground:  ',num2str(round(imageData{1}(xind(pt(2)),yind(pt(1)),curr_slice),2))];
    end
    
    % Add overlay text:
    if numel(imageData) > 1
        for i = 2:length(imageData)
            if ~isempty(imageData{i})
                output_txt = [output_txt, '\nOverlay ',num2str(i-1),':  ',...
                        num2str(round(imageData{i}(xind(pt(2)),yind(pt(1)),curr_slice),2))]; %#ok
            end
        end
    end
    output_txt = sprintf(output_txt);
end

function draw_callback(~,~,~)
    % Turn off and clear any remaining data tips
    datacursormode(handles.figure, 'off');
    delete(findall(handles.figure,'Type','hggroup'));
    
%     zoom off
    draw_on = true;
    pan_on = false;
    set(handles.figure,'pointer','cross'); 
    
    % Adjust checkmarks
    adjustToolsCheckmarks(2)
    
    % Make sure current figure
%     set(0, 'currentfigure', handles.figure); 
%     set(handles.figure, 'currentaxes', handles.axes{selectedImage}); 
end

% function zoom_callback(~,~,~)
%     % Turn off and clear any remaining data tips
%     datacursormode(handles.figure, 'off');
%     delete(findall(handles.figure,'Type','hggroup'));
%     
%     zoom on; 
%     draw_on = false;
%     pan_on = false;
%     
%     % Set pointer to hand (pan on only adjusts after motion):
%     % https://undocumentedmatlab.com/blog/undocumented-mouse-pointer-functions
%     setptr(handles.figure, 'glassplus');
%     
%     % Disable callback hijack
%     disableInteractiveModeHijack(3)
% 
%     % Adjust checkmarks
%     adjustToolsCheckmarks(3)
%     
% end

function pan_callback(~,~,~)
    % Turn off and clear any remaining data tips
    datacursormode(handles.figure, 'off');
    delete(findall(handles.figure,'Type','hggroup'));
    
%     zoom off;
    draw_on = false;
    pan_on = true;
    
    % Set pointer to hand (pan on only adjusts after motion):
    % https://undocumentedmatlab.com/blog/undocumented-mouse-pointer-functions
%     setptr(handles.figure, 'hand');
    set(handles.figure, 'Pointer', 'fleur');
    
    % Adjust checkmarks
    adjustToolsCheckmarks(4)
end

%% Undo / Redo

function cacheState(selectedImage, action, orientation, idx, color, alpha)

    % In case currently iterating undoCallback, new action resets forward
    % memory:
    redoCache = [];
    if numel(undoCache) > undoNum
        undoCache(undoNum + 1 : end) = [];
    end
    
    % Cache previous state
    if undoNum == 0
        undoCache(1) = struct('selectedImage', selectedImage, ...
            'action', action, 'orientation', orientation, 'idx', idx, ...
            'color', color, 'alpha', alpha);
    else
        undoCache(end + 1) = struct('selectedImage', selectedImage, ...
            'action', action, 'orientation', orientation, 'idx', idx, ...
            'color', color, 'alpha', alpha);
    end
    
    % If limit reached, delete first
    if length(undoCache) > undoLimit
        undoCache(1) = [];
    end
    
    undoNum = min(undoNum + 1, undoLimit);
    
end

% Undo Callback:
function undoCallback(~, ~, ~)

    if undoNum==0; return; end
  
    % Update or initialize redoCache
    if undoNum == length(undoCache) % first undo
        if undoCache(undoNum).selectedImage==1
            redoCache = struct('selectedImage', undoCache(undoNum).selectedImage, ...
                'action', undoCache(undoNum).action, ...
                'orientation', slice_orientation, ...
                'idx', undoCache(undoNum).idx, ...
                'color', imageData{undoCache(undoNum).selectedImage}(undoCache(undoNum).idx), ...
                'alpha', []);
        else
            redoCache = struct('selectedImage', undoCache(undoNum).selectedImage, ...
                'action', undoCache(undoNum).action, ...
                'orientation', slice_orientation, ...
                'idx', undoCache(undoNum).idx, ...
                'color', imageData{undoCache(undoNum).selectedImage}(undoCache(undoNum).idx), ...
                'alpha', alphaData{undoCache(undoNum).selectedImage}(undoCache(undoNum).idx));
        end
    else
        if undoCache(undoNum).selectedImage==1
            redoCache(end+1) = struct('selectedImage', undoCache(undoNum).selectedImage, ...
                'action', undoCache(undoNum).action, ...
                'orientation', slice_orientation, ...
                'idx', undoCache(undoNum).idx, ...
                'color', imageData{undoCache(undoNum).selectedImage}(undoCache(undoNum).idx), ...
                'alpha', []);
        else
            redoCache(end+1) = struct('selectedImage',undoCache(undoNum).selectedImage, ...
                'action', undoCache(undoNum).action, ...
                'orientation', slice_orientation, ...
                'idx', undoCache(undoNum).idx, ...
                'color', imageData{undoCache(undoNum).selectedImage}(undoCache(undoNum).idx), ...
                'alpha', alphaData{undoCache(undoNum).selectedImage}(undoCache(undoNum).idx));
        end
    end
    
    % Revert to most recent undo
    if strcmp(undoCache(undoNum).action, 'draw')
        imageData{undoCache(undoNum).selectedImage}(undoCache(undoNum).idx) = undoCache(undoNum).color;
        if undoCache(undoNum).selectedImage > 1
            alphaData{undoCache(undoNum).selectedImage}(undoCache(undoNum).idx) = undoCache(undoNum).alpha;
        end
    elseif strcmp(undoCache(undoNum).action, 'reorient')
        reorient_callback(menu_orientations(undoCache(undoNum).orientation), [], 1)
    end

    undoNum = undoNum - 1;
    updateImage
end

% Redo Callback:
function redoCallback(~,~,~)
    
    if numel(redoCache) == 0; return; end
    
    % Move forward to most recent redo
    if strcmp(redoCache(end).action, 'draw')
        imageData{redoCache(end).selectedImage}(redoCache(end).idx) = redoCache(end).color;
        if redoCache(end).selectedImage > 1
            alphaData{redoCache(end).selectedImage}(redoCache(end).idx) = redoCache(end).alpha;
        end
    elseif strcmp(redoCache(end).action, 'reorient')
        reorient_callback(menu_orientations(redoCache(end).orientation), [], 1)
    end
    
    redoCache(end) = [];
    undoNum = min(undoNum + 1, length(undoCache));
    updateImage
end

% Clear any overlay edits from undo/redo cache
function purgeUndoRedo(closeWhich)
    
    i = 1;
    while i <= length(undoCache)
       if (~isempty(undoCache(i).selectedImage) && undoCache(i).selectedImage == closeWhich) || ...
            (~isempty(undoCache(i).action) && ~strcmp(undoCache(i).action,'reorient')) % don't purge 'reorient' actions
          undoCache(i) = []; 
          undoNum = undoNum - 1;
       else 
           i = i + 1;
       end
    end
    
    i = 1;
    while i <= length(redoCache)
        if (~isempty(redoCache(i).selectedImage) && redoCache(i).selectedImage == closeWhich) || ...
            (~isempty(redoCache(i).action) && ~strcmp(redoCache(i).action,'reorient')) % don't purge 'reorient' actions
            redoCache(i) = []; 
        else 
            i = i + 1;
        end
    end
    
    undoNum = min(max(undoNum, 0), length(undoCache));
end

%% Other Editing

% Go to Slice Callback:
function goToSlice(~,~,~)
    prompt = {'Specify slice: '};
    dlg_title = 'Go to Slice'; num_lines = [1,15;]; defaultans = {num2str(zdim)};
    answer1 = inputdlg(prompt,dlg_title,num_lines,defaultans);
    if isempty(answer1); return; end
    slice_spec = str2double(answer1(1)); 
    if slice_spec <= zdim && slice_spec > 0
        curr_slice = slice_spec;
    else
        disp(['Slice specification out of range. Slice # must be between 1 and ',...
            num2str(zdim),'.']);
        return;
    end
    updateImage
end

function goToOrigin(~,~,~)
    % Go to origin slice
    curr_slice = round(origin(3));
    updateImage;
    
    % Clear data tips
    delete(findall(handles.figure,'Type','hggroup'));
    
    % Create datatip
    % https://undocumentedmatlab.com/articles/controlling-plot-data-tips
    dcm_obj = crosshair_callback;  % switch to crosshair tool
    hDatatip = dcm_obj.createDatatip(handles.images{1});
    set(hDatatip,'UIContextMenu',get(dcm_obj,'UIContextMenu'));
    set(hDatatip,'Position',[yind(round(origin(2))), xind(round(origin(1)))])
end

% Add Border Button Callback:
function add_border_callback(~,~,~)
    prompt = {'# voxels: '};
    dlg_title = 'Border'; num_lines = [1,25];
    defaultans = {'0'};
    answer1 = inputdlg(prompt,dlg_title,num_lines,defaultans,'on');
    if ~isempty(answer1)
        num_voxels_border = str2double(answer1);
        % Assure max is not beyond image dimensions
        size(imageData{selectedImage}(xind,yind,curr_slice))
        num_voxels_border = min(num_voxels_border, ...
            min(size(imageData{selectedImage}(xind,yind,curr_slice)))-1);
        if num_voxels_border ~= 0
            % Get border indices
            slice_data = false(size(imageData{selectedImage}(xind,yind,curr_slice)));
            slice_data(1:num_voxels_border,:) = true;
            slice_data(end-num_voxels_border:end,:) = true;
            slice_data(:,1:num_voxels_border) = true;
            slice_data(:,end-num_voxels_border:end) = true;
            [x, y] = ind2sub(size(slice_data), find(slice_data(:)));
            idx_border = sub2ind(size(imageData{selectedImage}), x, y, repmat(curr_slice, [numel(x), 1]));
            
            % Add to undoCache
            if selectedImage==1
                cacheState(selectedImage, 'draw', [], idx_border, ...
                    imageData{selectedImage}(idx_border), [])
            else
                cacheState(selectedImage, 'draw', [], idx_border, ...
                    imageData{selectedImage}(idx_border), ...
                    alphaData{selectedImage}(idx_border))
            end
            
            % Add border
            imageData{selectedImage}(idx_border) = draw_color;
            curr_drawing = imageData{selectedImage}(:,:,curr_slice);
            if selectedImage > 1
                if erase_draw
                    alphaData{selectedImage}(idx_border) = NaN;
                else
                    alphaData{selectedImage}(idx_border) = alphaValue{selectedImage};
                end
            end
            
            updateImage
        end
    end
end    

% "Apply to Whole Slice" Button Callback:
function apply2whole_slice_callback(~,~,~)

    % Get fill indices:
    slice = false(size(imageData{selectedImage}));
    slice(xind,yind,curr_slice) = true;
    idx_fill = find(slice);
    
    % Add to undoCache:
    if selectedImage==1
        cacheState(selectedImage, 'draw', [], idx_fill, ...
            imageData{selectedImage}(idx_fill), [])
    else
        cacheState(selectedImage, 'draw', [], idx_fill, ...
            imageData{selectedImage}(idx_fill), ...
            alphaData{selectedImage}(idx_fill))
    end
    
    % Fill slice w/ color
    imageData{selectedImage}(idx_fill) = draw_color;
    if selectedImage > 1
        if erase_draw
            alphaData{selectedImage}(idx_fill) = NaN;
        else
            alphaData{selectedImage}(idx_fill) = alphaValue{selectedImage};
        end
    end
    updateImage
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create User Interface Callbacks (clicks, arrow keys, mouse scrolls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Right-click Pan
function motion_callback_pan(~,~,~)
    set(handles.figure, 'WindowButtonUpFcn', @cursor_unclick_callback_pan);
    
    % Get mouse location:
    pt = get(handles.axes{1}, 'CurrentPoint');
    pt = round(pt(1,1:2));
    
    % Check that points are valid:
    pt_check = sum([(pt>0),(pt<=[length(yind),length(xind)])]);
    if pt_check~=4
        return
    end
    
    % Adjust axes based on pt comparison to save_points:
    diff_points_limit = 7; % higher #'s approximately mean faster/more sensitive panning
    diff_points = [pt(2), pt(1)] - save_points;
    diff_points(1) = max(diff_points(1),-diff_points_limit);
    diff_points(1) = min(diff_points(1),diff_points_limit);
    diff_points(2) = max(diff_points(2),-diff_points_limit);
    diff_points(2) = min(diff_points(2),diff_points_limit);
    
    % Determine initial zoom for panning (sagittal seems to require slightly less)
    %   -from fully zoomed out, panning requires slight zoom to enable (larger int means larger initial zoom)
    if slice_orientation==2 % sagittal
        pan_limits = 3;
    else
        pan_limits = 5; 
    end
    
    % Determine whether to stop panning in a given direction:
    if xmax == xdim && xmin >= pan_limits && diff_points(1)>0 % reached xmax limit
        diff_points(1) = 0;
    end
    if xmin == 1 && xmax <= (xdim - pan_limits) && diff_points(1)<0 % reached xmin limit
        diff_points(1) = 0;
    end
        
    xmax = min(xmax + diff_points(1), xdim);
    xmin = max(xmin + diff_points(1), 1);
    xind = xmax:-1:xmin;
    
    if slice_orientation==2
        if ymax == ydim && ymin >= pan_limits && diff_points(2)>0 % reached ymax limit
            diff_points(2) = 0;
        end
        if ymin == 1 && ymax <= (ydim - pan_limits) && diff_points(2)<0 % reached ymin limit
            diff_points(2) = 0;
        end
    
        ymax = min(ymax + diff_points(2), ydim);
        ymin = max(ymin + diff_points(2), 1);
        yind = ymax:-1:ymin;
    else
        if ymax == ydim && ymin >= pan_limits && diff_points(2)<0 % reached ymax limit
            diff_points(2) = 0;
        end
        if ymin == 1 && ymax <= (ydim - pan_limits) && diff_points(2)>0 % reached ymin limit
            diff_points(2) = 0;
        end
        ymax = min(ymax - diff_points(2), ydim);
        ymin = max(ymin - diff_points(2), 1);
        yind = ymin:ymax;
    end
    
    save_points = [pt(2),pt(1)];
    ax_xlim = [1,numel(yind)]; ax_ylim = [1,numel(xind)];
    refresh_img = true; 
    updateImage
    repositionAxes
end

% Arrow Keys 
function keypress_callback(varargin)
    % Clear any data tips
    delete(findall(handles.figure,'Type','hggroup'));
    
    keypress = varargin{2}.Key;
    switch keypress
        case 'uparrow' 
            if curr_slice < zdim
                curr_slice = curr_slice + 1; 
                updateImage; 
            end
        case 'downarrow'
            if curr_slice > 1
                curr_slice = curr_slice - 1; 
                updateImage;
            end
        case 'escape'
            closereq_callback
        case 'u'
            undoCallback
        case 'r'
            redoCallback
        case 'g'
            goToSlice
        case 'o'
            goToOrigin
        case 'c'
            crosshair_callback;
        case 'd'
            draw_callback
%         case 'z'
%             zoom_callback
        case 'p'
            pan_callback
    end     
end

% Unclick after pan
function cursor_unclick_callback_pan(~,~,~)
    set(handles.figure, 'WindowButtonMotionFcn', []);
    set(handles.figure, 'WindowButtonUpFcn', []);
    save_points = [0,0];
    
    % Reset pointer symbol from panning 'fleur' to given setting:
    if draw_on; set(handles.figure,'pointer','cross'); end     
end

% Cursor click callback:
function cursor_click_callback(src,~,~)
    if strcmpi(get(gco,'type'),'image') 
        
        % Get mouse location
        pt = get(handles.axes{1}, 'CurrentPoint');
        pt = round(pt(1,1:2));
        
        if pt(1) > 0 && pt(2) > 0
            
            if strcmp(src.SelectionType, 'alt') || pan_on % Right Click, or Ctrl + Left Click (pan)
                
                set(handles.figure, 'WindowButtonMotionFcn', @motion_callback_pan);
                save_points = [pt(2), pt(1)];
                set(handles.figure,'pointer','fleur'); % setptr(handles.figure, 'hand');
                
            elseif strcmp(src.SelectionType, 'normal') % Left Click (crosshair or draw)

                if draw_on
                    
                    in_motion = false;
                    
                    % Set Callbacks
                    set(handles.figure, 'WindowButtonUpFcn', @cursor_unclick_callback);
                    set(handles.figure, 'WindowButtonMotionFcn', @cursor_motion_callback);
                    
                    % Add color:
                    if shape==0 % non-shape trace
                        
                        save_points = [pt(2),pt(1)]; 
                        
                    elseif shape==1 || shape==2 % circle, rectangle
                        
                        shape_ind(1,:) = [xind(pt(2)), yind(pt(1))];
                        
                    elseif shape==3 % sphere
                        
                        s_center = [xind(pt(2)),yind(pt(1)),curr_slice];
                        s_center = s_center.*voxSize; % convert voxel to physical units
                        Min = s_center - sphere_radius;
                        Min = fix(max(voxSize,Min));
                        Max = s_center + sphere_radius;
                        Max = ceil(min([xdim, ydim, zdim].*voxSize,Max));
                        curr_drawing = zeros(size(imageData{selectedImage})); % used for propagating last drawing
                        for ix = Min(1):voxSize(1):Max(1) %Min(1):Max(1)
                            for jx = Min(2):voxSize(2):Max(2)  %Min(2):Max(2)
                                for kx = Min(3):voxSize(3):Max(3)  %Min(3):Max(3)
                                    if sqrt(sum(([ix,jx,kx] - s_center).^2 )) <= sphere_radius
                                        x_coord = ceil(ix/voxSize(1));
                                        y_coord = ceil(jx/voxSize(2));
                                        z_coord = ceil(kx/voxSize(3));
                                        if x_coord<1 || y_coord <1 || z_coord<1 || x_coord>xdim || y_coord>ydim || z_coord>zdim
                                            continue
                                        end
                                        curr_drawing(x_coord,y_coord,z_coord) = draw_color;
                                    end
                                end
                            end
                        end
                        idx_draw = find(curr_drawing(:)>0);
                    end 
                end
            end
        end
    end
end

% Mouse motion after click 
%   -this function interactively adds/removes shape drawings (circle, 
%   rectangle) as the mouse moves
%   -it also stores the traced outline for later interpolation if non-shape
%   drawing
%   -not used for spheres
function cursor_motion_callback(~,~,~)
    set(handles.figure, 'WindowButtonUpFcn', @cursor_unclick_callback);
    
    % Get mouse location:
    pt = get(handles.axes{1}, 'CurrentPoint'); 
    pt = round(pt(1,1:2));
    
    % Check Points:
    pt_check = sum([(pt>0), (pt<=[length(yind),length(xind)])]);
    if pt_check~=4; return; end
    
    in_motion = true;
    
    if shape==0 % non-shape tracing
        
        save_points = [save_points; pt(2), pt(1)];
        %     updateImage % this turns on tracing of drawings, but decreases performance slightly
        
    elseif shape==1 || shape==2 % Interactively add/remove circles or rectangle
        
        % Erase previous shape if still in-motion
        if ~isempty(idx_draw)
            imageData{selectedImage}(idx_draw) = prev_shape_colors;
            if selectedImage > 1
                alphaData{selectedImage}(idx_draw) = prev_alpha;
            end
        end
        
        shape_ind(2,:) = [xind(pt(2)), yind(pt(1))];
        
        if shape==1 % Circle
            
            center = sum(shape_ind,1).*.5;
            radius = .5*sqrt(sum(diff(shape_ind).^2));        
            [rr,cc] = meshgrid(1:ydim, 1:xdim);
            C = sqrt((rr-center(2)).^2+(cc-center(1)).^2)<=radius;
            [ix, iy] = ind2sub([xdim, ydim], find(C(:)));
            shape_coords = [ix, iy];
            
        elseif shape==2 % Rectangle

            [x,y] = meshgrid(min(shape_ind(:,1)):max(shape_ind(:,1)),...
                min(shape_ind(:,2)):max(shape_ind(:,2)));
            shape_coords = [x(:),y(:)];

        end
        
        % Draw new shape
        idx_draw = sub2ind(size(imageData{selectedImage}), ...
            shape_coords(:,1), ...
            shape_coords(:,2), ...
            repmat(curr_slice, [size(shape_coords,1),1]));
        
        % Save for subsequent erase if still in-motion
        prev_shape_colors = imageData{selectedImage}(idx_draw);
        imageData{selectedImage}(idx_draw) = draw_color;
        if selectedImage > 1
            prev_alpha = alphaData{selectedImage}(idx_draw);
            if erase_draw
                alphaData{selectedImage}(idx_draw) = NaN;
            else
                alphaData{selectedImage}(idx_draw) = alphaValue{selectedImage};
            end
        end
        
        % Update image to show current shape
        updateImage
    end
end

% Cursor unclick callback
function cursor_unclick_callback(~,~,~)
    
    % Remove callback functions
    set(handles.figure, 'WindowButtonMotionFcn', []);
    set(handles.figure, 'WindowButtonUpFcn', []);
 
    if shape==0 % non-shape tracing 

        if in_motion && size(save_points,1) > 2 % interpolate tracing
            
            % Find coordinates inside traced drawing:
            warning('off','all'); 
            inside = inpolygon(poss_ind(:,1),poss_ind(:,2),...
                save_points(:,1),save_points(:,2)); 
            warning('on','all')
            save_points = unique([save_points(:,1:2); poss_ind(inside,:)], 'rows');
            
            % Remove any out-of-bounds coordinates
            save_points(((save_points(:,1)>max(xind))+(save_points(:,2)>max(yind)))>0,:) = [];
            
        end
        
        % Find linear indices
        idx_draw = sub2ind(size(imageData{selectedImage}), ...
            xind(save_points(:,1))', yind(save_points(:,2))', ...
            repmat(curr_slice, [size(save_points,1),1]));
        
    elseif shape==1 || shape==2
        
       % Erase previous shape (must be redrawn here for undo/redo)
        if ~isempty(idx_draw)
            imageData{selectedImage}(idx_draw) = prev_shape_colors;
            if selectedImage > 1
                alphaData{selectedImage}(idx_draw) = prev_alpha;
            end
        end
        
        % Clear in-motion saved shape:
        prev_shape_colors = [];
        prev_alpha = [];
      
    end

    % Add to undoCache:
    if selectedImage==1
        cacheState(selectedImage, 'draw', [], ...
            idx_draw, imageData{selectedImage}(idx_draw), [])
    else
        cacheState(selectedImage, 'draw', [], ...
            idx_draw, imageData{selectedImage}(idx_draw), ...
            alphaData{selectedImage}(idx_draw))
    end

    % Cache drawing for propagation
    curr_drawing = zeros(size(imageData{selectedImage})); % used for propagating last drawing
    curr_drawing(idx_draw) = draw_color;

    % Add drawing
    imageData{selectedImage}(idx_draw) = draw_color;
    if selectedImage > 1
        if erase_draw
            alphaData{selectedImage}(idx_draw) = NaN;
        else
            alphaData{selectedImage}(idx_draw) = alphaValue{selectedImage};
        end
    end
    
    idx_draw = [];
    in_motion = false; 
    updateImage
end

function scroll_zoom_callback(~, eventdata, ~)
    
    % Get mouse location
    pt = get(handles.axes{1},'CurrentPoint');
    pt = round(pt(1,1:2));
    % Check if within Image
    if pt(2) <= xdim && pt(2) > 0 && pt(1) <= ydim && pt(1) > 0
        scroll_count = min(max(scroll_count + -eventdata.VerticalScrollCount,0),9);
        if scroll_count ~= 0
            zoom_factor = scroll_zoom_equiv(scroll_count+1);
            xlength = zoom_factor*xdim; ylength = zoom_factor*ydim;
            xmax = min(round(xind(pt(2))+.5*xlength),xdim);
            xmin = max(round(xind(pt(2))-.5*xlength),1);
            ymax = min(round(yind(pt(1))+.5*ylength),ydim);
            ymin = max(round(yind(pt(1))-.5*ylength),1);
            xind = xmax:-1:xmin; 
            if slice_orientation==2
                yind = ymax:-1:ymin;
            else
                yind = ymin:ymax;
            end
        else % if scroll_count == 0, reset to full image view
            xmin = 1; xmax = xdim;
            ymin = 1; ymax = ydim; 
            xind = xdim:-1:1; 
            if slice_orientation==2
                yind = ydim:-1:1;
            else
                yind = 1:ydim;
            end
        end
        
        ax_xlim = [1,numel(yind)]; ax_ylim = [1,numel(xind)];
        refresh_img = true; 
        updateImage
        repositionAxes
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overlay Callbacks:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function closeOverlay(~, ~, closeWhich) 
    
    if nargin < 3 || isempty(closeWhich)
        closeWhich = selectedImage;
    end
    
    if closeWhich > 1
        
        % Clear any overlay edits from undo/redo cache
        purgeUndoRedo(closeWhich)
        
        if length(handles.images) >= closeWhich && ~isempty(handles.images{closeWhich}) && ishandle(handles.images{closeWhich})
            delete(handles.images{closeWhich})
            handles.images{closeWhich} = [];
        end
        if length(handles.axes) >= closeWhich && ~isempty(handles.axes{closeWhich}) && isgraphics(handles.axes{closeWhich})
            delete(handles.axes{closeWhich})
            handles.axes{closeWhich} = [];
        end
        
        % Shift left:
        try
            delete(h_image(closeWhich))
            delete(h_close(closeWhich))
        catch
        end
        
        imageData{closeWhich} = []; 
        alphaData{closeWhich} = []; 
        fullPaths{closeWhich} = [];
        
        cmax{closeWhich} = []; 
        cmin{closeWhich} = [];  
       
        % Change selectedImage to nearest non-empty
        if selectedImage==closeWhich
            selectedImage = find(~isempty(imageData));
            selectedImage = max(selectedImage);
            changeSelection([], [], selectedImage)
        end
        
        updateImage
    else
        warning('No overlays are currently loaded.')
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Drawing Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Draw Color Spec Callback:
function change_color_callback(hObject, ~, whichColor)
    draw_callback
    
    % Uncheck all:
    for k = 1:length(h_color); set(h_color(k), 'Checked', 'off'); end
    set(h_color(whichColor + 1), 'Checked', 'on') 
    
    if whichColor == 0
        erase_draw = true;
    else
        erase_draw = false; 
    end
        
    % Adjust color:
    if whichColor == 10 % custom
        prompt = {'Image intensity value: '};
        dlg_title = 'Color'; num_lines = [1,25];
        defaultans = {num2str(draw_color)};
        answer1 = inputdlg(prompt, dlg_title, num_lines, defaultans, 'on');
        if ~isempty(answer1)
            draw_color = str2double(answer1);
            if isnan(draw_color)
              errordlg('Please enter a numeric value', ...
                  'Invalid Input', 'modal')
              uicontrol(hObject)
              return 
            end
        end
    else
        draw_color = whichColor;
    end
            
end

function draw_shapes_callback(~, ~, type)
    draw_callback
    
    % Uncheck all
    for k = 1:length(h_shapes); set(h_shapes(k), 'Checked', 'off'); end
    set(h_shapes(type + 1), 'Checked','on')
    
    % Set type
    shape = type;
    
    % If sphere, get radius
    if shape==3 % sphere
        prompt = {'Specify sphere radius in physical units of image:'};
        dlg_title = ''; num_lines = [1,25]; 
        answer = inputdlg(prompt,dlg_title,num_lines);
        if isempty(answer); return; end
        sphere_radius = str2double(answer(1));
    end
end

function propagate_draw_callback(~, ~, ~)
    draw_callback
    
    if ~isempty(curr_drawing)
        % Prompt for slices
        prompt = {sprintf('Slice Range: \n\nFirst:'),'Last:'};
        dlg_title = ''; num_lines = [1,15;1,15]; 
        answer = inputdlg(prompt,dlg_title,num_lines);
        if isempty(answer); return; end
        
        % Slices
        slice_first = str2double(answer(1)); 
        slice_last = str2double(answer(2));
        
        % Get indices
        [ix,iy,~] = ind2sub(size(curr_drawing), find(curr_drawing(:)>0));
        slices = (slice_first:slice_last)';
        iz = repmat(slices,1, numel(ix))';
        
        idx_propagate = sub2ind(size(imageData{selectedImage}), ...
            repmat(ix, [length(slices), 1]), ...
            repmat(iy, [length(slices), 1]), ...
            iz(:));
        
        % Add to undoCache
        if selectedImage==1
            cacheState(selectedImage, 'draw', [], idx_propagate, ...
                imageData{selectedImage}(idx_propagate), [])
        else
            cacheState(selectedImage, 'draw', [], idx_propagate, ...
                imageData{selectedImage}(idx_propagate), ...
                alphaData{selectedImage}(idx_propagate))
        end
        
        % Propagate
        curr_drawing(idx_propagate) = draw_color;
        imageData{selectedImage}(idx_propagate) = draw_color;
        if selectedImage > 1
            if erase_draw
                alphaData{selectedImage}(idx_propagate) = NaN;
            else
                alphaData{selectedImage}(idx_propagate) = alphaValue{selectedImage};
            end
        end
        
        updateImage
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% updateImage Functions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateImage
    
    % Update Main Image (Background):
    refresh_img = false;
    
    % Update Image
    for i = 1:length(imageData)

        if isempty(imageData{i}); continue; end
        
        slice = imageData{i}(xind, yind, curr_slice);
        
        % Calculate color limits
        if colorscaleType{i} == 1 % Slice min/max
            cmin{i} = min(slice(:)); 
            cmax{i} = max(slice(:)); 
            if cmin{i}==cmax{i}; cmax{i} = cmin{i} + 1; end
        elseif colorscaleType{i} == 2 % Global min/max 
            cmin{i} = min(imageData{i}(:)); 
            cmax{i} = max(imageData{i}(:)); 
            if cmin{i}==cmax{i}; cmax{i} = cmin{i} + 1; end
        end

        % Map slice data into color indices units:
        idx = min(colormap_n, round((colormap_n-1)*(slice-cmin{i})/(cmax{i}-cmin{i}))+1);

        % Check if image is present
        if numel(handles.images) >= i && (~isempty(handles.images{i}) && ishandle(handles.images{i}))

            % Update previously generated image
            if i == 1 % background
                handles.images{1}.CData = idx;
            else % overlay
                set(handles.images{i},'Parent',handles.axes{i},...
                    'CData',idx,... 
                    'AlphaData',alphaData{i}(xind,yind,curr_slice),...
                    'AlphaDataMapping','none');
            end
            
        else % image not already present

            if i == 1 % background
                
                refresh_img = true;
                if isempty(parsed_inputs.axes)
                    handles.images{1} = image(idx);
                    handles.axes{1} = gca;
                    colormap(handles.axes{1}, colormaps{1});
                    set(handles.axes{1},'Color',axes_color,...
                        'XColor',axes_color,'YColor',...
                        axes_color,'ZColor',axes_color,...
                        'GridColor',axes_color)
                else 
                    handles.images{1} = image('CData',idx,'Parent',handles.axes{1});
                    set(handles.axes{1},'YDir','reverse');
                    handles.images{1}.Parent = parsed_inputs.axes;
                end
                hold 'on';
                
            else % overlay

                % Create overlay image
                if ~isempty(handles.axes{i})
                    handles.images{i} = image(...
                        'Parent',handles.axes{i},...
                        'CData',idx,...
                        'AlphaData',alphaData{i}(xind,yind,curr_slice),...
                        'AlphaDataMapping','none'); 
                    refresh_img = 1;
                end
            
            end

        end
    end
    
    ax_xlim = [1,numel(yind)];
    ax_ylim = [1,numel(xind)];
    xticks = round(linspace(1, ax_ylim(2), n_ticks_x));
    yticks = round(linspace(1, ax_xlim(2), n_ticks_y));
%     disp(['x-range: [', num2str(min(xind)),', ', num2str(max(xind)),']'])
%     disp(['y-range: [', num2str(min(yind)),', ', num2str(max(yind)),']'])
    
    % Update Colorbar:
    updateColorbar
    
    % Update Title:
    updateTitle
    
    % Reposition Axes:
    if refresh_img
        repositionAxes
    end
    
end

function updateColorbar(~)
        
    % Add/Remove Colorbar(s)
    if colorbar_on
        
        % Remake colorbar, delete previous
        if nargin==1 && isgraphics(h_colorbar,'colorbar')
            delete(h_colorbar)
        end
        
        % Create or update colorbar 
        if ~isgraphics(h_colorbar,'colorbar') || nargin==1 % create
            h_colorbar = colorbar(handles.axes{selectedImage},...
                'Position',colorbar_pos,...
                'Color',axes_color);
        else % update pre-existing colorbar
            set(h_colorbar,'Visible','on','Color',axes_color)
        end
        
        % Set Colorbar Ticks
        h_colorbar.Limits = [0, colormap_n + 1]; % scaled indices
        h_colorbar.LimitsMode = 'manual'; 
        h_colorbar.Ticks = round(.05 * colormap_n) : colormap_n/n_colorbar_ticks : colormap_n;
        colorbar_cvec = linspace(cmin{selectedImage}, cmax{selectedImage}, colormap_n); % raw units
        if (cmax{selectedImage} - cmin{selectedImage}) >= n_colorbar_ticks % (.15 * colormap_n)
            h_colorbar.TickLabels = cellstr(sprintf('%1g\n', round(colorbar_cvec(round(h_colorbar.Ticks)))));                    
        else
            h_colorbar.TickLabels = cellstr(sprintf('%.2f\n',colorbar_cvec(round(h_colorbar.Ticks))));
        end    
        
    else % assure not visible
        
        if isgraphics(h_colorbar,'colorbar')
            set(h_colorbar,'Visible','off');
        end
        
    end
   
end

function updateTitle
    % Update Title:
    if title_on
        refresh_img = 1;
        if ~isgraphics(h_title)
            h_title = title(['Slice ',num2str(curr_slice)],...
                'FontName','Helvetica','FontSize',10,'FontWeight','Bold',...
                'Units','normalized','Visible','on','Color',axes_color); 
            title_pos = get(h_title,'Position');
            title_pos(2) = 1.003*title_pos(2);
            set(h_title,'Position',title_pos)
        else
            h_title.String = ['Slice ',num2str(curr_slice)];
            h_title.Visible = 'on';
        end
    elseif isgraphics(h_title)
        h_title.Visible = 'off';
    end
end

function repositionAxes(~)
    % Determine Positioning
    if title_on && colorbar_on && axis_tick_on % 1
        curr_axis_pos = ax_pos_all_on;
        prev_state = curr_state; curr_state = 1;
    elseif ~title_on && ~colorbar_on && ~axis_tick_on % 2
        curr_axis_pos = ax_pos_all_off;
        prev_state = curr_state; curr_state = 2;
    elseif title_on && ~colorbar_on && axis_tick_on % 3
        curr_axis_pos = ax_pos_no_colorbar;
        prev_state = curr_state; curr_state = 3;
    elseif title_on && colorbar_on && ~axis_tick_on % 4
        curr_axis_pos = ax_pos_no_tick;
        prev_state = curr_state; curr_state = 4;
    elseif ~title_on && colorbar_on && axis_tick_on % 5
        curr_axis_pos = ax_pos_no_title;
        prev_state = curr_state; curr_state = 5;
    elseif title_on && ~colorbar_on && ~axis_tick_on % 6
        curr_axis_pos = ax_pos_title_only;
        prev_state = curr_state; curr_state = 6;
    elseif ~title_on && colorbar_on && ~axis_tick_on % 7
        curr_axis_pos = ax_pos_colorbar_only;
        prev_state = curr_state; curr_state = 7;
    elseif ~title_on && ~colorbar_on && axis_tick_on % 8
        curr_axis_pos = ax_pos_tick_only;
        prev_state = curr_state; curr_state = 8;
    end
    if prev_state~=curr_state || nargin==1
        x_ax_percent = curr_axis_pos(3); 
        y_ax_percent = curr_axis_pos(4);
        resize1 = true;
    else
        resize1 = false;
    end
    colorbar_pos = [curr_axis_pos(1)+curr_axis_pos(3),...
        curr_axis_pos(2),colorbar_pos(3),curr_axis_pos(4)];
    % Check if axes were specified by user (in which case, don't reposition)
    if ~isempty(parsed_inputs.axes)
        curr_axis_pos = ax_pos_input; 
        resize1 = false;
    end
    % Apply Positioning:
    if strcmp(units,'physical')
        xticklabels = (yind(yticks)-origin(2)) .* voxSize(2);
        yticklabels = (xind(xticks)-origin(1)) .* voxSize(1);
        xticklabels = round(xticklabels,3,'significant');
        yticklabels = round(yticklabels,3,'significant');
    elseif strcmp(units,'voxel')
       xticklabels = yind(yticks);
       yticklabels = xind(xticks);
    end
    for i = 1:length(handles.axes)
        if ishandle(handles.axes{i})
            set(handles.axes{i},'XLim',ax_xlim,'YLim',ax_ylim, ...
                'XTick',yticks,'YTick',xticks, 'Position', curr_axis_pos, ...
                'XTickLabels', xticklabels, 'YTickLabels', yticklabels)
        end
    end
    if colorbar_on && ishandle(h_colorbar)
        set(h_colorbar,'Position',[colorbar_pos(1),...
            curr_axis_pos(2),colorbar_pos(3),curr_axis_pos(4)])
    end
    if resize1
        resizeFigure
    end
end
            
end % End Function