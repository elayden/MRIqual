function [handles] = neuroimage_select_color(cmap)

    handles = struct('figure',[],'axes',[],'image',[]);
    handles.figure = figure('menu','none','Name','Choose Color','NumberTitle','off');
    handles.figure.Position(4) = 70; 
    handles.axes = gca; set(handles.axes,'Position',[0,0,1,1],'Visible','off','nextplot','add');

    if nargin<1 || isempty(cmap)
        cmap = jet(600);
    end

    c_image = zeros(100,size(cmap,1),3);
    c_image(:,:,1) = repmat(cmap(:,1)',100,1);
    c_image(:,:,2) = repmat(cmap(:,2)',100,1);
    c_image(:,:,3) = repmat(cmap(:,3)',100,1);
    handles.image = image(c_image,'Parent',handles.axes,'ButtonDownFcn',{@select_color,cmap});
    set(handles.axes,'XLim',[0,size(cmap,1)],'YLim',[0,100])

    function select_color(~,~,~)
        pt = get(gca,'CurrentPoint');
        xpos = round(pt(1,1));
        use_color = cmap(xpos,:);
        disp(use_color)
    end

end