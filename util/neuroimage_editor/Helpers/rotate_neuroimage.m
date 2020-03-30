function rotate_neuroimage(ax,pitch,yaw,roll)

[azi,ele] = view(ax);
upvector  = camup(ax);

% Apply Pitch & Yaw:
if ~isempty(pitch) && isempty(yaw)
    view(ax,[azi,ele+pitch])
elseif isempty(pitch) && ~isempty(yaw)
    view(ax,[azi+yaw,ele])
elseif ~isempty(pitch) && ~isempty(yaw)
    view(ax,[azi+yaw,ele+pitch])
end
set(ax, 'cameraupvector', upvector);

% Apply Roll:
if ~isempty(roll) && roll~=0
    roll_neuroimage(ax,roll)
    view(ax,[azi,ele])
end


% Note: roll_neuroimage has right idea, need to replicate that for pitch
% and yaw

% xlimits = get(ax,'XLim');
% ylimits = get(ax,'YLim');
% zlimits = get(ax,'ZLim');
% origin = [.5*diff(xlimits)+xlimits(1),.5*diff(ylimits)+ylimits(1),.5*diff(zlimits)+zlimits(1)];

end