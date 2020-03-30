function upvector = roll_neuroimage(arg1,arg2)
% ROLL_NEUROIMAGE
% 
% Adapted from internal Matlab function camroll.m
% -Instead of rolling along axis campos -> camtarget, this allows rolling
% along the Y-axis

% Parse Inputs:
if nargin==1
  ax = gca;
  theta = arg1;
elseif nargin==2
  ax = arg1;
  theta = arg2;
else error(message('MATLAB:camroll:TooManyInputs'))
end

% cpSave  = get(ax, 'cameraposition' );
% ctSave  = get(ax, 'cameratarget'   );

xlimits = get(ax,'XLim');
ylimits = get(ax,'YLim');
zlimits = get(ax,'ZLim');
origin = [.5*diff(xlimits)+xlimits(1),.5*diff(ylimits)+ylimits(1),.5*diff(zlimits)+zlimits(1)];
cam_pos = [origin(1),origin(2)+.5*origin(2),origin(3)];

upvector  = get(ax,'cameraupvector');

if ~righthanded(ax), theta = -theta; end

v = (origin-cam_pos);
v = v/norm(v);
v = [0,0,1];
alph = (theta)*pi/180;
cosa = cos(alph);
sina = sin(alph);
vera = 1 - cosa;
x = v(1);
y = v(2);
z = v(3);
rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
       x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
       x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';

r = crossSimple(v, upvector);
u = crossSimple(r, v);

newUp = u*rot;
upvector = newUp/norm(newUp);
set(ax, 'cameraupvector', upvector);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c=crossSimple(a,b)
c(1) = b(3)*a(2) - b(2)*a(3);
c(2) = b(1)*a(3) - b(3)*a(1);
c(3) = b(2)*a(1) - b(1)*a(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val=righthanded(ax)

dirs=get(ax, {'xdir' 'ydir' 'zdir'}); 
num=length(find(lower(cat(2,dirs{:}))=='n'));

val = mod(num,2);
