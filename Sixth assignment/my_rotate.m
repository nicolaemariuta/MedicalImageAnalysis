%my function for interpolation, takes as input the original 2d image and angle
%in degrees
function Irotate =  my_rotate(I,pi) 

%create the original grid
[Xgrid, Ygrid] = meshgrid(1:size(I,1),1:size(I,2));
pts = [Xgrid(:) Ygrid(:)] - size(I,1)/2;

%transform into radians
theta = 90*pi/180;

%calculate the new coordinates of the points
pts2 = ([cos(theta) -sin(theta) ; sin(theta) cos(theta)]*pts')'+size(I,1)/2;

%use interpolation to obtain the finale image
Irotate = interp2(I,reshape(pts2(:,1),256,256),reshape(pts2(:,2),256,256));

end