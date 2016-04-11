%my function to make rotation and translation where 
%tx,ty,tz are the distances for translation
%Qx,Qy,Qz are the angles for rotations
function P =  my_3d_affine(I,tx,ty,tz,Qx,Qy,Qz)

%create the original grid
[Xgrid, Ygrid, Zgrid] = meshgrid(1:size(I,1),1:size(I,2),1:size(I,3));
pts = [Xgrid(:) Ygrid(:) Zgrid(:)] - size(I,1)/2;

%take sizes of the image for later use
dx = size(I,1);
dy = size(I,2);
dz = size(I,3);

%transform into radians
% xtheta = 90*Qx/180;
% ytheta = 90*Qy/180;
% ztheta = 90*Qz/180;

%make rotation
rx = [1 0 0 ; 0 cos(xtheta) -sin(xtheta) ; 0 sin(xtheta) cos(xtheta)];
ry = [cos(ytheta) 0 sin(ytheta) ; 0 1 0 ; -sin(ytheta) 0 cos(ytheta)];
rz = [cos(ztheta) -sin(ztheta) 0 ; sin(ztheta) cos(ztheta) 0 ; 0 0 1];
rotate = rz*ry*rx;

pts_rotate = rotate * pts'; 

%make translation
row = ones(1, size(pts,1));
pts2 = [pts_rotate ; row];
translate = [1 0 0 tx ; 0 1 0 ty ; 0 0 1 tz ; 0 0 0 1];
pts_translate = (translate*pts2)' + size(I,1)/2 ;

%create the final image
P = interp3(I, reshape(pts_translate(:,1),dx,dy,dz), reshape(pts_translate(:,2),dx,dy,dz), reshape(pts_translate(:,3),dx,dy,dz));
end