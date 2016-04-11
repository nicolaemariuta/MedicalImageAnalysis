function [rx ry rz] = rotationMatrices(xtheta,ytheta,ztheta)

rx = [1 0 0 ; 0 cos(xtheta) -sin(xtheta) ; 0 sin(xtheta) cos(xtheta)];
ry = [cos(ytheta) 0 sin(ytheta) ; 0 1 0 ; -sin(ytheta) 0 cos(ytheta)];
rz = [cos(ztheta) -sin(ztheta) 0 ; sin(ztheta) cos(ztheta) 0 ; 0 0 1];

end