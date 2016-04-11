%my cost function; returns the ssd and derivative of ssd
function [f, df] =  my_cost_function(xtheta,ytheta,ztheta,tx,ty,tz, I1, I2) 
I1(isnan(I1)) = 0;
I2(isnan(I2)) = 0;
%calculate ssd
f = my_ssd(I1,I2);

%calculate gradients of I2
[gx,gy,gz] = gradient(I2);

%affine transformation matrices:
rx = [1 0 0 ; 0 cos(xtheta) -sin(xtheta) ; 0 sin(xtheta) cos(xtheta)];
ry = [cos(ytheta) 0 sin(ytheta) ; 0 1 0 ; -sin(ytheta) 0 cos(ytheta)];
rz = [cos(ztheta) -sin(ztheta) 0 ; sin(ztheta) cos(ztheta) 0 ; 0 0 1];
t = [1 0 0 tx ; 0 1 0 ty ; 0 0 1 tz ; 0 0 0 1];
Tx = [1 0 0 tx ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];
Ty = [1 0 0 0 ; 0 1 0 ty ; 0 0 1 0 ; 0 0 0 1];
Tz = [1 0 0 0 ; 0 1 0 0 ; 0 0 1 tz ; 0 0 0 1];
%gradients for transformation matrices:
drx = [0 0 0 ; 0 -sin(xtheta) -cos(xtheta) ; 0 cos(xtheta) -sin(xtheta)];
dry = [-sin(ytheta) 0 cos(ytheta) ; 0 0 0 ; -cos(ytheta) 0 -sin(ytheta)];
drz = [-sin(ztheta) -cos(ztheta) 0 ; cos(ztheta) -sin(ztheta) 0 ; 0 0 0];
dtx = [ 0 0 0 1 ; 0 0 0 0 ; 0 0 0 0 ; 0 0 0 0 ];
dty = [ 0 0 0 0 ; 0 0 0 1 ; 0 0 0 0 ; 0 0 0 0 ];
dtz = [ 0 0 0 0 ; 0 0 0 0 ; 0 0 0 1 ; 0 0 0 0 ];


%gradient of ssd over the 6 degrees of freedom
dQrx = (1/numel(I2))*(-2)*sum(sum(sum(I1-I2)))*(gx+gy+gz)*(drx*ry*rz*t);
dQry = (1/numel(I2))*(-2)*sum(sum(sum(I1-I2)))*(gx+gy+gz)*(rx*dry*rz*t);
dQrz = (1/numel(I2))*(-2)*sum(sum(sum(I1-I2)))*(gx+gy+gz)*(rx*ry*drz*t);
dQtx = (1/numel(I2))*(-2)*sum(sum(sum(I1-I2)))*(gx+gy+gz)*(rx*ry*rz*dtx*Ty*Tz);
dQty = (1/numel(I2))*(-2)*sum(sum(sum(I1-I2)))*(gx+gy+gz)*(rx*ry*rz*Tx*dty*Tz);
dQtz = (1/numel(I2))*(-2)*sum(sum(sum(I1-I2)))*(gx+gy+gz)*(rx*ry*rz*Tx*Ty*dtz);


%%HINT
G = [gx(:) gy(:) gz(:)];

e = [I1-I2]; e = e(:);

d = -2*e*G;

dQrx = (1/numel(I2))*sum(sum((-2)*rz*ry*drx*t.*d,2));




%matrix with all gradients
df = [dQrx dQry dQrz dQtx dQty dQtz];




end