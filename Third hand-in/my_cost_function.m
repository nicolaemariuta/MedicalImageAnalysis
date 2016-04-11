%my cost function; returns the ssd and derivative of ssd
function df =  my_cost_function(P, I1, I2) 
xtheta = P(1);
ytheta = P(2);
ztheta = P(3);
tx = P(4);
ty = P(5);
tz = P(6);

I1(isnan(I1)) = 0;
I2(isnan(I2)) = 0;

%calculate gradients of I2
[gx,gy,gz] = gradient(I2);
G = [gx(:) gy(:) gz(:)];

%affine transformation matrices:
rx = [1 0 0 0 ; 0 cos(xtheta) -sin(xtheta) 0 ; 0 -sin(xtheta) cos(xtheta) 0 ; 0 0 0 1];
ry = [cos(ytheta ) 0 sin(ytheta) 0 ; 0 1 0 0 ; -sin(ytheta) 0 cos(ytheta) 0 ; 0 0 0 1];
rz = [cos(ztheta) sin(ztheta) 0 0 ; -sin(ztheta) cos(ztheta) 0 0 ; 0 0 1 0 ; 0 0 0 1];
t = [1 0 0 tx ; 0 1 0 ty ; 0 0 1 tz ; 0 0 0 1];
Tx = [1 0 0 tx ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];
Ty = [1 0 0 0 ; 0 1 0 ty ; 0 0 1 0 ; 0 0 0 1];
Tz = [1 0 0 0 ; 0 1 0 0 ; 0 0 1 tz ; 0 0 0 1];
%gradients for transformation matrices:
drx = [0 0 0 0 ; 0 -sin(xtheta) -cos(xtheta) 0 ; 0 -cos(xtheta) -sin(xtheta) 0 ; 0 0 0 1];
dry = [-sin(ytheta) 0 cos(ytheta) 0 ; 0 1 0 0 ; -cos(ytheta) 0 -sin(ytheta) 0 ; 0 0 0 1];
drz = [-sin(ztheta) cos(ztheta) 0 0 ; -cos(ztheta) -sin(ztheta) 0 0 ; 0 0 1 0 ; 0 0 0 1];
dtx = [ 0 0 0 1 ; 0 0 0 0 ; 0 0 0 0 ; 0 0 0 0 ];
dty = [ 0 0 0 0 ; 0 0 0 1 ; 0 0 0 0 ; 0 0 0 0 ];
dtz = [ 0 0 0 0 ; 0 0 0 0 ; 0 0 0 1 ; 0 0 0 0 ];

%gradient of ssd over the 6 degrees of freedom
diff = I1-I2;
diff = diff(:);
d = diff*sum(G);


row = ones(1, size(d,1))';
disp(size(d));
disp(size(row));
d = [d  row];

dQrx = (1/numel(I2))*sum(sum((-2)*drx*ry*rz*t*d',2));
dQry = (1/numel(I2))*sum(sum((-2)*rx*dry*rz*t*d',2));
dQrz = (1/numel(I2))*sum(sum((-2)*rx*ry*drz*t*d',2));
dQtx = (1/numel(I2))*sum(sum((-2)*rx*ry*rz*dtx*d',2));
dQty = (1/numel(I2))*sum(sum((-2)*rx*ry*rz*dty*d',2));
dQtz = (1/numel(I2))*sum(sum((-2)*rx*ry*rz*dtz*d',2));

%matrix with all gradients
df = dQrx + dQry + dQrz + dQtx + dQty + dQtz;

end