function J = my_jacobian(dx,dy,dz)

%calculate gradients
[gx_x,gx_y,gx_z] = gradient(dx);
[gy_x,gy_y,gy_z] = gradient(dy);
[gz_x,gz_y,gz_z] = gradient(dz);

%determinant calculation
J = gx_x.*gy_y.*gz_z + ...
    gy_x.*gz_y.*gx_z + ...
    gz_x.*gx_y.*gy_z - ...
    gz_x.*gy_y.*gx_z - ...
    gx_x.*gx_y.*gx_z - ...
    gy_x.*gy_y.*gy_z;

%add identity matrix
I3 = zeros(size(dx));
I3(1:1+size(dx,1)+size(dx,2)*size(dx,3):end) = 1;
J = I3 + J;

end