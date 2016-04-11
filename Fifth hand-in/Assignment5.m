clear;
%compile the code files
% mex -O 'SplineInterpolation.cpp'
% mex -O 'PNorm_det.cpp'
% mex -O 'dDdPFunc.cpp'
%load the images

%load the images and rescale so I can compute
B1 = my_load_mgh('nu1.mgz');
B1 = imresize3d(B1,0.25,[],'linear','bound');
B2 = my_load_mgh('nu2.mgz');
B2 = imresize3d(B2,0.25,[],'linear','bound');
I1 = B1;
I2 = B2;

% %apply rigid registration with Akshay cost function
options.Method = 'lbfgs';
options.MaxIter = 100;
P = [0 0 0 0 0 0 0 0 0];
[gridX,gridY,gridZ] = ndgrid(1:64,1:64,1:64);
points = [gridX(:) gridY(:) gridZ(:)];
rTrivial = interpn(gridX,gridY,gridZ,I1, gridX,gridY,gridZ,'cubic',0);
affinePoints = minFunc(@rigidCostFunc,P(:),options,I2,points,[32 32 32],rTrivial(:),[1 1 1],'ssd',0,0);


d1 = my_ssd(I1,I2);
regI2 = my_3d_affine(I2,affinePoints(4),affinePoints(5),affinePoints(6),affinePoints(1),affinePoints(2),affinePoints(3));
d2 = my_ssd(I1,regI2);
I2 = regI2;

%display the results
subplot(2,3,1),imagesc((squeeze(I1(:,:,32)))),colormap('gray'),title('slice from baseline image');
subplot(2,3,2),imagesc((squeeze(I2(:,:,32)))),colormap('gray'),title('slice from follow-up image');
subplot(2,3,3),imagesc((squeeze(regI2(:,:,32)))),colormap('gray'),title('slice from registered image');
subplot(2,3,4),imagesc(abs((squeeze(I1(:,:,32))) - (squeeze(I2(:,:,32))))),colormap('gray'),title('baseline - followup');
subplot(2,3,5),imagesc(abs((squeeze(I1(:,:,32))) - (squeeze(regI2(:,:,32))))),colormap('gray'),title('baseline - registered');

save('registeredI2.tiff','I2');

%%
clear;
%load the images and rescale so I can compute
B1 = my_load_mgh('nu1.mgz');
B1 = imresize3d(B1,0.25,[],'linear','bound');
tmp =load('registeredI2');
B2 = tmp.I2;
%B2 = imresize3d(B2,0.25,[],'linear','bound');
I1 = B1;
I2 = B2;

I1(isnan(I1)) = 0;
I2(isnan(I2)) = 0;


% %apply rigid registration with my cost function but it is not working
% P = [0 0 0 0 0 0];
% fdf = @(P)my_cost_function(P,I1,I2);
% odf = fminunc(fdf,P);
% I2 = my_3d_affine (I2, odf(4), odf(5),odf(6), odf(1), odf(2), odf(3));


%define spacing for the spline interpolation
spacing = [2 2 2];
offset = -spacing;
%create the evaluation grid
[gridX,gridY,gridZ] = ndgrid(1:2:64,1:2:64,1:2:64);
pts = [gridX(:) gridY(:) gridZ(:)];
%create the control points grid
[cpts] = ndgrid(1:4:64,1:4:64,1:4:64);
p = zeros([size(cpts)+2 3]);
%rTrivial values in the baseline grid
rTrivial = interpn(I1, gridX,gridY,gridZ,'cubic',0);
%make initial interpolation
[~,dfdp,idx] = SplineInterpolation(pts,p,offset,spacing);


szP = size(p);

%nonrigid cost function mostly for testing
d = my_ssd(I1,I2);
[f,df] = non_rigid_cost(p,pts,rTrivial, I2, dfdp, idx, szP, spacing);
disp('here');
%make optimization
options.Method = 'lbfgs';
options.MaxIter = 100;
controlPoints = minFunc(@non_rigid_cost,p(:),options,pts,rTrivial,I2,dfdp,idx,szP,spacing);


%use optimization control points to spline interpolate the follow-up image
[x y z] = ndgrid(1:64,1:64,1:64);
ptsfinal = [x(:) y(:) z(:)];
controlPoints = reshape(controlPoints,szP);
deltaPoints = SplineInterpolation(ptsfinal,controlPoints,offset,spacing);
newPts = ptsfinal + deltaPoints;
registeredImage = interpn(I2, reshape(newPts(:,3),size(I2)),reshape(newPts(:,2),size(I2)),reshape(newPts(:,1),size(I2)));

%calculate sum of square difference at the end
dfinal = my_ssd(I1,registeredImage);
% 
% %display the results
% subplot(2,3,1),imagesc((squeeze(I1(:,:,32)))),colormap('gray'),title('slice from baseline image');
% subplot(2,3,2),imagesc((squeeze(I2(:,:,32)))),colormap('gray'),title('slice from follow-up image');
% subplot(2,3,3),imagesc((squeeze(registeredImage(:,:,32)))),colormap('gray'),title('slice from registered image');
% subplot(2,3,4),imagesc(abs((squeeze(I1(:,:,32))) - (squeeze(I2(:,:,32))))),colormap('gray'),title('baseline - followup');
% subplot(2,3,5),imagesc(abs((squeeze(I1(:,:,32))) - (squeeze(registeredImage(:,:,32))))),colormap('gray'),title('baseline - registered');


%compute Jacobian
dx = deltaPoints(:,1);
dx = interpn(I2, reshape(dx,size(I2)),reshape(dx,size(I2)),reshape(dx,size(I2)));

dy = deltaPoints(:,2);
dy = interpn(I2, reshape(dy,size(I2)),reshape(dy,size(I2)),reshape(dy,size(I2)));

dz = deltaPoints(:,3);
dz = interpn(I2, reshape(dz,size(I2)),reshape(z,size(I2)),reshape(dz,size(I2)));

J = my_jacobian(dx,dy,dz);

%display the results
subplot(2,3,1),imagesc((squeeze(I1(:,:,32)))),colormap('gray'),title('slice from baseline image');
subplot(2,3,2),imagesc((squeeze(I2(:,:,32)))),colormap('gray'),title('slice from follow-up image');
subplot(2,3,3),imagesc((squeeze(registeredImage(:,:,32)))),colormap('gray'),title('slice from registered image');
subplot(2,3,4),imagesc(abs((squeeze(I1(:,:,32))) - (squeeze(I2(:,:,32))))),colormap('gray'),title('baseline - followup');
subplot(2,3,5),imagesc(abs((squeeze(I1(:,:,32))) - (squeeze(registeredImage(:,:,32))))),colormap('gray'),title('baseline - registered');
subplot(2,3,6),imagesc((squeeze(J(:,:,32)))),colormap('gray'),title('jacobian');



