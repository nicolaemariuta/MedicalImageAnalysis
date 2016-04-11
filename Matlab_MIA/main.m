close all; clear; clc;
load('RigTransformed.mat');

% Im1 = load_mgh('nu_1.mgz');
% Im1 = imresize3d(Im1,0.5,[],'linear','fill');
% Im2 = load_mgh('nu_2.mgz');
% Im2 = imresize3d(Im2,0.5,[],'linear','fill');
% %%
% % base, warp, weight, iterations
% [~, B ] = RigReg(Im1, Im2, 20, 50);
% 
%save('RigTransformed.mat','B', 'Im1', 'Im2');

%%
figure;
subplot(2,3,1), imshow(uint8(squeeze(Im1(60,:,:))));
subplot(2,3,2), imshow(uint8(squeeze(Im2(60,:,:))));
subplot(2,3,3), imshow(uint8(squeeze(B(60,:,:))));
subplot(2,3,5), imshow(uint8(squeeze(Im1(60,:,:)-Im2(60,:,:))));
subplot(2,3,6), imshow(uint8(squeeze(Im1(60,:,:)-B(60,:,:))));

%%
baseSSD = SSD(Im1, Im2);
rigSSD = SSD(Im1, B);

%%
mex -O dDdPFunc.cpp;

mex -O SplineInterpolation.cpp;

mex -O PNorm_det.cpp;

%% Baseline grid
[x, y, z] = ndgrid(1:8:128,1:8:128,1:8:128);
pts = [x(:) y(:) z(:)];
% rTrivial basline values in pts grid
rTrivial = interp3(Im1, x, y, z,'linear',0);

clear x y z
%% FFD Grid
[x] = ndgrid(1:16:128,1:16:128,1:16:128);
p = zeros([size(x)+2 3]);

clear x
%%
% phi : delta x
% dfdp : delta phi/ delta p
% idx : controllpoint mapping
% p FFD N*N*N*3 image value
offset = [-16, -16, -16];
spacing = [16, 16, 16];
[~, dfdp, idx] = SplineInterpolation(pts, p, offset, spacing); 
szP = size(p);
%% Minimizing CF

options.Method = 'lbfgs';
options.MaxIter = 100;
options.MaxFunEvals = 300;

controlpts = minFunc(@nonRigCost, p(:), options, pts,rTrivial,B,dfdp,idx,szP,spacing);

% %% Adding FFD movment to image grid
% deltaPoints = SplineInterpolation(pts(:), controlpts, offset, spacing);
% 
% ptsMoved = pts(:) + deltaPoints(:);
% ptsMoved = reshape(ptsMoved, size(pts));
% %%
% save('values.mat', 'controlpts');
% clear;

%% Interpolating image from moved grid
controlpts=importdata('values.mat');
[x,y,z]=ndgrid(1:128,1:128,1:128);
pts = [x(:) y(:) z(:)];

spacing = [16 16 16];
offset = -spacing;

deltaPoints = SplineInterpolation(pts(:), controlpts, offset, spacing);

newPts = pts+deltaPoints;


newImage = interp3(B, reshape(newPts(:,1),size(FloatingImage)), reshape(newPts(:,2),size(FloatingImage)), reshape(newPts(:,3),size(FloatingImage)));



% 
% 
% %% Interpolating image from moved grid
% load('values.mat');
% offset = [0, 0, 0];
% spacing = [8, 8, 8];
% R = SplineInterpolation(B(:), ptsMoved, offset, spacing);
% 
