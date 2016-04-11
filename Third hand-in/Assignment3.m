%%  rotate 2d image
clearvars;
%load image and set paramters
B = my_load_mgh('orig_1.mgz');
I = squeeze(B(:,:,100));
theta = 90*pi/180;

%use my function for rotate
Irotate = my_rotate(I,10);

%show the results
subplot(1,2,1),imagesc(I),colormap('gray'),title('original image');
subplot(1,2,2),imagesc(Irotate),colormap('gray'),title('rotate image');

%% translate 2d image 

%load image and set paramters
B = my_load_mgh('orig_1.mgz');
I = squeeze(B(:,:,100));
tx = 20;
ty = 140;

%use my function for translate
Itranslate = my_translate(I,tx,ty);

%show the results
subplot(1,2,1),imagesc(I),colormap('gray'),title('original image');
subplot(1,2,2),imagesc(Itranslate),colormap('gray'),title('translate image');

%% Complete affine transformation including translation and rotation
clearvars;

B = my_load_mgh('orig_1.mgz');
%trim the image so I can obtain a matrix that I can process
Bpart = B(97:160,97:160,97:160);


Btransformed = my_3d_affine (Bpart, 1, 1, 1, 10, 10, 10);
Btransformed(isnan(Btransformed)) = 0;

%display a slice to observe the changes
I = squeeze(Bpart(:,32,:));
Iaffine = squeeze(Btransformed(:,32,:));

%show the results
subplot(1,2,1),imagesc(I),colormap('gray'),title('original image');
subplot(1,2,2),imagesc(Iaffine),colormap('gray'),title('translate image');


%% Mutual information
clearvars;
%load the image
B = my_load_mgh('orig_1.mgz');
I = squeeze(B(:,:,100));

%rotate image multiple times and calculate the rigid transformation
mutual_array = squeeze(zeros(360,1));
for i = 1 : 360
    I2 = my_rotate(I,i);
    mutual_array(i) = my_mutual_information(I,I2);
end 

%plot the result
plot(1:360,mutual_array);



%% sum of square difference
clearvars;

B = my_load_mgh('orig_1.mgz');
I = squeeze(B(:,:,100));

%translate image multiple times and calculate the rigid transformation
ssd_array = squeeze(zeros(200,1));
for i = 1 : 200
    I2 = my_translate(I,i,i);
    ssd_array(i) = my_ssd(I,I2);
end 

%plot the result
plot(1:200,ssd_array);


%% registration
%load the image
clearvars;
B = my_load_mgh('orig_1.mgz');
%trim the image so I can obtain a matrix that I can process
B1 = B(97:160,97:160,97:160);
B2 = my_3d_affine (B1, 15, 25, 30, 1, 10, 20);

%df = my_cost_function(0,0,0,0,0,0,B1,B2);


P = [0 0 0 0 0 0];
fdf = @(P)my_cost_function(P,B1,B2);
odf = fminunc(fdf,P);

disp(odf);



