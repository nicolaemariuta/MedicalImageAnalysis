%Testing Load image
B = my_load_mgh('orig_1.mgz');

I = squeeze(B(100,:,:));

imagesc(I),colormap('gray');

%% Exercise 1: Magnitude and phase of image 
clear;
B = my_load_mgh('orig_1.mgz');

%read the desired slice of the image and calculate the Fourier Transform
I = squeeze(B(100,:,:));
Ifft = fft2(I);

%calculate magnitude and phase for the image in frequency domain
Mag = abs(Ifft);
Phase = angle(Ifft);

%reconstruct images 
ImMag = abs(ifft2(Mag));
ImPhase = abs((ifft2(Phase)));

%plot the images
subplot(2,2,1),imagesc(I),colormap('gray'),title('original image');
subplot(2,2,2),imagesc(log(1+ImMag)),colormap('gray'),title('magnitude reconstruction');
subplot(2,2,3),imagesc(I),colormap('gray'),title('original image');
subplot(2,2,4),imagesc(log(1+ImPhase)),colormap('gray'),title('phase reconstruction');




%% Exercise2: Low pass filter using Gaussian
clear;
%load the image
B = my_load_mgh('orig_1.mgz');
I = squeeze(B(100,:,:));

%apply low pass filter using different fitlers
Ig1 = my_gaussian_filter(I,0.1);
Ig2 = my_gaussian_filter(I,0.5);
Ig3 = my_gaussian_filter(I,1.0);
Ig4 = my_gaussian_filter(I,1.5);
Ig5 = my_gaussian_filter(I,10.0);

%plot the results
subplot(3,2,1);imagesc(I),colormap('gray'),title('Original image');
subplot(3,2,2);imagesc(Ig1),colormap('gray'),title('Filter with std = 0.1');
subplot(3,2,3);imagesc(Ig2),colormap('gray'),title('Filter with std = 0.5');
subplot(3,2,4);imagesc(Ig3),colormap('gray'),title('Filter with std = 1.0');
subplot(3,2,5);imagesc(Ig4),colormap('gray'),title('Filter with std = 1.5');
subplot(3,2,6);imagesc(Ig5),colormap('gray'),title('Filter with std = 10.0');

%% Exercise 3: Histogram Equalization
clear;
%load the image
B = my_load_mgh('orig_1.mgz');
I = squeeze(B(100,:,:));

%calculate histogram
px = my_histogram(I);
%calculate cumulative histogram
cdf = imChist(I);
%perform histogram equalization
Ieq = my_histogram_equalization(I);

%calculate histogram of the equalized image
pxeq = my_histogram(Ieq);
cdfeq = imChist(Ieq);


%displaying results
subplot(2,3,1),imagesc(I),colormap('gray'),title('Original image');
subplot(2,3,2),plot(1:256,px),title('histogram of the original image');
subplot(2,3,3),plot(1:256, cdf),title('cumulative histogram of the original image');
subplot(2,3,4),imagesc(Ieq),colormap('gray'),title('Image after histogram equalization');
subplot(2,3,5),plot(1:256,pxeq),title('histogram of the processed image');
subplot(2,3,6),plot(1:256, cdfeq),title('cumulative histogram of the processed image');

%% Exercise 4: Median and average filter
clear;
%load the image
B = my_load_mgh('orig_1.mgz');
I = squeeze(B(100,:,:));

%median filtering
Im3 = medfilt2(I, [3 3]);
Im5 = medfilt2(I, [5 5]);
Im7 = medfilt2(I, [7 7]);
Im11 = medfilt2(I, [11 11]);
Im19 = medfilt2(I, [19 19]);

%average filtering
h3 = ones(3,3)/9;Ia3 = imfilter(I,h3);
h5 = ones(5,5)/25;Ia5 = imfilter(I,h5);
h7 = ones(7,7)/49;Ia7 = imfilter(I,h7);
h11 = ones(11,11)/121;Ia11 = imfilter(I,h11);
h19 = ones(19,19)/361;Ia19 = imfilter(I,h19);

%displaying results
subplot(2,6,1),imagesc(I),colormap('gray'),title('Original image');
subplot(2,6,2),imagesc(Im3),colormap('gray'),title('3x3 median filter');
subplot(2,6,3),imagesc(Im5),colormap('gray'),title('5x5 median filter');
subplot(2,6,4),imagesc(Im7),colormap('gray'),title('7x7 median filter');
subplot(2,6,5),imagesc(Im11),colormap('gray'),title('11x11 median filter');
subplot(2,6,6),imagesc(Im19),colormap('gray'),title('19x19 median filter');

subplot(2,6,7),imagesc(I),colormap('gray'),title('Original image');
subplot(2,6,8),imagesc(Ia3),colormap('gray'),title('3x3 average filter');
subplot(2,6,9),imagesc(Ia5),colormap('gray'),title('5x5 average filter');
subplot(2,6,10),imagesc(Ia7),colormap('gray'),title('7x7 average filter');
subplot(2,6,11),imagesc(Ia11),colormap('gray'),title('11x11 average filter');
subplot(2,6,12),imagesc(Ia19),colormap('gray'),title('19x19 average filter');




%%

I = imread('Syn1sec170.tif');
subplot(1,2,1),imshow(I);

Ifft = fft2(I);
Ifftshift = fftshift(Ifft);
subplot(1,2,2),imagesc(abs(log(Ifftshift)));


%% Exercise 5: Piecewise Linear Interpolation
clear;
%load the image
B = my_load_mgh('orig_1.mgz');
I = squeeze(B(100,:,:));


%create the original grid
[Xgrid, Ygrid] = meshgrid(1:(size(I,1)-2),1:(size(I,2)-2));
Xgrid = Xgrid+1;
Ygrid = Ygrid+1;
pts = [Xgrid(:) Ygrid(:)];

%create the new grid for interpolate
XgridShift1 = Xgrid +0.2;
YgridShift1 = Ygrid +0.2;
ptsShift1 = [XgridShift1(:) YgridShift1(:)];

XgridShift2 = Xgrid +0.5;
YgridShift2 = Ygrid +0.5;
ptsShift2 = [XgridShift2(:) YgridShift2(:)];

XgridShift3 = Xgrid +0.8;
YgridShift3 = Ygrid +0.8;
ptsShift3 = [XgridShift3(:) YgridShift3(:)];


%apply interpolation
Interpolate1 = my_interpolation(I,pts,ptsShift1);
Interpolate2 = my_interpolation(I,pts,ptsShift2);
Interpolate3 = my_interpolation(I,pts,ptsShift3);

subplot(2,2,1),imagesc(I),colormap('gray'),title('Original image');
subplot(2,2,2),imagesc(Interpolate1),colormap('gray'),title('Interpolation after shift (0.2,0.2)');
subplot(2,2,3),imagesc(Interpolate1),colormap('gray'),title('Interpolation after shift (0.5,0.5)');
subplot(2,2,4),imagesc(Interpolate1),colormap('gray'),title('Interpolation after shift (0.8,0.8)');


%% Exercise 6: Gradient of an image by using Sobel operator
clear;
%load the image
B = my_load_mgh('orig_1.mgz');
I = squeeze(B(100,:,:));

%define the Sobel operator kernels
Gx = [-1 0 1 ; -2 0 2 ; -1 0 1];
Gy = [-1 -2 -1 ; 0 0 0 ; 1 2 1];

%fourier transform of the image and the kernels
Ifft = fft2(I);

Gxfft = fft2(Gx);
Gxfft2 = zeros(size(I));
xpos = 1; ypos = 1;
Gxfft2(xpos:xpos+3-1,ypos:ypos+3-1) = Gxfft;

Gyfft = fft2(Gy);
Gyfft2 = zeros(size(I));
xpos = 1; ypos = 1;
Gyfft2(xpos:xpos+3-1,ypos:ypos+3-1) = Gyfft;

%apply convolution in fourier domain 
Gradxfft = Gxfft2.*Ifft;
Gradyfft = Gyfft2.*Ifft;

%ifft for each gradient
Gradx = real(ifft2(Gradxfft));
Grady = real(ifft2(Gradyfft));

Igrad = sqrt(Gradx.^2 + Grady.^2);

subplot(1,2,1),imagesc(I),colormap('gray'),title('Original image');
subplot(1,2,2),imagesc(Igrad),colormap('gray'),title('Gradient');

%% Exercise 7: DBC implementation
clear;
%load the image
v1 = my_load_mgh('images/nu_1.mgz');
Iv1 = squeeze(v1(:,100,:));

v2 = my_load_mgh('images/nu_2.mgz');
Iv2 = squeeze(v2(:,100,:));

%log transform of the images
lv1 = log(Iv1);
lv2 = log(Iv2);

%calculate the differential bias field
diff1 = lv1-lv2;
diff2 = lv2-lv1;
med_diff1 = medfilt2(diff1, [11 11]);
med_diff2 = medfilt2(diff2, [11 11]);
ratio1 = exp(med_diff1);
ratio2 = exp(med_diff2);

%plot the results
subplot(2,3,1),imagesc(Iv1),colormap('gray'),title('v1');
subplot(2,3,4),imagesc(Iv2),colormap('gray'),title('v2');
subplot(2,3,2),imagesc(ratio1),colormap('gray'),title('b1/b2');
subplot(2,3,5),imagesc(ratio2),colormap('gray'),title('b2/b1');
subplot(2,3,3),imagesc(Iv1-ratio2),colormap('gray'),title('v1 - b1/b2');
subplot(2,3,6),imagesc(Iv2-ratio1),colormap('gray'),title('v2 - b2/b1');




%% Question 9.4.d.
clear;
%load the image
B = my_load_mgh('orig_1.mgz');
I = squeeze(B(100,:,:));

%window size
M = 7;
N = 7;
mid_val=round((M*N)/2);

%find the number of rows and columns to be padded with zero
in=0;
for i=1:M
    for j=1:N
        in=in+1;
        if(in==mid_val)
            PadM=i-1;
            PadN=j-1;
            break;
        end
    end
end

%padding the image with zero on all sides
B=padarray(I,[PadM,PadN]);


for i= 1:size(B,1)-((PadM*2)+1)
    
    for j=1:size(B,2)-((PadN*2)+1)
        cdf=zeros(256,1);
        inc=1;
        for x=1:M
            for y=1:N
  %find the middle element in the window    
                if(inc==mid_val)
                    ele=B(i+x-1,j+y-1)+1;
                end
                    pos=B(i+x-1,j+y-1)+1;
                    cdf(pos)=cdf(pos)+1;
                   inc=inc+1;
            end
        end
                      
        %compute the cdf for the values in the window
        for l=2:256
            cdf(l)=cdf(l)+cdf(l-1);
        end
            Img(i,j)=round(cdf(ele)/(M*N)*255);
     end
end


%calculate histogram
px = my_histogram(I);
%calculate cumulative histogram
cdf = imChist(I);

%calculate histogram of the equalized image
pxeq = my_histogram(B);
cdfeq = imChist(B);


%displaying results
subplot(2,3,1),imagesc(I),colormap('gray'),title('Original image');
subplot(2,3,2),plot(1:256,px),title('histogram of the original image');
subplot(2,3,3),plot(1:256, cdf),title('cumulative histogram of the original image');
subplot(2,3,4),imagesc(B),colormap('gray'),title('Image after local histogram equalization');
subplot(2,3,5),plot(1:256,pxeq),title('histogram of the processed image');
subplot(2,3,6),plot(1:256, cdfeq),title('cumulative histogram of the processed image');






