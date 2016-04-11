function Ig = my_gaussian_filter(I,d)
%function to apply lowpass filter to an image using a Gaussian kernel
% I -> the original image
% d -> standard deviation of the gaussian

% generate the gaussian filter
g = fspecial('gaussian',size(I,1),d);
g1 = mat2gray(g);

%calculate fourier transform of the image
If = fftshift(fft2(I));
Ifg = If.*g1;

%inverse fourier tranform for obtaining the filtered image
Ig = real(ifft2(fftshift(Ifg)));