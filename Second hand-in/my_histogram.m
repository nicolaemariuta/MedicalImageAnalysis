%my function to calculate histogram of an image with pixels intensities
%from 0 to 255
function H = my_histogram(I)

%first element in H is count for pixel intensities 0; 
%the last one is for pixel intensity 255
H = zeros(256);
for R = 1: size(I,1)
    for C = 1: size(I,2)
        H(I(R,C)+1) = H(I(R,C)+1) +1;
    end 
end 

end 