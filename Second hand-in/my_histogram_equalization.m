%my function for performing histogram equalization
function Ieq = my_histogram_equalization(I)

%calculate the cumulative histogram
cdf = imChist(I);

%Histogram equalization
Ieq = zeros(size(I));
for R = 1: size(I,1)
    for C = 1: size(I,2)
        Ieq(R,C) = floor(255*cdf(I(R,C)+1));
    end 
end 

end