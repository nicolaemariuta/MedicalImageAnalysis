%my function to calculate cumulative histogram of an image
function C = imChist(I)

%calculate histogram
H = my_histogram(I);
C = cumsum(H);
count = size(I,1)*size(I,2);

for i = 1 : size(C)
    C(i) = C(i)/count;
end

end

