%function for calculating the sum of square difference
%just apply the formula
function d =  my_ssd(I1, I2) 
%replace all nans because otherwise some matlab functions will not work
I1(isnan(I1)) = 0;
I2(isnan(I2)) = 0;

%calculation of ssd
d =sum((I1(:) - I2(:)).^2)/numel(I1);

end
