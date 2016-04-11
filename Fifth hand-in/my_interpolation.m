%my function for interpolation. Input parameters: original image, original
%grid and the new points
function Interpolate = my_interpolation(I, pts, ptsShift)
% the interpolate image
Interpolate = zeros(size(I));

for i = 1:size(pts,1)
   xi = pts(i,1);yi = pts(i,2);
   x = ptsShift(i,1);y = ptsShift(i,2);
   
   %calculate weights for pixels in original image:
    w11 = ((xi+1) - x)*((yi+1) - y);
    w21 = (x - xi)*((yi+1) - y);
    w12 = ((xi+1) - x)*(y - yi);
    w22 = (x - xi)*(i - yi);
    
    %calculate interpolation at point (xi,yi) in the new grid
    Interpolate(xi,yi) = I(xi,yi)*w11 + I(xi+1,yi)*w21 +  I(xi,yi+1)*w12 + I(xi+1,yi+1)*w22/(((xi+1)-xi)*((yi+1)-yi));
end

end