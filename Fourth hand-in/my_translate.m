%my function for translation, takes as input the original 2d image and the
%number of pixels to move along x and y 
function Itranslate =  my_translate(I,tx, ty) 

%translation matrix
A = [1 0 tx ; 0 1 ty ; 0 0 1];


%create the original grid
[Xgrid, Ygrid] = meshgrid(1:size(I,1),1:size(I,2));
pts = [Xgrid(:) Ygrid(:)];

%calculate the new coordinates of the points
column = ones(size(pts,1),1);
pts = [pts column];
pts2 = (A * pts')';

%create the final image
Itranslate = interp2(I,reshape(pts2(:,1),256,256),reshape(pts2(:,2),256,256));

end
