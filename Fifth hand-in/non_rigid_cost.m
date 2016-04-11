function [f,df] = non_rigid_cost(p,pts,rTrivial, ImMove, dfdp, idx, szP, spacing)

%reshape the control points grid
p = reshape(p,szP);
Np = numel(p)/3;
offset = -spacing;

%apply spline interpolate and add delta to the evaluation points
deltaPoints = SplineInterpolation(pts, p, offset, spacing);
pts = pts + deltaPoints;

%dummy variables
det = ones(length(pts),1);
sourceVoxelSize = [1,1,1];
mx = 10;

%calculate the gradients
[f, d(:,1), d(:,2), d(:,3)] = PNorm_det(pts, rTrivial+2,ImMove+2,[0 130 0 130], [mx mx], [0 0 0], sourceVoxelSize, 2, double(det));
df = dDdPFunc(d,dfdp,idx,Np);

df = df(:);

end








