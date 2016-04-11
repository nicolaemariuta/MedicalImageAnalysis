function [points, s, r, t] = applyRigid(p, points, center, allowScaling)
  % function [pts s r t] = applyRigid(p, points, center, allowScaling)
  %
  % Applies a rigid transformation to a set of points.
  %
  % INPUT:
  % - p           : Vector of (9x1) scalars defining the rigid transformation
  %                 (rotX, rotY, rotZ, tX, tY, tZ, s1, s2, s3).
  % - points      : A (3xN) matrix of points for rigid transformation.
  % - center      : A [x,y,z] vector indicating the center of rotation.
  % - allowScaling: A binary variable indicating whether scaling of the points is allows.
  %
  % OUTPUT:
  % - points: The rigidly transformed points.
  % - s  : Scaling matrix.
  % - r  : Rotation matrix.
  % - t  : Translation vector.
  
  
  % Initialize variables
  rotX = p(1);
  rotY = p(2);
  rotZ = p(3);
  
  t11 = p(4);
  t12 = p(5);
  t13 = p(6);
  
  if (allowScaling)
    s1 = p(7);
    s2 = p(8);
    s3 = p(9);
  else
    s1 = 0;
    s2 = 0;
    s3 = 0;
  end
  
  % Rotation matrices
  [rX rY rZ] = rotationMatrices(rotX, rotY, rotZ);
  r = rZ*rY*rX;
  
  % Scaling matrix
  s = diag([1+s1,1+s2,1+s3]);
  
  % Translation vector
  t = [t11; t12; t13];
  
  % Apply transformation
  dim = size(points,2);
  
  % The transformation below could be done in one-liner, but memory jumps quite a lot like this:
  % sizePoints = size(points,1);
  % pts = (points-repmat(center,sizePoints,1))*s*r + repmat(t',[sizePoints 1]) + repmat(center,sizePoints,1);
  
  % Translate to rotation center
  for d = 1:dim
    points(:,d) = points(:,d)-center(d);
  end

  % Scale and rotate
  points = points*s;
  points = points*r;
  
  % Translate back from center and then apply translation
  for d = 1:dim
    points(:,d) = points(:,d)+center(d)+t(d);
  end
end