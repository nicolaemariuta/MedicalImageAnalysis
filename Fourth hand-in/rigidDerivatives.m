function [dRx, dRy, dRz, dSx, dSy,dSz] = rigidDerivatives(p)
  % function [dRx, dRy, dRz, dSx, dSy,dSz] = rigidDerivatives(p)
  %
  % Calculates the partial derivatives of a rigid transformation.
  %
  % INPUT:
  % - p: A vector of 9 scalars defining the rigid transformation
  %      (rotX, rotY, rotZ, tX, tY, tZ, s1, s2, s3).
  %
  % OUTPUT:
  % - dRx: x-derivative of rotation matrix.
  % - dRy: y-derivative of rotation matrix.
  % - dRz: z-derivative of rotation matrix.
  % - dSx: x-derivative of scaling matrix.
  % - dSy: y-derivative of scaling matrix.
  % - dSz: z-derivative of scaling matrix.
  
  
  % Initialize
  rotX = p(1);
  rotY = p(2);
  rotZ = p(3);
  
  % Rotation matrices
  [rX rY rZ] = rotationMatrices(rotX, rotY, rotZ);
  
  % Derived rotation matrices
  dRx = rZ*rY*[0 0 0; 0 -sin(rotX) -cos(rotX); 0 cos(rotX) -sin(rotX)];
  dRy = rZ*[-sin(rotY) 0 cos(rotY); 0 0 0; -cos(rotY) 0 -sin(rotY)]*rX;
  dRz = [-sin(rotZ) -cos(rotZ) 0; cos(rotZ) -sin(rotZ) 0; 0 0 0]*rY*rX;
  
  % Derived scaling matrices
  dSx = [1 0 0;0 0 0; 0 0 0];
  dSy = [0 0 0; 0 1 0; 0 0 0];
  dSz = [0 0 0; 0 0 0; 0 0 1];
end