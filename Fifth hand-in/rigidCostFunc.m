function [f, df] = rigidCostFunc(p, I, points, center, rTrival, scale, similarity, allowScaling, constrainRotation)
  % function [f df] = rigidCostFunc(p, I, points, center, rTrival, scale, similarity, allowScaling)
  %
  % Cost function for the rigid transformation.
  %
  % INPUT:
  % - p           : Vector of (9x1) scalars defining the rigid transformation
  %                 (rotX, rotY, rotZ, tX, tY, tZ, s1, s2, s3).
  % - I           : The source series.
  % - points      : A (3xN) matrix of points for rigid transformation.
  % - center      : A [x,y,z] vector indicating the center of rotation.
  % - rTrival     : The intensities in the evaluation points (target series).
  % - scale       : The voxel size.
  % - similarity  : The similarity measure (string).
  % - allowScaling: A binary variable indicating whether scaling of the points is allows.
  % - constrainRotation: Whether excessive rotation is allowed or angles>1 are penalized heavily
  %
  % OUTPUT:
  % - f : Function value (cost).
  % - df: Derivatives.
  
  % Apply rigid transformation to the given points
  
  
  [pts, S, R] = applyRigid(p, points, center, allowScaling);
  
  % Calcuate the derivatives of the transformation
  [dRx, dRy, dRz, dSx, dSy, dSz] = rigidDerivatives(p);
  
  % Do not weight the points
  det = ones(size(pts,1),1);

  % Compute similarity

      % b-splines with support two yield the +2 to ensure support.
      % [0 130 0 130] is interpreted as 0=smallest allowed voxel intensity. 130 is not used to limit upper range
      % similarly and is in essence a dead parameter.
      % [140 140] is the number of bins in the 2d histogram that the function is guaranteed to be able to handle.
      % Images should be rescaled to a max of 128 (to allow for overhead) when 140 is used. If rescaled to anything
      % smaller only a subset of the histogram is populated but this is otherwise handled correctly.
      % [0 0 0] offset for b-spline interpolation - ignore.
      % scale gives voxel (an)isotrpphy - should not be other than [1 1 1] without thinking.
      % det gives weights of local estimates - ignore.
      [res d(:,1) d(:,2) d(:,3)] = PNorm_det(pts, rTrival+2, I+2, [0 130 0 130], [140 140], [0 0 0], [scale(1) scale(2) scale(3)],2,det);

  
  % Normalize derivatives
  d = d./repmat(scale,size(d,1),1);
  
  % Calculate derivites of similarity measure wrt. parameters
  df          = zeros(numel(p),1);
  sizePoints  = size(points,1);
  centeredPts = points - repmat(center,sizePoints,1);
  
  % Translation
  df(4) = sum(d(:,1));
  df(5) = sum(d(:,2));
  df(6) = sum(d(:,3));
  % Rotation
  df(1) = sum(sum( centeredPts*S*dRx.*d, 2));
  df(2) = sum(sum( centeredPts*S*dRy.*d, 2));
  df(3) = sum(sum( centeredPts*S*dRz.*d, 2));
  % Scaling
  if (~allowScaling)
    df(7) = 0;
    df(8) = 0;
    df(9) = 0;
  else
    df(7) = sum(sum( centeredPts*dSx*R.*d, 2));
    df(8) = sum(sum( centeredPts*dSy*R.*d, 2));
    df(9) = sum(sum( centeredPts*dSz*R.*d, 2));
  end
  
  switch similarity
    case 'ssd'
      f = res;
    otherwise
      error('Unknown similarity measure: %s', similarity)
  end
  df=df(:);
  
  
  if constrainRotation
    threshold = pi/4; % In radians
    for i = 1:3
      if abs(p(i))>threshold 
        pen = abs(p(i))-threshold;
        % There is an implicit factor of "1" in front of the penalty term. This matches the typical cost function
        % dynamic range between 5 and 4.5. So an added penalty around 1 is effectively "STOP!"
        f = f + pen^2;
        df(i) = df(i) + 2*pen;
      end
    end
  end
end
