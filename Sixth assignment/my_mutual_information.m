function jointEntropy = my_mutual_information(I1, I2)
%replace all nans because otherwise some matlab functions will not work
I1(isnan(I1)) = 0;
I2(isnan(I2)) = 0;
%calculate the joint histogram
indrow = floor(double(I1(:))) + 1;
indcol = floor(double(I2(:))) + 1; %// Should be the same size as indrow
jointHistogram = accumarray([indrow indcol], 1);
jointProb = jointHistogram / numel(indrow);

%joint entropy
indNoZero = jointHistogram ~= 0;
jointProb1DNoZero = jointProb(indNoZero);
jointEntropy = -sum(jointProb1DNoZero.*log2(jointProb1DNoZero));

end