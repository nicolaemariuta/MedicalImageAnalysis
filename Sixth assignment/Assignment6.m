% kmeans clusttering k = 2
clear;
%load the images and rescale so I can compute
B = my_load_mgh('nu1.mgz');

I = imresize3d(B,1,[],'linear','bound');
I(isnan(I)) = 0;

%calculate histogram
pixels = floor(double(I(:)))+1;
histogram = accumarray(pixels,1);
data = [(1:256)'/256  histogram/max(histogram)];

%get initial random centroids
nrCentroids = 2;
centroids = [];
for i = 1:nrCentroids
    r = floor(rand*256)+1;
    centroids = [centroids ; data(r,1)  data(r,2)];
end 

%initialize the clusters
cluster1 = [];
cluster2 = [];

for i = 1:100
    %fill the clusters
    for j = 1 : size(data,1)
        point = squeeze(data(j,:));
        d1 = sqrt((centroids(1,1) - point(1))^2 + (centroids(1,2) - point(2))^2 );
        d2 = sqrt((centroids(2,1) - point(1))^2 + (centroids(2,2) - point(2))^2 );
        
        if(d1<d2)
            cluster1 = [cluster1 ; point];
        else 
            cluster2 = [cluster2 ; point];
        end 
     end 
    
    %calculate new means
    %calcualte mean for each cluster and find next centroids
    m1 = mean(cluster1);
    m2 = mean(cluster2);
    
    distances1 = sqrt(sum((cluster1 - ones(size(cluster1))*diag(m1)),2).^ 2);
    distances2 = sqrt(sum((cluster2 - ones(size(cluster2))*diag(m2)),2).^ 2);
    
    c1 = find(distances1 == min(distances1(:)));
    c2 = find(distances2 == min(distances2(:)));
    
    centroids(1,:) = cluster1(c1(1),:);
    centroids(2,:) = cluster2(c2(1),:);
    
    disp(i);
    
end

%find actual values of pixels in each cluster
[sharedVals1,idxs1] = intersect(data(:,1),cluster1(:,1));
[sharedVals2,idxs2] = intersect(data(:,1),cluster2(:,1));

%segment the image
segment1 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs1)));
segment1(lind) = I(lind);

segment2 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs2)));
segment2(lind) = I(lind);;

subplot(1,3,1),imagesc((squeeze(I(:,:,128)))),colormap('gray'),title('slice from baseline image');
subplot(1,3,2),imagesc((squeeze(segment1(:,:,128)))),colormap('gray'),title('slice from segment1');
subplot(1,3,3),imagesc((squeeze(segment2(:,:,128)))),colormap('gray'),title('slice from segment2');


%%
% kmeans clusttering k = 5
clear;
%load the images and rescale so I can compute
B = my_load_mgh('nu1.mgz');

I = imresize3d(B,1,[],'linear','bound');
I(isnan(I)) = 0;

%calculate histogram
pixels = floor(double(I(:)))+1;
histogram = accumarray(pixels,1);
data = [(1:256)'/256  histogram/max(histogram)];

%get initial random centroids
nrCentroids = 5;
centroids = [];
for i = 1:nrCentroids
    r = floor(rand*256)+1;
    centroids = [centroids ; data(r,1)  data(r,2)];
end 

%initialize the clusters
cluster1 = [];
cluster2 = [];
cluster3 = [];
cluster4 = [];
cluster5 = [];

for i = 1:100
    %fill the clusters
    for j = 1 : size(data,1)
        point = squeeze(data(j,:));
        d = [];
        d = [d sqrt((centroids(1,1) - point(1))^2 + (centroids(1,2) - point(2))^2 )];
        d = [d sqrt((centroids(2,1) - point(1))^2 + (centroids(2,2) - point(2))^2 )];
        d = [d sqrt((centroids(3,1) - point(1))^2 + (centroids(3,2) - point(2))^2 )];
        d = [d sqrt((centroids(4,1) - point(1))^2 + (centroids(4,2) - point(2))^2 )];
        d = [d sqrt((centroids(5,1) - point(1))^2 + (centroids(5,2) - point(2))^2 )];
        [m,im] = min(d);
       
        if(im == 1)
            cluster1 = [cluster1 ; point];
        elseif(im == 2)
            cluster2 = [cluster2 ; point];
        elseif(im == 3)
            cluster3 = [cluster3 ; point]; 
        elseif(im == 4)
            cluster4 = [cluster4 ; point];
        else
            cluster5 = [cluster5 ; point];
        end 
     end 
    
    %calculate new means
    %calcualte mean for each cluster and find next centroids
    m1 = mean(cluster1);
    m2 = mean(cluster2);
    m3 = mean(cluster3);
    m4 = mean(cluster4);
    m5 = mean(cluster5);
    
    distances1 = sqrt(sum((cluster1 - ones(size(cluster1))*diag(m1)),2).^ 2);
    distances2 = sqrt(sum((cluster2 - ones(size(cluster2))*diag(m2)),2).^ 2);
    distances3 = sqrt(sum((cluster3 - ones(size(cluster3))*diag(m3)),2).^ 2);
    distances4 = sqrt(sum((cluster4 - ones(size(cluster4))*diag(m4)),2).^ 2);
    distances5 = sqrt(sum((cluster5 - ones(size(cluster5))*diag(m5)),2).^ 2);
    
    c1 = find(distances1 == min(distances1(:)));
    c2 = find(distances2 == min(distances2(:)));
    c3 = find(distances3 == min(distances3(:)));
    c4 = find(distances4 == min(distances4(:)));
    c5 = find(distances5 == min(distances5(:)));
    
    centroids(1,:) = cluster1(c1(1),:);
    centroids(2,:) = cluster2(c2(1),:);
    centroids(3,:) = cluster3(c3(1),:);
    centroids(4,:) = cluster4(c4(1),:);
    centroids(5,:) = cluster5(c5(1),:);
    
    disp(i);
    
end

%find actual values of pixels in each cluster
[sharedVals1,idxs1] = intersect(data(:,1),cluster1(:,1));
[sharedVals2,idxs2] = intersect(data(:,1),cluster2(:,1));
[sharedVals3,idxs3] = intersect(data(:,1),cluster3(:,1));
[sharedVals4,idxs4] = intersect(data(:,1),cluster4(:,1));
[sharedVals5,idxs5] = intersect(data(:,1),cluster5(:,1));

%segment the image
segment1 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs1)));
segment1(lind) = I(lind);

segment2 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs2)));
segment2(lind) = I(lind);

segment3 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs3)));
segment3(lind) = I(lind);

segment4 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs4)));
segment4(lind) = I(lind);

segment5 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs5)));
segment5(lind) = I(lind);

subplot(2,3,1),imagesc((squeeze(I(:,:,128)))),colormap('gray'),title('slice from baseline image');
subplot(2,3,2),imagesc((squeeze(segment1(:,:,128)))),colormap('gray'),title('slice from segment1');
subplot(2,3,3),imagesc((squeeze(segment2(:,:,128)))),colormap('gray'),title('slice from segment2');
subplot(2,3,4),imagesc((squeeze(segment3(:,:,128)))),colormap('gray'),title('slice from segment3');
subplot(2,3,5),imagesc((squeeze(segment4(:,:,128)))),colormap('gray'),title('slice from segment4');
subplot(2,3,6),imagesc((squeeze(segment5(:,:,128)))),colormap('gray'),title('slice from segment5');

%%
% kmeans clusttering k = 9
clear;
%load the images and rescale so I can compute
B = my_load_mgh('nu1.mgz');

I = imresize3d(B,1,[],'linear','bound');
I(isnan(I)) = 0;

%calculate histogram
pixels = floor(double(I(:)))+1;
histogram = accumarray(pixels,1);
data = [(1:256)'/256  histogram/max(histogram)];

%get initial random centroids
nrCentroids = 9;
centroids = [];
for i = 1:nrCentroids
    r = floor(rand*256)+1;
    centroids = [centroids ; data(r,1)  data(r,2)];
end 

%initialize the clusters
cluster1 = [];
cluster2 = [];
cluster3 = [];
cluster4 = [];
cluster5 = [];
cluster6 = [];
cluster7 = [];
cluster8 = [];
cluster9 = [];

for i = 1:100
    %fill the clusters
    for j = 1 : size(data,1)
        point = squeeze(data(j,:));
        d = [];
        d = [d sqrt((centroids(1,1) - point(1))^2 + (centroids(1,2) - point(2))^2 )];
        d = [d sqrt((centroids(2,1) - point(1))^2 + (centroids(2,2) - point(2))^2 )];
        d = [d sqrt((centroids(3,1) - point(1))^2 + (centroids(3,2) - point(2))^2 )];
        d = [d sqrt((centroids(4,1) - point(1))^2 + (centroids(4,2) - point(2))^2 )];
        d = [d sqrt((centroids(5,1) - point(1))^2 + (centroids(5,2) - point(2))^2 )];
        d = [d sqrt((centroids(6,1) - point(1))^2 + (centroids(6,2) - point(2))^2 )];
        d = [d sqrt((centroids(7,1) - point(1))^2 + (centroids(7,2) - point(2))^2 )];
        d = [d sqrt((centroids(8,1) - point(1))^2 + (centroids(8,2) - point(2))^2 )];
        d = [d sqrt((centroids(9,1) - point(1))^2 + (centroids(9,2) - point(2))^2 )];
        [m,im] = min(d);
       
        if(im == 1)
            cluster1 = [cluster1 ; point];
        elseif(im == 2)
            cluster2 = [cluster2 ; point];
        elseif(im == 3)
            cluster3 = [cluster3 ; point]; 
        elseif(im == 4)
            cluster4 = [cluster4 ; point];
        elseif(im == 5)
            cluster5 = [cluster5 ; point];
        elseif(im == 6)
            cluster6 = [cluster6 ; point];
        elseif(im == 7)
            cluster7 = [cluster7 ; point];    
        elseif(im == 8)
            cluster8 = [cluster8 ; point];  
        else
            cluster9 = [cluster9 ; point];
        end 
     end 
    
    %calculate new means
    %calcualte mean for each cluster and find next centroids
    m1 = mean(cluster1);
    m2 = mean(cluster2);
    m3 = mean(cluster3);
    m4 = mean(cluster4);
    m5 = mean(cluster5);
    m6 = mean(cluster6);
    m7 = mean(cluster7);
    m8 = mean(cluster8);
    m9 = mean(cluster9);
    
    distances1 = sqrt(sum((cluster1 - ones(size(cluster1))*diag(m1)),2).^ 2);
    distances2 = sqrt(sum((cluster2 - ones(size(cluster2))*diag(m2)),2).^ 2);
    distances3 = sqrt(sum((cluster3 - ones(size(cluster3))*diag(m3)),2).^ 2);
    distances4 = sqrt(sum((cluster4 - ones(size(cluster4))*diag(m4)),2).^ 2);
    distances5 = sqrt(sum((cluster5 - ones(size(cluster5))*diag(m5)),2).^ 2);
    distances6 = sqrt(sum((cluster6 - ones(size(cluster6))*diag(m6)),2).^ 2);
    distances7 = sqrt(sum((cluster7 - ones(size(cluster7))*diag(m7)),2).^ 2);
    distances8 = sqrt(sum((cluster8 - ones(size(cluster8))*diag(m8)),2).^ 2);
    distances9 = sqrt(sum((cluster9 - ones(size(cluster9))*diag(m9)),2).^ 2);
    
    c1 = find(distances1 == min(distances1(:)));
    c2 = find(distances2 == min(distances2(:)));
    c3 = find(distances3 == min(distances3(:)));
    c4 = find(distances4 == min(distances4(:)));
    c5 = find(distances5 == min(distances5(:)));
    c6 = find(distances6 == min(distances6(:)));
    c7 = find(distances7 == min(distances7(:)));
    c8 = find(distances8 == min(distances8(:)));
    c9 = find(distances9 == min(distances9(:)));
    
    centroids(1,:) = cluster1(c1(1),:);
    centroids(2,:) = cluster2(c2(1),:);
    centroids(3,:) = cluster3(c3(1),:);
    centroids(4,:) = cluster4(c4(1),:);
    centroids(5,:) = cluster5(c5(1),:);
    centroids(6,:) = cluster6(c6(1),:);
    centroids(7,:) = cluster7(c7(1),:);
    centroids(8,:) = cluster8(c8(1),:);
    centroids(9,:) = cluster9(c9(1),:);
    
    disp(i);
    
end

%find actual values of pixels in each cluster
[sharedVals1,idxs1] = intersect(data(:,1),cluster1(:,1));
[sharedVals2,idxs2] = intersect(data(:,1),cluster2(:,1));
[sharedVals3,idxs3] = intersect(data(:,1),cluster3(:,1));
[sharedVals4,idxs4] = intersect(data(:,1),cluster4(:,1));
[sharedVals5,idxs5] = intersect(data(:,1),cluster5(:,1));
[sharedVals6,idxs6] = intersect(data(:,1),cluster6(:,1));
[sharedVals7,idxs7] = intersect(data(:,1),cluster7(:,1));
[sharedVals8,idxs8] = intersect(data(:,1),cluster8(:,1));
[sharedVals9,idxs9] = intersect(data(:,1),cluster9(:,1));

%segment the image
segment1 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs1)));
segment1(lind) = I(lind);

segment2 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs2)));
segment2(lind) = I(lind);

segment3 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs3)));
segment3(lind) = I(lind);

segment4 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs4)));
segment4(lind) = I(lind);

segment5 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs5)));
segment5(lind) = I(lind);

segment6 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs6)));
segment6(lind) = I(lind);

segment7 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs7)));
segment7(lind) = I(lind);

segment8 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs8)));
segment8(lind) = I(lind);

segment9 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs9)));
segment9(lind) = I(lind);

%
subplot(2,5,1),imagesc((squeeze(I(:,:,128)))),colormap('gray'),title('slice from baseline image');
subplot(2,5,2),imagesc((squeeze(segment1(:,:,128)))),colormap('gray'),title('slice from segment1');
subplot(2,5,3),imagesc((squeeze(segment2(:,:,128)))),colormap('gray'),title('slice from segment2');
subplot(2,5,4),imagesc((squeeze(segment3(:,:,128)))),colormap('gray'),title('slice from segment3');
subplot(2,5,5),imagesc((squeeze(segment4(:,:,128)))),colormap('gray'),title('slice from segment4');
subplot(2,5,6),imagesc((squeeze(segment5(:,:,128)))),colormap('gray'),title('slice from segment5');
subplot(2,5,7),imagesc((squeeze(segment6(:,:,128)))),colormap('gray'),title('slice from segment6');
subplot(2,5,8),imagesc((squeeze(segment7(:,:,128)))),colormap('gray'),title('slice from segment7');
subplot(2,5,9),imagesc((squeeze(segment8(:,:,128)))),colormap('gray'),title('slice from segment8');
subplot(2,5,10),imagesc((squeeze(segment9(:,:,128)))),colormap('gray'),title('slice from segment9');

%% Otsu's optimal value threshold for 2 segments
clear;
% load the images and rescale so I can compute
B = my_load_mgh('nu1.mgz');

I = imresize3d(B,1,[],'linear','bound');
I(isnan(I)) = 0;

% calculate normalized histogram
pixels = floor(double(I(:)))+1;
histogram = accumarray(pixels,1);
data = [(0:255)'  histogram/numel(I)];

%compute the cumulative sums
cumSums = zeros(size(data));
for i = 0:255
    cumSums(i+1,:) = [i sum(data(1:(i+1),2))];
end

%compute the cumulative means
cumMeans = zeros(size(data));
for i = 0:255
    cumMeans(i+1,:) = [i sum(data(1:(i+1),2).*(0:i)')];
end

%compute global intensity mean
mG = sum(data(:,2).*(0:255)');

%compute between-class variance
sigmab = ((mG*cumSums(:,2) - cumMeans(:,2)).^2)./(cumSums(:,2).*(1-cumSums(:,2)));
sigmab(isnan(sigmab)) = 0;
%optimal threshold
[ksigmab, kOptimal] = max(sigmab);

%apply segmentation
segment1 = zeros(size(I));
lind = sub2ind(size(I),find(I<kOptimal));
segment1(lind) = I(lind); 


segment2 = zeros(size(I));
lind = sub2ind(size(I),find(I>=kOptimal));
segment2(lind) = I(lind); 

%evaluate separability measure
sigmag = sum(data(:,2).*(((0:255) - mG)').^2);
measure = ksigmab/sigmag;



%display results
subplot(1,3,1),imagesc((squeeze(I(:,:,128)))),colormap('gray'),title('slice from baseline image');
subplot(1,3,2),imagesc((squeeze(segment1(:,:,128)))),colormap('gray'),title('slice from segment1');
subplot(1,3,3),imagesc((squeeze(segment2(:,:,128)))),colormap('gray'),title('slice from segment2');

%% Otsu's optimal value threshold for 5 segments
clear;
% load the images and rescale so I can compute
B = my_load_mgh('nu1.mgz');

I = imresize3d(B,1,[],'linear','bound');
I(isnan(I)) = 0;

% calculate normalized histogram
pixels = floor(double(I(:)))+1;
histogram = accumarray(pixels,1);
data = [(0:255)'  histogram/numel(I)];

%compute the cumulative sums
cumSums = zeros(size(data));
for i = 0:255
    cumSums(i+1,:) = [i sum(data(1:(i+1),2))];
end

%compute the cumulative means
cumMeans = zeros(size(data));
for i = 0:255
    cumMeans(i+1,:) = [i sum(data(1:(i+1),2).*(0:i)')];
end

%compute global intensity mean
mG = sum(data(:,2).*(0:255)');

%compute between-class variance
sigmab = ((mG*cumSums(:,2) - cumMeans(:,2)).^2)./(cumSums(:,2).*(1-cumSums(:,2)));
sigmab(isnan(sigmab)) = 0;
%optimal threshold
[ksigmab, kOptimal] = sort(sigmab);

kOptimal = sort(kOptimal(252:256,:));
%apply segmentation
segment1 = zeros(size(I));
lind = sub2ind(size(I),find(I<kOptimal(1)));
segment1(lind) = I(lind); 


segment2 = zeros(size(I));
lind = sub2ind(size(I),find(I>=kOptimal(1) & I<kOptimal(2)));
segment2(lind) = I(lind); 

segment3 = zeros(size(I));
lind = sub2ind(size(I),find(I>=kOptimal(2) & I<kOptimal(3)));
segment3(lind) = I(lind); 

segment4 = zeros(size(I));
lind = sub2ind(size(I),find(I>=kOptimal(4) & I<kOptimal(5)));
segment4(lind) = I(lind); 

segment5 = zeros(size(I));
lind = sub2ind(size(I),find(I>=kOptimal(5)));
segment5(lind) = I(lind); 

%display results
subplot(2,3,1),imagesc((squeeze(I(:,:,128)))),colormap('gray'),title('slice from baseline image');
subplot(2,3,2),imagesc((squeeze(segment1(:,:,128)))),colormap('gray'),title('slice from segment1');
subplot(2,3,3),imagesc((squeeze(segment2(:,:,128)))),colormap('gray'),title('slice from segment2');
subplot(2,3,4),imagesc((squeeze(segment3(:,:,128)))),colormap('gray'),title('slice from segment3');
subplot(2,3,5),imagesc((squeeze(segment4(:,:,128)))),colormap('gray'),title('slice from segment4');
subplot(2,3,6),imagesc((squeeze(segment5(:,:,128)))),colormap('gray'),title('slice from segment5');

%% Region growing
clear;
%load the images and rescale so I can compute
B = my_load_mgh('nu1.mgz');

I = imresize3d(B,1,[],'linear','bound');
I(isnan(I)) = 0;

%calculate histogram
pixels = floor(double(I(:)))+1;
histogram = accumarray(pixels,1);
data = [(1:256)'/256  histogram/max(histogram)];

%get initial random centroids
nrCentroids = 5;
centroids = [];
for i = 1:nrCentroids
    r = floor(rand*256)+1;
    centroids = [centroids ; data(r,1)  data(r,2)];
end 

%initialize the clusters
cluster1 = [];
cluster2 = [];
cluster3 = [];
cluster4 = [];
cluster5 = [];

for i = 1:100
    %fill the clusters
    for j = 1 : size(data,1)
        point = squeeze(data(j,:));
        d = [];
        d = [d sqrt((centroids(1,1) - point(1))^2 + (centroids(1,2) - point(2))^2 )];
        d = [d sqrt((centroids(2,1) - point(1))^2 + (centroids(2,2) - point(2))^2 )];
        d = [d sqrt((centroids(3,1) - point(1))^2 + (centroids(3,2) - point(2))^2 )];
        d = [d sqrt((centroids(4,1) - point(1))^2 + (centroids(4,2) - point(2))^2 )];
        d = [d sqrt((centroids(5,1) - point(1))^2 + (centroids(5,2) - point(2))^2 )];
        [m,im] = min(d);
       
        if(im == 1)
            cluster1 = [cluster1 ; point];
        elseif(im == 2)
            cluster2 = [cluster2 ; point];
        elseif(im == 3)
            cluster3 = [cluster3 ; point]; 
        elseif(im == 4)
            cluster4 = [cluster4 ; point];
        else
            cluster5 = [cluster5 ; point];
        end 
     end 
    
    %calculate new means
    %calcualte mean for each cluster and find next centroids
    m1 = mean(cluster1);
    m2 = mean(cluster2);
    m3 = mean(cluster3);
    m4 = mean(cluster4);
    m5 = mean(cluster5);
    
    distances1 = sqrt(sum((cluster1 - ones(size(cluster1))*diag(m1)),2).^ 2);
    distances2 = sqrt(sum((cluster2 - ones(size(cluster2))*diag(m2)),2).^ 2);
    distances3 = sqrt(sum((cluster3 - ones(size(cluster3))*diag(m3)),2).^ 2);
    distances4 = sqrt(sum((cluster4 - ones(size(cluster4))*diag(m4)),2).^ 2);
    distances5 = sqrt(sum((cluster5 - ones(size(cluster5))*diag(m5)),2).^ 2);
    
    c1 = find(distances1 == min(distances1(:)));
    c2 = find(distances2 == min(distances2(:)));
    c3 = find(distances3 == min(distances3(:)));
    c4 = find(distances4 == min(distances4(:)));
    c5 = find(distances5 == min(distances5(:)));
    
    centroids(1,:) = cluster1(c1(1),:);
    centroids(2,:) = cluster2(c2(1),:);
    centroids(3,:) = cluster3(c3(1),:);
    centroids(4,:) = cluster4(c4(1),:);
    centroids(5,:) = cluster5(c5(1),:);
    
    disp(i);
    
end

%find actual values of pixels in each cluster
[sharedVals1,idxs1] = intersect(data(:,1),cluster1(:,1));
[sharedVals2,idxs2] = intersect(data(:,1),cluster2(:,1));
[sharedVals3,idxs3] = intersect(data(:,1),cluster3(:,1));
[sharedVals4,idxs4] = intersect(data(:,1),cluster4(:,1));
[sharedVals5,idxs5] = intersect(data(:,1),cluster5(:,1));

%segment the image
segment1 = zeros(size(I));
lind1 = sub2ind(size(I),find(ismember(I,idxs1)));
segment1(lind1) = I(lind1);

segment2 = zeros(size(I));
lind2 = sub2ind(size(I),find(ismember(I,idxs2)));
segment2(lind2) = I(lind2);

segment3 = zeros(size(I));
lind3 = sub2ind(size(I),find(ismember(I,idxs3)));
segment3(lind3) = I(lind3);

segment4 = zeros(size(I));
lind4 = sub2ind(size(I),find(ismember(I,idxs4)));
segment4(lind4) = I(lind4);

segment5 = zeros(size(I));
lind5 = sub2ind(size(I),find(ismember(I,idxs5)));
segment5(lind5) = I(lind5);

%calculating the seeds by taking the mean of each segmentation
c1 = find(ismember(data,centroids(1,1)));
[x1,y1,z1] = ind2sub(size(I),find(ismember(I,c1)));
seed1 =[x1(100),y1(100),z1(100)];

c2 = find(ismember(data,centroids(2,1)));
[x2,y2,z2] = ind2sub(size(I),find(ismember(I,c2)));
seed2 = [x2(100),y2(100),z2(100)];

c3 = find(ismember(data,centroids(3,1)));
[x3,y3,z3] = ind2sub(size(I),find(ismember(I,c3)));
seed3 = [x3(100),y3(100),z3(100)];

c4 = find(ismember(data,centroids(4,1)));
[x4,y4,z4] = ind2sub(size(I),find(ismember(I,c4)));
seed4 = [x4(100),y4(100),z4(100)];

c5 = find(ismember(data,centroids(5,1)));
[x5,y5,z5] = ind2sub(size(I),find(ismember(I,c5)));
seed5 = [x5(100),y5(100),z5(100)];

[~,m1] = regionGrowing(squeeze(I),seed1,20, Inf, [], true, false);
[~,m2] = regionGrowing(squeeze(I),seed2,20, Inf, [], true, false);
[~,m3] = regionGrowing(squeeze(I),seed3,20, Inf, [], true, false);
[~,m4] = regionGrowing(squeeze(I),seed4,20, Inf, [], true, false);
[~,m5] = regionGrowing(squeeze(I),seed5,20, Inf, [], true, false);

%apply obtained masks to see the segmented regions
I1 = I.*m1;
I2 = I.*m2;
I3 = I.*m3;
I4 = I.*m4;
I5 = I.*m5;

%display results
subplot(2,3,1),imagesc((squeeze(I(:,:,128)))),colormap('gray'),title('slice from baseline image');
subplot(2,3,2),imagesc((squeeze(I1(:,:,128)))),colormap('gray'),title('slice from segment1');
subplot(2,3,3),imagesc((squeeze(I2(:,:,128)))),colormap('gray'),title('slice from segment2');
subplot(2,3,4),imagesc((squeeze(I3(:,:,128)))),colormap('gray'),title('slice from segment3');
subplot(2,3,5),imagesc((squeeze(I4(:,:,128)))),colormap('gray'),title('slice from segment4');
subplot(2,3,6),imagesc((squeeze(I5(:,:,128)))),colormap('gray'),title('slice from segment5');



%% Maximum likelihood segmentation

clear;
%load the images and rescale so I can compute
B = my_load_mgh('nu1.mgz');

I = imresize3d(B,1,[],'linear','bound');
I(isnan(I)) = 0;

%calculate histogram
pixels = floor(double(I(:)))+1;
histogram = accumarray(pixels,1);
data = [(1:256)'/256  histogram/max(histogram)];

%get initial random centroids
nrCentroids = 2;
centroids = [];
for i = 1:nrCentroids
    r = floor(rand*256)+1;
    centroids = [centroids ; data(r,1)  data(r,2)];
end 

%initialize the clusters
cluster1 = [];
cluster2 = [];

for i = 1:100
    %fill the clusters
    for j = 1 : size(data,1)
        point = squeeze(data(j,:));
        d1 = sqrt((centroids(1,1) - point(1))^2 + (centroids(1,2) - point(2))^2 );
        d2 = sqrt((centroids(2,1) - point(1))^2 + (centroids(2,2) - point(2))^2 );
        
        if(d1<d2)
            cluster1 = [cluster1 ; point];
        else 
            cluster2 = [cluster2 ; point];
        end 
     end 
    
    %calculate new means
    %calcualte mean for each cluster and find next centroids
    m1 = mean(cluster1);
    m2 = mean(cluster2);
    
    distances1 = sqrt(sum((cluster1 - ones(size(cluster1))*diag(m1)),2).^ 2);
    distances2 = sqrt(sum((cluster2 - ones(size(cluster2))*diag(m2)),2).^ 2);
    
    c1 = find(distances1 == min(distances1(:)));
    c2 = find(distances2 == min(distances2(:)));
    
    centroids(1,:) = cluster1(c1(1),:);
    centroids(2,:) = cluster2(c2(1),:);
    
    disp(i);
    
end

%find actual values of pixels in each cluster
[sharedVals1,idxs1] = intersect(data(:,1),cluster1(:,1));
[sharedVals2,idxs2] = intersect(data(:,1),cluster2(:,1));

%segment the image
segment1 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs1)));
segment1(lind) = I(lind);

segment2 = zeros(size(I));
lind = sub2ind(size(I),find(ismember(I,idxs2)));
segment2(lind) = I(lind);;

%apply maximum likelihood on the resulting regions

%number of elements in each region
n1 = numel(segment1);
n2 = numel(segment2);

histogram1 = accumarray(floor(double(n1(:)))+1,1);
data1 = [(0:255)'  histogram1/numel(n1)];

histogram2 = accumarray(floor(double(n2(:)))+1,1);
data2 = [(0:255)'  histogram2/numel(n2)];



