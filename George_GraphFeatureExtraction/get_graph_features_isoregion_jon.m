%% this function consider the global graph feature in the case where region of interest are isolated, e.g., epithelium regions are serveral isolated components in the binary mask

function [vfeature,Description] = get_graph_features_isoregion_jon(x,y,mask)

% graphfeats    Calculates graph-based features for nuclear centroids
% located at (x,y) in the image. 
% 
% Necessary input:
% x,y: x and y coordinates of points that will be used for graph
% construction (i.e. nuclear centroids). should be column vector
% mask is a binary map indicates the region of interest

% Output Description: vfeature contains the following:

% Voronoi Features
% 1: Area Standard Deviation
% 2: Area Average
% 3: Area Minimum / Maximum
% 4: Area Disorder
% 5: Perimeter Standard Deviation
% 6: Perimeter Average
% 7: Perimeter Minimum / Maximum
% 8: Perimeter Disorder
% 9: Chord Standard Deviation
% 10: Chord Average
% 11: Chord Minimum / Maximum
% 12: Chord Disorder

% Delaunay Triangulation
% 13: Side Length Minimum / Maximum
% 14: Side Length Standard Deviation
% 15: Side Length Average
% 16: Side Length Disorder
% 17: Triangle Area Minimum / Maximum
% 18: Triangle Area Standard Deviation
% 19: Triangle Area Average
% 20: Triangle Area Disorder

 load('GraphFeatureDescription.mat')
for i=1:20
    Description{i}=[GraphFeatureDescription{i} '-iso']; % append iso to the name
end
x_old=x;y_old=y;
%% need to split the whole graph/data points set into parts
 %% Cheng's old stuff
 
% cc= bwconncomp(mask);
% stats = regionprops(cc, 'Area');
% idx = find([stats.Area] > 100);
% mask = ismember(labelmatrix(cc), idx); %show( mask);
% %%%% project the mask to original dimension
% 
% all_idx=sub2ind([size(mask,1) size(mask,2)],round(y),round(x));
% setBasicGraph=[]; % store graph vetex
% 
% cc= bwconncomp(mask);

%%

% figure;
% axes('units','normalized','position',[0 0 1 1]);
% imagesc(I); axis ij; axis off; hold on; axis image;
% set(gcf,'Color',[1 1 1]);
%    plot(VX,VY,'g-','linewidth',2)
% hold on
dataPts=[x y];
centroid_list=[round(x),round(y)];


%% Actually, lets try doing an imdilate to identify clusters
SE = strel('disk',5,0);
groupmask=imdilate(mask,SE);
imshow(groupmask);
title('dilated mask');
groupmask=imerode(groupmask,SE);
imshow(groupmask);
title('re-eroded mask');
stats = regionprops(groupmask,'Centroid');
cc= bwconncomp(groupmask);



%% New nuclei identification program - assigns cells into neighborhoods based on who is close to them
% ncell.threshold=50;
% ncell.total=[centroid_list,zeros(size(centroid_list,1))];
% ncell.remaining=[centroid_list,1:size(centroid_list,1)];
% ncell.neighborhood_number=zeros(size(centroid_list,1));
% % ncell.neighborhood_number=1;
% % m=1;
% while ~isempty(ncell.remaining(1))
% %     for k1=2:size(ncell.remaining,1)
%         ncell.run=ncell.remaining(1,:);
%         ncell=for_cells_get_neighbors(ncell);
%         if max
%         ncell.total(ncell.remaining(1,3),3)=m;
%         ncell.remaining=ncell.remaining(2:end,:);
%         ncell=for_cells_get_neighbors(ncell);
%         ncell.total(ncell.remaining(ncell.in_neighborhood,3))=m; %takes the cells identified as neighbors, and sets the ncell total to their group number "m"
%         ncell.remaining=ncell.remaining(~ncell.in_neighborhood);
% 
% %     end
% end
%% Old setup

% for k=1:cc.NumObjects  %For each object in the mask - trying to find the index in cc that corresponds to the xy coodinates in dataPts
%     masktemp=zeros(size(mask,1),size(mask,2));
%     masktemp(cc.PixelIdxList{k})=1; %make mask of that specific object
% %     tempbounds=mask2bounds(masktemp);
% %     for k1=1:length(x)
% %         poly_list=inpoly(centroid_list, [tempbounds.c', tempbounds.r']);
% %     end
% 
% 
% dataPts_ID=0;
% for k1=1:size(dataPts,1)
%     if masktemp(centroid_list(k1,1),centroid_list(k1,2))==1
%         dataPts_ID=k1;
%         k1=size(dataPts,1);
% %         figure(1);
% %         imshow(masktemp);
% %         hold on;
% %         scatter(centroid_list(k,2), centroid_list(k,1));
% %         hold off;
%     end
% end
%     
% %     curmask_list=masktemp(:);% Linear vector of object pixels
% %     curnuclei_label=curmask_list(all_idx); 
% %     clustCent=dataPts(find(curnuclei_label),:);
% %     del = delaunay(clustCent(1,:), clustCent(2,:));
% %     triplot(del,clustCent(1,:), clustCent(2,:),'g','LineWidth',1);
% 
% if dataPts_ID==0;
% %     figure(2);
% %     imshow(masktemp);
% %     hold on;
% %     scatter(centroid_list(:,2), centroid_list(:,1));
% masktemp_centroid = regionprops(masktemp,'centroid');
%     for k1=1:size(dataPts,1)
%         d(k1) = pdist([centroid_list(k1,1),centroid_list(k1,2);masktemp_centroid.Centroid],'euclidean');  
%     end
%     [~,dataPts_ID]=min(d);
%     
% end
% 
%      setBasicGraph{k}=dataPts(dataPts_ID,:);
% end



%%
% for k=1:cc.NumObjects  %For each object in the mask - trying to find the index in cc that corresponds to the xy coodinates in dataPts
%     masktemp=zeros(size(mask,1),size(mask,2));
%     masktemp(cc.PixelIdxList{k})=1; %grab that specific object
%     curmask_list=masktemp(:);% Linear vector of object pixels
%     curnuclei_label=curmask_list(all_idx); % set the 
%     clustCent=dataPts(find(curnuclei_label),:);
% %     del = delaunay(clustCent(1,:), clustCent(2,:));
% %     triplot(del,clustCent(1,:), clustCent(2,:),'g','LineWidth',1);
%     setBasicGraph{k}=clustCent;


% end
% hold off;

%% Calculate the Voronoi diagram.
% [VX,VY] = voronoi(x,y);
[V, C] = voronoin([x(:),y(:)]);

% Okay, so:
% VX, VY - These guys contain the vertices in a way such that
% plot(VX,VY,'-',x,y,'.') creates the voronoi diagram. I don't actually
% think these are used later on, but I'm keeping them here just in case.

% C - This variable is an m by 1 cell array, where m is the number of cell
% centroids in your image. Each element in C is a vector
% with the coordinates of the vertices of that row's voronoi polygon. 
% V - This is a q by 2 matrix, where q is the number of vertices and 2 is
% the number of dimensions of your image. Each element in V contains the
% location of the vertex in 2D space.
% The idea here is that if you want to see the coordinates for the vertices
% of polygon 5, for example, you would go:
%     X = V(C{5},:)
% which would display five rows, each with the 2D coordinates of the vertex
% of polygon 5.

%% Get the delaunay triangulation...
% del = delaunay(x,y);    

% Returns a set of triangles such that no data points are contained in any 
% triangle's circumcircle. Each row of the numt-by-3 matrix TRI defines 
% one such triangle and contains indices into the vectors X and Y. When 
% the triangles cannot be computed (such as when the original data is 
% collinear, or X is empty), an empty matrix is returned.
    
%% Record indices of inf and extreme values to skip these cells later
Vnew        = V;
Vnew(1,:)   = [];

% Find the data points that lie far outside the range of the data
[Vsorted,I]     = sort([Vnew(:,1);Vnew(:,2)]);
N               = length(Vsorted);
Q1              = round(0.25*(N+1));
Q3              = round(0.75*(N+1));
IQR             = Q3 - Q1;
highrange       = Q3 + 1.5*IQR;
lowrange        = Q1 - 1.5*IQR;
Vextreme        = [];
Vextreme        = [Vextreme; V(find(V > highrange))];
Vextreme        = [Vextreme; V(find(V < lowrange))];

banned = [];
for i = 1:length(C)
    if(~isempty(C{i}))
        
    if(max(any(isinf(V(C{i},:)))) == 1 || max(max(ismember(V(C{i},:),Vextreme))) == 1)
        banned = [banned, i];
    end
    end
end
% If you've eliminated the whole thing (or most of it), then only ban 
% indices that are infinity (leave the outliers)
if(length(banned) > length(C)-2)
    banned = [];
    for i = 1:length(C)
        if(max(any(isinf(V(C{i},:)))) == 1)
            banned = [banned, i];
        end
    end
end

%% Voronoi Diagram Features
% Area
c = 1;
d = 1;
e = d;
for k=1:cc.NumObjects
    curBG=setBasicGraph{k};
    x=curBG(:,1);y=curBG(:,2);
    del=delaunay(x,y);
    for i = 1:length(C)
        if(~ismember(i,banned) && ~isempty(C{i}))
            X = V(C{i},:);
            chord(1,:) = X(:,1);
            chord(2,:) = X(:,2);
            % Calculate the chord lengths (each point to each other point)
            for ii = 1:size(chord,2)
                for jj = ii+1:size(chord,2)
                    chorddist(d) = sqrt((chord(1,ii) - chord(1,jj))^2 + (chord(2,ii) - chord(2,jj))^2);
                    d = d + 1;
                end
            end
            
            % Calculate perimeter distance (each point to each nearby point)
            for ii = 1:size(X,1)-1
                perimdist(e) = sqrt((X(ii,1) - X(ii+1,1))^2 + (X(ii,2) - X(ii+1,2))^2);
                e = e + 1;
            end
            perimdist(size(X,1)) = sqrt((X(size(X,1),1) - X(1,1))^2 + (X(size(X,1),2) - X(1,2))^2);
            
            % Calculate the area of the polygon
            area(c) = polyarea(X(:,1),X(:,2));
            c = c + 1;
            clear chord X
        end
    end
    if(~exist('area','var'))
        vfeature = zeros(1,20);
        return;
    end
end

vfeature(1) = std(area); 
vfeature(2) = mean(area);
vfeature(3) = min(area) / max(area);
vfeature(4) = 1 - ( 1 / (1 + (vfeature(1) / vfeature(2))) );

vfeature(5) = std(perimdist);
vfeature(6) = mean(perimdist);
vfeature(7) = min(perimdist) / max(perimdist);
vfeature(8) = 1 - ( 1 / (1 + (vfeature(5) / vfeature(6))) );

vfeature(9) = std(chorddist);
vfeature(10) = mean(chorddist);
vfeature(11) = min(chorddist) / max(chorddist);
vfeature(12) = 1 - ( 1 / (1 + (vfeature(9) / vfeature(10))) );

%% Delaunay 
% Edge length and area
c = 1;
d = 1;
for k=1:cc.NumObjects
    curBG=setBasicGraph{k};
    x=curBG(:,1);y=curBG(:,2);
    for i = 1:size(del,1)
        t = [x(del(i,:)),y(del(i,:))];
        
        sidelen(c:c+2) = [sqrt( ( t(1,1) - t(2,1) )^2 + (t(1,2) - t(2,2))^2 ), ...
            sqrt( ( t(1,1) - t(3,1) )^2 + (t(1,2) - t(3,2))^2 ), ...
            sqrt( ( t(2,1) - t(3,1) )^2 + (t(2,2) - t(3,2))^2 )];
        dis(i,1:3) = sum( sidelen(c:c+2) );
        c = c + 3;
        triarea(d) = polyarea(t(:,1),t(:,2));
        d = d + 1;
    end
end

vfeature(13) = min(sidelen) / max(sidelen);
vfeature(14) = std(sidelen); 
vfeature(15) = mean(sidelen);
vfeature(16) = 1 - (1 / (1 + (vfeature(14) / vfeature(15)) ) );

vfeature(17) = min(triarea) / max(triarea);
vfeature(18) = std(triarea);
vfeature(19) = mean(triarea);
vfeature(20) = 1 - (1 / (1 + (vfeature(18) / vfeature(19))) );

end

% function [ncell]=for_cells_get_neighbors(ncell)
%         ncell.in_neighborhood=zeros(size(ncell.run,1));
%         for k1=1:size(ncell.run,1)
%             distance=pdist(ncell.remaining(:,1:2));
%             distance=distance(1:size(ncell.remaining,1));
%             ncell.in_neighborhood=ncell.in_neighborhood || distance <= ncell.threshold;
%         end       
% end
