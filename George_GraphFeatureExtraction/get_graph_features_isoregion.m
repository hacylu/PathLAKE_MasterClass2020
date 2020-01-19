%% this function consider the global graph feature in the case where region of interest are isolated, e.g., epithelium regions are serveral isolated components in the binary mask

function [vfeature,Description] = get_graph_features_isoregion(x,y,mask,para)

% graphfeats    Calculates graph-based features for nuclear centroids
% located at (x,y) in the image.
%
% Necessary input:
% x,y: x and y coordinates of points that will be used for graph
% construction (i.e. nuclear centroids). should be column vector
% mask is a binary map indicates the region of interest
% para.T_smallregion: threshold to remove small region


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
cc= bwconncomp(mask);
stats = regionprops(cc, 'Area');
idx = find([stats.Area] > para.T_smallregion);
mask = ismember(labelmatrix(cc), idx); %show( mask);
%%%% project the mask to original dimension

all_idx=sub2ind([size(mask,1) size(mask,2)],round(y),round(x));
setBasicGraph=[]; % store graph vetex

cc= bwconncomp(mask);

if para.show
    figure;
    axes('units','normalized','position',[0 0 1 1]);
    imagesc(para.I); axis ij; axis off; hold on; axis image;
    set(gcf,'Color',[1 1 1]);
        plot(x,y,'go','linewidth',2)
    hold on
end

dataPts=[x y];
for k=1:cc.NumObjects
    masktemp=zeros(size(mask,1),size(mask,2));
    masktemp(cc.PixelIdxList{k})=1; %show(masktemp);
    curmask_list=masktemp(:);% sum(masktemp)
    curnuclei_label=curmask_list(all_idx); % sum(curnuclei_label)
    clustCent=dataPts(logical(curnuclei_label),:);
    if size(clustCent,1)~=2
        clustCent=clustCent';
    end
    if para.show
        if ~isempty(clustCent)
            %%% delaunay
            del = delaunay(clustCent(1,:), clustCent(2,:));
            triplot(del,clustCent(1,:), clustCent(2,:),'g','LineWidth',1);
            %%% voronoi
%             figure;
%             show(para.I);hold on;
%             voronoi(clustCent(1,:), clustCent(2,:));
%             
        end
    end
    setBasicGraph{k}=clustCent;
end
if para.show
    hold off;
%     return;
end
%% Calculate the Voronoi diagram.
% [VX,VY] = voronoi(x,y);
% [V, C] = voronoin([x(:),y(:)]);

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
% Vnew        = V;
% Vnew(1,:)   = [];
% 
% % Find the data points that lie far outside the range of the data
% [Vsorted,I]     = sort([Vnew(:,1);Vnew(:,2)]);
% N               = length(Vsorted);
% Q1              = round(0.25*(N+1));
% Q3              = round(0.75*(N+1));
% IQR             = Q3 - Q1;
% highrange       = Q3 + 1.5*IQR;
% lowrange        = Q1 - 1.5*IQR;
% Vextreme        = [];
% Vextreme        = [Vextreme; V(find(V > highrange))];
% Vextreme        = [Vextreme; V(find(V < lowrange))];
% 
% banned = [];
% for i = 1:length(C)
%     if(~isempty(C{i}))
%         
%         if(max(any(isinf(V(C{i},:)))) == 1 || max(max(ismember(V(C{i},:),Vextreme))) == 1)
%             banned = [banned, i];
%         end
%     end
% end
% % If you've eliminated the whole thing (or most of it), then only ban
% % indices that are infinity (leave the outliers)
% if(length(banned) > length(C)-2)
%     banned = [];
%     for i = 1:length(C)
%         if(max(any(isinf(V(C{i},:)))) == 1)
%             banned = [banned, i];
%         end
%     end
% end

%% Voronoi Diagram Features
% Area
c = 1;% index that count for different measurements
d = 1;
e = d;
for k=1:cc.NumObjects
    curBG=setBasicGraph{k};
    x=curBG(1,:)';y=curBG(2,:)';
    if length(x)>3
        %% Calculate the Voronoi diagram.
        %         [VX,VY] = voronoi(x,y);
        [V, C] = voronoin([x(:),y(:)]);
        %% Get the delaunay triangulation...
        del = delaunay(x,y);
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
        
        %         del=delaunay(x,y);
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
                Varea(c) = polyarea(X(:,1),X(:,2));
                c = c + 1;
                clear chord X
            end
        end
    end
end
% if there is no feature can be computed, return 0
if(~exist('Varea','var'))
    vfeature = zeros(1,20);
    return;
end

vfeature(1) = std(Varea);
vfeature(2) = mean(Varea);
vfeature(3) = min(Varea) / max(Varea);
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
    x=curBG(1,:)';y=curBG(2,:)';
    if length(x)>3
        %% Calculate the Voronoi diagram.
        %         [VX,VY] = voronoi(x,y);
        [V, C] = voronoin([x(:),y(:)]);
        %% Get the delaunay triangulation...
        del = delaunay(x,y);
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
end

vfeature(13) = min(sidelen) / max(sidelen);
vfeature(14) = std(sidelen);
vfeature(15) = mean(sidelen);
vfeature(16) = 1 - (1 / (1 + (vfeature(14) / vfeature(15)) ) );

vfeature(17) = min(triarea) / max(triarea);
vfeature(18) = std(triarea);
vfeature(19) = mean(triarea);
vfeature(20) = 1 - (1 / (1 + (vfeature(18) / vfeature(19))) );
