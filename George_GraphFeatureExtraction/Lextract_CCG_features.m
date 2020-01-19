% para.alpha_min= min parameter for the graph forming
% para.alpha_res= resolution between min-max parameter
% para.alpha_max= max parameter for the graph forming
% bounds.CellClusterC_c=column of vetex
% bounds.CellClusterC_r=row of vetex
function [feats,description,CCGinfo] = Lextract_CCG_features(bounds,para)
%% Extracts features from nuclei/glandular bounds
% bounds is a struct containing the centroid and boundary info
% The structure of bounds is as follows:
% {bounds.r,bounds.c,bounds.centroid_r,bounds.centroid_c]
% NOTE!!! the r = x, c=y coordinate

% img is a color or grayscale image used for computed for haralick features
% img can be omitted if not using haralick featuers

% if i really want to fix those names
% if strcmp(type, 'nuclei')
     load NewImgDescription.mat
%     temp = regexp(o.Description.ImageFeatures,'CGT');
%     temp2 = cellfun(@isempty,temp,'UniformOutput',false);
%
%     %description{3} = o.Description.ImageFeatures([temp2{:}] == 0);
% else
%     load NewImgDescription.mat
% end

% Dimred2 = @(x) Dimred(x, 'GE', param);
% EigV = cellfun(Dimred2,data,'UniformOutput',false);

fprintf('\nExtracting Cluster Graph Features...')
% note that this should be build on the cell clusters!!!
% bounds.CellCluster
CCbounds.centroid_c=bounds.CellClusterC_c;
CCbounds.centroid_r=bounds.CellClusterC_r;

CCGinfo=[];
set_alpha=[para.alpha_min:para.alpha_res:para.alpha_max];
for f=1:length(set_alpha)
    curpara.alpha=set_alpha(f); curpara.radius=para.radius;
    
    [CCGfeats,~,CCGinfo{f}]=extract_cluster_graph_feats(CCbounds,curpara);
    feats{f}=CCGfeats;

%     feats{f} = extract_cluster_graph_feats(CCbounds,curpara);
    temp = regexp(o.Description.ImageFeatures,'GSG');
    temp2 = cellfun(@isempty,temp,'UniformOutput',false);
    temp = [o.Description.ImageFeatures([temp2{:}] == 0)];
    for i=1:length(temp)
        cur=temp{i};
        cur(1:3)='CCG';
        str=sprintf('a=%.2f',curpara.alpha);
        cur(end+1:end+length(str))=str;
        temp{i}=cur;
    end
    description{f} =temp;
end





function [graphfeats] = extract_graph_feats(bounds)
gb_r = [bounds.centroid_r];
gb_c = [bounds.centroid_c];

[graphfeats] = get_graph_features([gb_r]',[gb_c]');


function [morphfeats]=extract_morph_feats(bounds)
%% Morph
% gb_r = {bounds.r};
% gb_c = {bounds.c};

badglands = [];
for j = 1:length(bounds.nuclei)
    try
        curB=bounds.nuclei{j};
%         [feat] = morph_features([gb_r{j}]',[gb_c{j}]');
                [feat] = morph_features(curB(:,1),curB(:,2));
        feats(j,:) = feat;
        
    catch ME
        badglands = [badglands j]; %track bad glands
    end
end

feats(badglands,:) = []; %remove bad glands

morphfeats = [mean(feats) std(feats) median(feats) min(feats)./max(feats)];



function [CGTfeats] = extract_CGT_feats(bounds,para)
%% CGT
% addpath(genpath(['Z:\2012_10_GlandOrientation']))

%[CGTfeats,c_matrix,info] = get_CGT_features(bounds);
a = para.alpha;
r = para.radius;

%[CGTfeats, c_matrix, info, feats, network, edges] = get_CGT_features_networks_weight(bounds,a,r);
[CGTfeats, c_matrix, info, feats, network, edges] = extract_CGT_features(bounds,a,r);

function [CCGfeats,feature_list,CCGinfo] = extract_cluster_graph_feats(bounds,para)
%CCG

%addpath(genpath(['Z:\2012_10_GlandOrientation']))

% build graph
alpha = para.alpha;
r = para.radius;
[VX,VY,x,y,edges] = construct_ccgs(bounds,alpha, r);
CCGinfo.VX=VX;
CCGinfo.VY=VY;
CCGinfo.x=x;
CCGinfo.y=y;
CCGinfo.edges=edges;
% [VX,VY,x,y,edges] = Lconstruct_ccgs_optimized(bounds,alpha, r);

%[CCGfeats,feature_list] = cluster_graph_features_networkboost(bounds, edges);
[CCGfeats,feature_list] = cluster_graph_features_optimized(bounds, edges);

function [Texturefeats] = extract_texture_feats(img)

%% texture
%addpath(genpath(['Z:\2012_10_GlandOrientation']))

% alternatively 13 haralick
%addpath(genpath('Z:\Datasets\SatishData\haralick'))

info.dist = 1;
info.win = 1;
info.grays = 256;

n = 0;

if ndims(img) < 3
    gimg = img;
elseif ndims(img) == 3
    gimg = rgb2gray(img); % assume img is rgb
else
    fprintf('Unidentified image format')
end

mask(:,:,1) = (gimg ~= max(max(gimg)));
%for grays = [64,128,256]
%    for win = [1,2,3]
%        n = n + 1
grays = info.grays;
win = info.win;
dist = info.dist;

himg = uint16(rescale_range(gimg,0,grays-1));
%himg = rescale_range(gimg,0,grays-1);
%f = haralick_img(himg,mask,grays,win,dist,1);
f = haralick_img(himg,mask,grays,win,dist,1);
Hf = f.img3;
%        HaralickFeat{n} = single(Hf);
%    end
%end

for i = 1:size(Hf,3)
    feat = Hf(:,:,i);
    img_mask = mask(:,:,1);
    roi = feat(img_mask ==1);
    
    Texturefeats(n+1) = mean(roi(:));
    Texturefeats(n+2) = std(roi(:));
    %Texturefeats(n+3) = mode(roi(:));
    n = n + 2;
end


% count = 1;
% modifier = [{'mean '} {'standard deviation '}]
% for j = 1:numel(feat.names)
% for i = 1:numel(modifier)
%     HaralickTextureFeatureDescription{count} = [modifier{i} 'intensity ' feat.names{j}];
%     count = count + 1;
% end
% end
