function [feats,description,CCGinfo,CGinfo] = Lextract_all_features(bounds,para,img)
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
% 
idxF=1;
fprintf('\nExtracting Graph Features...')
feats{idxF} = extract_graph_feats(bounds);
description{idxF} = [o.Description.ImageFeatures(1:51)];
idxF=idxF+1;

fprintf('\nExtracting Morph Features...')
feats{idxF} = extract_morph_feats(bounds);
temp = regexp(o.Description.ImageFeatures,'Morph');
temp2 = cellfun(@isempty,temp,'UniformOutput',false);
description{idxF} = [o.Description.ImageFeatures([temp2{:}] == 0)];
% idxF=idxF+1; % don't need this, because the later process will take into
% account this idx

curpara.curName=para.curName;
curpara.CCGlocation=para.CCGlocation;

fprintf('\nExtracting CGT Features...')
CGinfo=[]; % cell graph information
set_alpha=[para.CGalpha_min:para.alpha_res:para.CGalpha_max];
for f=1:length(set_alpha)
    curpara.alpha=set_alpha(f); curpara.radius=para.radius;
    
    [feats{f+idxF},CGinfo{f+idxF}]= extract_CGT_feats(bounds,curpara);
    temp = regexp(o.Description.ImageFeatures,'CGT');
    temp2 = cellfun(@isempty,temp,'UniformOutput',false);
    temp= [o.Description.ImageFeatures([temp2{:}] == 0)];
    
    for i=1:length(temp)
        cur=temp{i};
        str=sprintf('_a=%.2f',curpara.alpha);
        cur(end+1:end+length(str))=str;
        temp{i}=cur;
    end
    description{f+idxF} =temp;
end
idxF=idxF+length(set_alpha);

curpara.curName=para.curName;
curpara.CCGlocation=para.CCGlocation;

fprintf('\nExtracting Cluster Graph Features...');
% note that this should be build on the cell clusters!!!
% bounds.CellCluster
CCbounds.centroid_c=bounds.CellClusterC_c;
CCbounds.centroid_r=bounds.CellClusterC_r;

CCGinfo=[];% cell cluster graph information
set_alpha=[para.CCGalpha_min:para.alpha_res:para.CCGalpha_max];
for f=1:length(set_alpha)
    curpara.alpha=set_alpha(f); curpara.radius=para.radius;
    [CCGfeats,~,CCGinfo{f}]=extract_cluster_graph_feats(CCbounds,curpara);
    feats{f+idxF}=CCGfeats;

%     feats{f} = extract_cluster_graph_feats(CCbounds,curpara);
    temp = regexp(o.Description.ImageFeatures,'GSG');
    temp2 = cellfun(@isempty,temp,'UniformOutput',false);
    temp = [o.Description.ImageFeatures([temp2{:}] == 0)];
    for i=1:length(temp)
        cur=temp{i};
        cur(1:3)='CCG';
        str=sprintf('_a=%.2f',curpara.alpha);
        cur(end+1:end+length(str))=str;
        temp{i}=cur;
    end
    description{f+idxF} =temp;
end
idxF=idxF+length(set_alpha);



% if nargin > 2
%     fprintf('\nExtracting Texture Features...')
%     feats{5} = extract_texture_feats(img);
%     temp = regexp(o.Description.ImageFeatures,'Haralick');
%     temp2 = cellfun(@isempty,temp,'UniformOutput',false);
%     description{5} = [o.Description.ImageFeatures([temp2{:}] == 0)];
% end

%feat_description


function [graphfeats] = extract_graph_feats(bounds)
gb_r = [bounds.centroid_r];
gb_c = [bounds.centroid_c];

if length(gb_r)<3
    graphfeats=zeros(1,51);
else
    [graphfeats] = get_graph_features([gb_r]',[gb_c]');
end

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



function [CGTfeats, CGinfo] = extract_CGT_feats(bounds,para)
%% CGT
% addpath(genpath(['Z:\2012_10_GlandOrientation']))

%[CGTfeats,c_matrix,info] = get_CGT_features(bounds);
a = para.alpha;
r = para.radius;

%[CGTfeats, c_matrix, info, feats, network, edges] = get_CGT_features_networks_weight(bounds,a,r);
[CGTfeats, ~, ~, ~, ~, ~,VX,VY,~,~] = extract_CGT_features(bounds,a,r,para);
CGinfo.VX=VX;
CGinfo.VY=VY;
% CGinfo.x=x;
% CGinfo.y=y;
% CGinfo.edges=edges;

function [CCGfeats,feature_list,CCGinfo] = extract_cluster_graph_feats(bounds,para)
%CCG

%addpath(genpath(['Z:\2012_10_GlandOrientation']))

% build graph
alpha = para.alpha;
r = para.radius;

% alpha = curpara.alpha;
% r = curpara.radius;

%% buid ccg
if ~(isfield(para,'curName')||isfield(para,'CCGlocation'))
    error('please specify the image name for ccg');
else    
    str=sprintf('%s_CCG_%s_a=%.2f.mat',para.curName,para.CCGlocation,alpha);
    if exist(str,'file')~=2      
        % build graph
        [VX,VY,x,y,edges] = construct_ccgs(bounds,alpha, r);
        save(str,'VX','VY','x','y','edges','-v7.3');
    else
        load(str);
    end
end

% [VX,VY,~,~,edges] = construct_ccgs(bounds,alpha, r);
CCGinfo.VX=VX;
CCGinfo.VY=VY;
% CCGinfo.x=x;
% CCGinfo.y=y;
% CCGinfo.edges=edges;
% [VX,VY,x,y,edges] = Lconstruct_ccgs_optimized(bounds,alpha, r);

%[CCGfeats,feature_list] = cluster_graph_features_networkboost(bounds, edges);
% if length(bounds.centroid_c)~=length(edges)
%     disp('d');
% end
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
