% this example was made by Cheng Lu, Case Western Reserve University
% 2020 Jan 11.
%% read image and pre-segmented nuclear mask
if ispc
    strPath=[pwd '\imgs\'];
else
    strPath=[pwd '/imgs/'];
end

strI=sprintf('%sNon_120.mat',strPath);
load(strI);

strNu=sprintf('%sNon_120.mat_nuclei.mat',strPath);
load(strNu);
%
% strCluster=sprintf('%sNon_MSCellCluster_BW30_IM%d',strPath,i);
% load(strCluster);
%
%     load(sprintf('%s%d_CCGInfo_a=0.5.mat',strPath,i));
maskES=imread(sprintf('%sNon_120_mask.png',strPath));
maskES=imresize(maskES,[size(I,1),size(I,2)]);
%imshow(maskES);
%% display nuclear contours in epi/stroma regions
%%% some of the nuclear contours are not correct, please ignore that, this
%%% is just a show case of the concept :)

str_nuclei_label=sprintf('%snuclei_label_in_epistroma_IM120.mat',strPath);
load(str_nuclei_label);

[m,n,k]=size(I);
bwNuclei=zeros(m,n);
for kk = 1:length(nuclei)
    nuc=nuclei{kk};
    for kki=1:length(nuc)
        bwNuclei(nuc(kki,1),nuc(kki,2))=1;
    end
end
bwNuclei = logical(bwNuclei);% imshow(bwNuclei);
% bwNuclei = logical(imfill(bwNuclei,'holes'));% imshow(bwNuclei);
s=regionprops(bwNuclei,'Area','Solidity','Orientation','Centroid','MajorAxisLength','MinorAxisLength','ConvexHull');

figure(2);imshow(I);title('Nuclei on image');
hold on;
for k = 1:length(s)
    if nuclei_label(k)
        plot(nuclei{k}(:,2), nuclei{k}(:,1), 'g-', 'LineWidth', 1);
        curC=properties(k).Centroid;
        plot(curC(1),curC(2),'b.');
    end
    if ~nuclei_label(k)
        plot(nuclei{k}(:,2), nuclei{k}(:,1), 'c-', 'LineWidth', 1);
        curC=properties(k).Centroid;
        plot(curC(1),curC(2),'r.');
    end
    text(curC(1),curC(2),sprintf('%d',properties(k).Area),'FontSize',20);
end
hold off;
%% display the best fit of ellipses of nuclei, and the axes ratio of a best fit ellipse 
figure;
imshow(I);title('Best fit of ellipses of nuclei on image');

t = linspace(0,2*pi,50);

hold on
for k = 1:length(s)
    a = s(k).MajorAxisLength/2;
    b = s(k).MinorAxisLength/2;
    Xc = s(k).Centroid(1);
    Yc = s(k).Centroid(2);
    phi = deg2rad(-s(k).Orientation);
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    plot(x,y,'g','Linewidth',1);
    
    text(Xc,Yc,sprintf('%.2f',a/b),'FontSize',20);
end
hold off
%% display the convex hull of nuclei, and solidity 

% peroperties=regionprops(,'Orientation','Centroid','MajorAxisLength','MinorAxisLength');

figure;imshow(I);title('Nuclear convehull on image with solidity');
hold on;
for k = 1:length(s)
    
    curHull=s(k).ConvexHull;
    
    if nuclei_label(k)
        plot(nuclei{k}(:,2), nuclei{k}(:,1), 'g-', 'LineWidth', 1);
        curC=properties(k).Centroid;
        plot(curC(1),curC(2),'b.');
    end
    if ~nuclei_label(k)
        plot(nuclei{k}(:,2), nuclei{k}(:,1), 'y-', 'LineWidth', 1);
        curC=properties(k).Centroid;
        plot(curC(1),curC(2),'r.');
    end
    plot(curHull(:,1),curHull(:,2),'r-','LineWidth', 4);
    text(curC(1),curC(2),sprintf('%.2f',properties(k).Solidity),'FontSize',20);
end
hold off;
%% display Cell Cluster Graphs on image
ctemp=[properties.Centroid];
bounds.centroid_c=ctemp(1:2:end);
bounds.centroid_r=ctemp(2:2:end);

bounds.nuclei=nuclei(~nuclei_label);
[VX,VY,x,y,edges] = construct_ccgs(bounds,0.43, 0.2);

figure(3);imshow(I);title('CCG on image');hold on;
plot(VY', VX', 'y-', 'LineWidth', 2);
for k=1:length( bounds.centroid_r)
    plot(bounds.centroid_c(k),bounds.centroid_r(k),'b.','MarkerSize',10);
end
hold off;
%% show global graph (delaunay trianglation)
if ispc
    addpath([pwd '\George_GraphFeatureExtraction']);

else
    addpath([pwd '/George_GraphFeatureExtraction']);
end
para.show=1;
para.T_smallregion=1;
para.I=I;
maskES(:)=1;
get_graph_features_isoregion_forshowonly(bounds.centroid_c',bounds.centroid_r',maskES,para);

%% show global graph (Voronoi diagram)
figure(6);imshow(I);axis([1 size(I,1) 1 size(I,2)]);
hold on;
x = [bounds.centroid_r];
y = [bounds.centroid_c];

%plot(y,x,'rs','linewidth',2); hold on %plot centroids
plot(y,x,'rs','linewidth',1); hold on %plot centroids
[VY,VX] = voronoi(x,y);
%plot(VY,VX,'g-','linewidth',2)
plot(VX,VY,'g-','linewidth',1)
hold off;
%% show nuclear polarity on image
%%% turn nuclei to bounds strcture
bounds_np=[];
for k = 1:length(nuclei)
    bounds_np(k).c=nuclei{k}(:,2)';
    bounds_np(k).r=nuclei{k}(:,1)';
    curpro=properties(k).Centroid;
    bounds_np(k).centroid_c=curpro(1);
    bounds_np(k).centroid_r=curpro(2);
end
%%% display nuclear polariy
figure;imshow(I);hold on;
scale=1;
directionality_colormap(bounds_np,scale,0);
%% show texture features of nuclear surface 
addpath([pwd '\George_GraphFeatureExtraction\haralick']);
bwNuclei= false(size(I,1),size(I,2));
for tt=1:length(nuclei)
    curN=nuclei{tt};
    for ttt=1:length(curN)
        bwNuclei(curN(ttt,1),curN(ttt,2))=1;
    end
end
bwNuclei= imfill(bwNuclei,'holes');    %         imshow(bwNuclei);
% parameters of haralick feature calculation
info.dist = 1;
info.win = 1;
info.grays = 256;
n = 0;
gimg = rgb2gray(I); % assume img is rgb

grays = info.grays;
win = info.win;
dist = info.dist;
himg = uint16(rescale_range(gimg,0,grays-1));
% this process maybe time comsuming :)
f = haralick_img(himg,bwNuclei,grays,win,dist,1); 


haralick_im=f.img3;
colormap jet;
temp=haralick_im(:,:,8);% the 8th one is entropy of intesity 
% adjust the value to show on the image
haralick_showim=haralick_im(:,:,8);
temp(temp==0)=3;
Fmin=min(temp(logical(bwNuclei)));
haralick_showim(haralick_showim==0)=Fmin;
hmax=max(max(haralick_showim));
haralick_showim(haralick_showim~=hmax)=haralick_showim(haralick_showim~=hmax)*0.9;
figure;imshow(haralick_showim,[]);
colormap(jet(20));