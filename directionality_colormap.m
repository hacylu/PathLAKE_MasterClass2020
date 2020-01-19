function [angle_degrees,angle_color] = directionality_colormap(bounds,scale, flip)
% plot directionality colormap for all boundary segmentations
% bounds is a struct with bounds().r storing the x-axis location of points
% and bounds().c storing the y-axis location of points.(should be row
% vector)
% George Lee (2013)

%            imshow(img)
%             hold on
%             directionality_colormap(bounds,scale,0);
% 
% "bounds" is a struct with bounds.r and bounds.c holding the x and y coordinates 
% of the vertices of the nuclei edge points and bounds.centroid_r and bounds.centroid_c 
% holding the coordinates of each nuclei center.   "scale" is the factor by which to scale up 
% the resulting boundary.  I don't think this variable actually matters for these figures so scale = 1 
% should work. I've attached what I think is all the code needed to create these figures and some example 
% bounds structs contained in a larger variable.  Sometimes the arrows or nuclei outlines require flipping...
% it seems to vary image to image and I haven't figured out any order to the flipping.
%     
if nargin < 2
    scale = 1;
end

% if nargin < 2
%     flip = 0;
% end

lwidth = 1;
S = 0.4;
hold on

[principalaxis] = fitEllipseToBoundary(bounds);

% %plot axis on bounds
% [handle]=  plotAxis(bounds(i), axis(i,:))

if nargin > 2
    
    %plot vector field - note that quiver is weird and requires flipping
    for i = 1:length(principalaxis(:,1))
        if principalaxis(i,1) >0
            principalaxis(i,:) = -principalaxis(i,:);
        end
    end
    
    %quiver([bounds(:).centroid_c].*scale, [bounds(:).centroid_r].*scale, [principalaxis(:,2)]', [principalaxis(:,1)]',S,'k','linewidth',2)
    quiver([bounds(:).centroid_c].*scale, [bounds(:).centroid_r].*scale, [principalaxis(:,2)]', [principalaxis(:,1)]',S,'g','linewidth',lwidth)
else
    %adjust axis
    for i = 1:length(principalaxis(:,1))
        if principalaxis(i,1) >0
            principalaxis(i,:) = -principalaxis(i,:);
        end
    end
    
    %quiver([bounds(:).centroid_r].*scale, [bounds(:).centroid_c].*scale, [principalaxis(:,2)]', [principalaxis(:,1)]',S,'k','linewidth',2)
    quiver([bounds(:).centroid_r].*scale, [bounds(:).centroid_c].*scale, [principalaxis(:,2)]', [principalaxis(:,1)]',S,'g','linewidth',lwidth)
    
end

angle_degrees = 180/pi *atan(principalaxis(:,2)./principalaxis(:,1))+90;

for i = 1:length(bounds)
    
    
    bin = 18;
    w = 180/bin;
    angle_color(i) = ceil(angle_degrees(i)/w);
    
    cmap = colormap(jet(bin));
    
    % plot boundary
    %plot([bounds(i).c bounds(i).c(1)],[bounds(i).r bounds(i).r(1)],'LineWidth',1.5, 'Color', [.0 .9 .0]);
    %plot([bounds(i).c bounds(i).c(1)],[bounds(i).r bounds(i).r(1)],'LineWidth',1.5, 'Color', [.0 .0 1]);
    %plot([bounds(i).c bounds(i).c(1)],[bounds(i).r bounds(i).r(1)],'LineWidth',1.5, 'Color',  cmap(angle_color,:));
    plot([bounds(i).c(:)' bounds(i).c(1)].*scale,[bounds(i).r(:)' bounds(i).r(1)].*scale,'LineWidth',2, 'Color',  cmap(angle_color(i),:));
%     plot([bounds(i).r(:)' bounds(i).r(1)].*scale,[bounds(i).c(:)' bounds(i).c(1)].*scale,'LineWidth',lwidth, 'Color',  cmap(angle_color(i),:));
    
    hold on
    % fill in color
    %fill([bounds(i).c bounds(i).c(1)],[bounds(i).r bounds(i).r(1)], cmap(angle_color,:));
    
    %     % colored vectors
    %     %S = 0.3
    %     S = 20
    %     hold on
    %         quiver([bounds(i).centroid_r], [bounds(i).centroid_c], [principalaxis(i,2)]', [principalaxis(i,1)]',S, 'Color', cmap(angle_color,:),'linewidth',3)
    
    
end

% function [major_axis minor_axis]=  fitEllipseToBoundary(bounds)
% %FITELLIPSETOBOUNDARY calculates the axes of an ellipse fit to a boundary
% % for a set of objects
% %
% % [major_axis minor_axis]=  fitEllipseToBoundary(bounds).
% % bounds is a struct with bounds().r storing the x-axis location of points
% % and bounds().c storing the y-axis location of points.
% %
% % Returns the major and minor axis of the object as determined by fitting an
% % ellipse to the boundary of the object
% %
% % Author: Rachel E. Sparks
% % Date:   October, 20, 2012
%
% major_axis = zeros([ length(bounds) 2]);
% minor_axis = zeros ([length(bounds) 2]);
%
% for i=1:length(bounds);
%     components = princomp([bounds(i).r; bounds(i).c]');
%     major_axis(i,:) = components(:,1);
%     minor_axis(i,:) = components(:,2);
% end
