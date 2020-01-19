function [major_axis minor_axis eigval]=  fitEllipseToBoundary(bounds)
%FITELLIPSETOBOUNDARY calculates the axes of an ellipse fit to a boundary
% for a set of objects
%
% [major_axis minor_axis]=  fitEllipseToBoundary(bounds).
% bounds is a struct with bounds().r storing the x-axis location of points
% and bounds().c storing the y-axis location of points.
%
% Returns the major and minor axis of the object as determined by fitting an
% ellipse to the boundary of the object
%
% Author: Rachel E. Sparks 
% Date:   October, 20, 2012

    major_axis = zeros([ length(bounds) 2]);
    minor_axis = zeros ([length(bounds) 2]);
    eigval = zeros ([length(bounds) 2]);

    for i=1:length(bounds)
%         [components,~,e] = princomp([bounds(i).r; bounds(i).c]');
%         major_axis(i,:) = components(:,1);
%         minor_axis(i,:) = components(:,2);

        
       [components,~,e] = pca([bounds(i).r; bounds(i).c]'); % changed to different version of pca
        major_axis(i,:) = components(:,1)';
        minor_axis(i,:) = components(:,2)';

        eigval = e(1);
        
    end

end

%debug
%plot(bounds(i).r, bounds(i).c)
%hold on; 
%quiver(bounds(i).r, bounds(i).c, components(i,1), components(i,2))