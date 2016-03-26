function [ num ] = hole( image )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%image= bwmorph(image,'skel',Inf);
filled_image = imfill(image, 'holes');
b=filled_image-image;
% % Count how many unique objects there are
 [L,num] = bwlabel(b);

end

