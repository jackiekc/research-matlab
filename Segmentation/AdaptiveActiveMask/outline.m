% COPYRIGHT
%       This file is part of the Matlab code provided for the following paper:
%
%		Kuan-Chieh Jackie Chen, Yiyi Yu, Ruiqin Li, Hao-Chih Lee, Ge Yang, Jelena Kovacevic,
%		"Adaptive active-mask image segmentation for quantitative characterization of 
%		mitochondrial morphology,"
%		2012 19th IEEE International Conference on Image Processing (ICIP), pp.2033-2036, Sept. 30 2012-Oct. 3 2012
%
%       Authors: Kuan-Chieh Jackie Chen
%		Data Created: 2011
% 		Last Modified: 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function img = outline(im,seg,adjust)
% img = outline(im,seg,adjust)
%
% PURPOSE:
% Nicely highlights an image with a given segmentation and adjusts it's contrast

if nargin == 4
    seg = bwperim(seg);
end
if nargin == 2
    adjust = [0;1];
end

num = max(unique(seg(:)));

c = rgb2hsv(ind2rgb(seg+1,hsv2rgb([0 0 0; linspace(0,1,num+1)' ones(num+1,2)])));
value =  double(im);
value =  value/max(value(:));
value = imadjust(value,adjust);
% value = zeros(size(im));
value(seg~=0) = 1;

sat = ones(size(im));
sat(seg==0) = 0;
hue = c(:,:,1);

labeled(:,:,1) = hue;
labeled(:,:,2) = sat;
labeled(:,:,3) = value;

img = hsv2rgb(labeled);