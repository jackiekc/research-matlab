function img = outline(im,seg,adjust,perm)
if nargin == 4
    seg = bwperim(seg);
end
if nargin == 2
    adjust = [0;1];
end

% seg = bwlabel(seg);
% % num = length(unique(seg(:)));
num = max(unique(seg(:)));
% 
% perm = randperm(num);
% for i = 1:num
%     seg(seg==i) = perm(i);
% end

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