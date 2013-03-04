function highlight(im,seg)
% nicely highlights an image with a given segmentation
%
% Mike McCann 2011

% num = length(unique(seg(:)));
num = max(unique(seg(:)));

c = rgb2hsv(ind2rgb(seg+1,hsv2rgb([0 0 0; linspace(0,1,num+1)' ones(num+1,2)])));
% value =  double(im-min(im(:)))/double(max(im(:)));
value =  double(im)/double(max(im(:)));
% value = ones(size(im));
sat = ones(size(im));
sat(seg==0) = 0;
hue = c(:,:,1);

labeled(:,:,1) = hue;
labeled(:,:,2) = sat;
labeled(:,:,3) = value;

imshow(hsv2rgb(labeled));