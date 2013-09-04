% COPYRIGHT
%       This file is part of the Matlab code provided for the following paper:
%
%		Kuan-Chieh Jackie Chen, Yiyi Yu, Ruiqin Li, Hao-Chih Lee, Ge Yang, Jelena Kovacevic,
%		"Adaptive active-mask image segmentation for quantitative characterization of 
%		mitochondrial morphology,"
%		2012 19th IEEE International Conference on Image Processing (ICIP), pp.2033-2036, Sept. 30 2012-Oct. 3 2012
%
%       Authors: Kuan-Chieh Jackie Chen
% 		Last Modified: 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source_foldername = 'ctl_20110927_larva1_1_a3_1';

fimgmat = 'example.mat';
if exist(fimgmat,'file')
    load(fimgmat);
else
    imgs = {};
    imglist = dir(fullfile(source_foldername, '*.tif'));
    for i = 1:length( imglist )
        imgfn = fullfile(source_foldername,imglist(i).name);
        imgs{i} = struct('name',imglist(i).name(1:end-4),'folder',source_foldername,'img',imread(imgfn));
    end
    imgs = [imgs{:}];
    save(fimgmat, 'imgs');
end

masks = {};
for i = 1:length(imgs)
    i
    savepath = fullfile( 'example_psis', imgs(i).folder );
    
    result = {};
    gammas = [-15];
    ad_sigma_maxs = [5];
    if( ~exist( savepath, 'dir' ) )
        mkdir( savepath );
    end

    masks{i} = Batch_Mito_adaptive(imgs(i).img, savepath, imgs(i).name,gammas,ad_sigma_maxs);

end