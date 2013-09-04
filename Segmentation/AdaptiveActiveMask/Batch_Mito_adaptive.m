% COPYRIGHT
%       This file is part of the Matlab code provided for the following paper:
%
%		Kuan-Chieh Jackie Chen, Yiyi Yu, Ruiqin Li, Hao-Chih Lee, Ge Yang, Jelena Kovacevic,
%		"Adaptive active-mask image segmentation for quantitative characterization of 
%		mitochondrial morphology,"
%		2012 19th IEEE International Conference on Image Processing (ICIP), pp.2033-2036, Sept. 30 2012-Oct. 3 2012
%
%       Authors: Kuan-Chieh Jackie Chen
%		Data Created: Jan 26, 2012
% 		Last Modified: Feb 22, 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi = Batch_Mito_adaptive( f_mito, savepath, fname, gammas, ad_sigma_maxs )
% psi = Batch_Mito_adaptive( f_mito, savepath, fname, gammas, ad_sigma_maxs )
%
% PURPOSE:
% A wrapper for batch processing using Adaptive Active Masks
%
% REQUIRE:
% AdaptiveActiveMasks.m
% outline.m
% 
% INPUT(S):
% (1) f_mito = A single image
% (2) savepath = A path of the directory to store the results including masks as
% matfile and visualization of the results
% (3) fname = A file name prefix for the results
% (4) gammas = A list or a single value of the gamma parameter ( suggested
% values to try: [-5:-5:-30] )
% (5) ad_sigma_maxs = A list or a single value of the adaptive sigma
% parameter ( suggested values to try: [5] )
%
% OUTPUT(S):
% (1) psi = The resulting masks
%
% NOTES:
% None

% Parameters
P = 256; % initial # of Mask
K = 1; % initial resolution level 
K_0 = 0; % stop resolution level (lower, more refinement)
a = 3.5;
b = 8;
alpha = -0.9;
beta=1;
       
% pad image for ActiveMask
org_size = size(f_mito);
new_size = ceil(org_size./(2^K)).*(2^K);
pre_pad = floor((new_size - org_size) / 2);
post_pad = (new_size-org_size)-pre_pad;
paded_f_mito = padarray(f_mito, pre_pad, 'replicate', 'pre');
paded_f_mito = padarray(paded_f_mito, post_pad, 'replicate', 'post');

% loop over different parameters of gamma and adaptive sigma
for gamma = gammas
    for ad_sigma_max = ad_sigma_maxs
        AAM_start=tic;
        psi = AdaptiveActiveMasks(paded_f_mito, P, K, K_0, a, b, alpha, beta, gamma, ad_sigma_max);
        disp(sprintf('Approximate time taken: %2.2f seconds',toc(AAM_start)));
        
%       reassign the masks to have index from 0 to # of masks
        idx = 0;
        unq_label = unique(psi);
        for seg_label = 1:length(unq_label)
            psi(psi==unq_label(seg_label)) = seg_label-1;
            idx = idx+1;
        end
        
%       resize the psi to match the size of f_mito
        k_ratio = size(paded_f_mito,1) / size(psi,1)
        psi = kron(psi, ones(k_ratio, k_ratio) );
        psi = psi(pre_pad(1)+1:end-post_pad(1),pre_pad(2)+1:end-post_pad(2));
        % figure,highlight(f_mito, new_psi);

%       save the masks as mat  
        matname = sprintf('%s/%s_psi_%d_%d',savepath,fname,ceil(gamma),ceil(ad_sigma_max));
        save( [matname '.mat'], 'psi' );
        
%       generate visualization of the results
        imwrite(outline(f_mito,psi,[0,.5]), [matname '.png'], 'png');
    end
end