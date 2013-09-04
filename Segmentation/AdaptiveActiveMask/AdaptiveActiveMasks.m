% COPYRIGHT
%       This file is part of the Matlab code provided for the following paper:
%
%		Kuan-Chieh Jackie Chen, Yiyi Yu, Ruiqin Li, Hao-Chih Lee, Ge Yang, Jelena Kovacevic,
%		"Adaptive active-mask image segmentation for quantitative characterization of 
%		mitochondrial morphology,"
%		2012 19th IEEE International Conference on Image Processing (ICIP), pp.2033-2036, Sept. 30 2012-Oct. 3 2012
%
%       Authors: Kuan-Chieh Jackie Chen
%       Last Modified: Jan 26, 2012.
%
%       This file is modified from the code provided in the following paper:
%
%       Gowri Srinivasa, Matthew C. Fickus, Yusong Guo, Adam D. Linstedt, Jelena Kovacevic,
%       "Active mask segmentation of fluorescence microscope cell images",
%       IEEE Transactions on Image Processing. 
%                                                                           
%       This program is free software; you can redistribute it and/or modify it  
%       under the terms of the GNU General Public License as published by the     
%       Free Software Foundation; either version 2 of the License, or (at your    
%       option) any later version. This software is distributed in the hope that  
%       it will be useful, but without any warranty; without even the implied     
%       warranty of merchantability or fitness for a particular purpose.          
%       See the GNU General Public License for more details                       
%       (enclosed in the file COPYING).   
% 
%       GNU General Public License,
%       Copyright (C) 2008,
%       bimagicLab, Center for Bioimage Informatics (CBI)
%       Carnegie Mellon University,
%       Pittsburgh, PA 15213, USA.
%                                                                           
% COMMENTS
%
%       Authors: Gowri Srinivasa and Matthew C. Fickus
%       Latest modifications: April 22, 2008.
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psi] = AdaptiveActiveMasks(f,M,K,K_0,a_max,b_in,alpha,beta,gamma,ad_sigma_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%f - image to be segmented
%
%M - (initial) number of masks (we rename this initP in the code as M,N are
%    used to denote the dimensions of f)
%
%K - number of decomposition levels; it is an integer in {0,...,Kmax}
%   Kmax=3 for DS-1 and DS-2
%   K=0 indicates the coarsest resolution and
%   K=3 is one-eighth the original resolution f
%
%K_0 - number of levels to which we want to refine the segmentation 
%         K_0=0 implies the segmentation is successively lifted to the 
%                orignial resolution; 
%         K_0=K implies that segmentation is not refined at all
%
%a - cell array of scale parameters for the lowpass filter for region-based 
%    distributing function; the dimensions of a are (K+1)XJ_k where J_k refers 
%    to the number of scales (annealing) at resolution k 
%    a{1} has J_K numbers to be used as scale parameters for h at the coarsest scale, K
%    a{K+1-K_0} has the scale parameters for h at the desired original resolution (K_0) of f
%
% b - scale parameter used in g, the lowpass filter for voting-based
% distributing fuction
%
%alpha - weight for the region-based distributing function; lies in (-1,0)
%
%beta - harshness of the threshold 
%
%gamma - (average) value of pixels on the foreground-background border
%
%
%Note - In the scripts:
%K_0 is given by "k" and "shrink": 
%  k is the desired degree of refinement (assuming a maximum decomposition level of Kmax=3)
%  shrink modifies the original image to the desired (lower) resolution
%beta is obtained from "4/(high-low)" 
%gamma is obtained from "(high+low)/2" 

%%%%%%%%%%%%%%%%%%%%%%%%%%%Output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%psi - collection of M masks, each number in psi represents a region in f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format compact;
epsilon=0.0001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization of the masks; more the randomness, greater the chance of
%achieving a good segmentation: lesser number of unwanted splits and merges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initP=M; %P will denote the number of masks 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img=f;%the original image 
org_img=f;
%zero pad the image to avoid the masks from wrapping over
% obsolete by gunhee
% xPad=16 ;
% yPad=16 ;
shrink = K - K_0 ;

a{1} = [a_max-3 a_max-4];
a{2} = [a_max-2 a_max-3];
a{3} = [a_max-1 a_max-2];
a{4} = [a_max-0 a_max-1];

% for l=1:shrink
%     img=ComputeBlur(img);
% end

% a{4}=[24 22 20];
% a{3}=[20 18 16];
% a{2}=[16 14 12];
% a{1}=20; %Originally fixed at [10 8 6], but was too slow. Then changed to 24, but was too smooth; 12 was too sharp (we got too many small fragments) Settled on 20 as it is faster than 16.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create an image mask (characteristic function for the image) and zero pad
%the original image and set the time counter on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make each image dimension power to 2. (zero-padding)
% imgMask = padarray(ones(size(img)),[xPad,yPad],0,'both');
% img=padarray(img,[xPad,yPad],0,'both');
size_new = get_size_power2(size(img)) ;
pre_pad_i = floor((size_new-size(img))/2) ;
post_pad_i = round((size_new-size(img))/2) ;
imgMask = padarray(ones(size(img)),pre_pad_i,0,'pre') ;
imgMask = padarray(imgMask,post_pad_i,0,'post') ;
img=padarray(img,pre_pad_i,0,'pre') ;
img=padarray(img,post_pad_i,0,'post') ;

fft_imgMask = fft2(imgMask) ;
fft_img = fft2(img) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Multiresolution loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k_iter=0;

for k=K:-1:K_0 %from K, the coarsest resolution to K_0, the desired resolution
%     k_iter=k_iter+1;
    dil=a{3-k+1};
%     dil=(2^(-shrink))*dil; %If the dilation factors are not scaled with resolution then the image gets oversmoothed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Multiscale loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:length(dil)
        a_k_j=dil(j);
        [Mpad,Npad]=size(img);
        [X,Y]=ndgrid([0:floor(Mpad/2),-floor(Mpad/2)+1:-1],[0:floor(Npad/2),-floor(Npad/2)+1:-1]);
        h=exp(-(1/(a_k_j^2))*(X.^2+Y.^2)); % gaussian window
        fft_h = fft2(h) ;
        hScale=ifft2(fft_imgMask.*fft_h);
        f=(ifft2(fft_img.*fft_h))./(hScale+epsilon);
%         f=f(xPad+1:end-xPad,yPad+1:end-yPad);
        f=f(pre_pad_i(1)+1:end-post_pad_i(1),pre_pad_i(2)+1:end-post_pad_i(2));
        org_f = org_img;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %MR (blur and downsample) the image to the desired extent (size==k)
        %Initilize a MXNXP dimensional function, chi into which the contour is
        %embedded. The corresponding mask of interest is extracted in Psi.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for l=1:k
            f=ComputeBlur(f);
            org_f=ComputeBlur(org_f);
        end
        [M,N]=size(f); %used to compute the desired dimension of f
        
        % fPad=padarray(ones(size(f)),[xPad,yPad],0,'both');
        % make each image dimension power to 2. (zero-padding)
        size_new = get_size_power2(size(f)) ;
        pre_pad_f = floor((size_new-size(f))/2) ;
        post_pad_f = round((size_new-size(f))/2) ;
        fPad = padarray(ones(size(f)),pre_pad_f,0,'pre') ;
        fPad = padarray(fPad,post_pad_f,0,'post') ;

        if (k==K) && (a_k_j==max(dil))%we have to intialize chi and psi for the very first run
            P=initP;
            maskP=P;
            maskIndex=1:P;
            chi=zeros(M,N,P); %multivalued function initlized at all zeros            
            %%%%%%%%%%%%%%%%%%%%%%%Random seeding%%%%%%%%%%%%%%%%%%%%%%%%%% 
            psi=ceil(P*rand(M,N)); %we introduce P numbers into the MXN mask Psi
%             psi=load('random_psi.txt');
%             psi = psi + 1 ;
%             load tmp_rand_psi.mat ;
            for p=1:P
                temp=chi(:,:,p);
                temp(psi==p)=1; %a characteristic function that identifies where psi has the value p
                chi(:,:,p)=temp; %the values in the pth dimension are set to the indicator function (above)
            end
            [maxval, psi] = max(chi,[],3); %the updated psi is computed as the maxarg of chi (in the pth dimension)
        elseif (k~=K) && (a_k_j==max(dil)) % we continue working with the chi and psi that exist at the previous scale
            newchi=zeros(M,N,P); %initialize a newchi to match the new scale
            for p=1:P
                newchi(:,:,p)=kron(chi(:,:,p),[1 1; 1 1]); %when the size of the image is doubled along each dimension
                %chi is dilated by dilating each of its values
            end
            clear chi;
            chi=newchi;
            newPsi=kron(psi,[1 1; 1 1]); %values of psi are dilated in the same way
            clear psi;
            psi=newPsi;
            clear newchi ;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adaptive Region-based distributing function 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %R is a coarse separation between foreground and
        %background regions with cusps between touching foreground objects
        ad_sigma = ad_sigma_max / 2^k;
%       Guassian filter
        estimated_bg = imfilter(double(org_f),fspecial('gaussian',ceil(ad_sigma*6), ad_sigma ), 'replicate' );

        R = alpha.*erf( beta*(f - (estimated_bg - gamma) ) );
%         matname = sprintf('R_%d_%d_%d',ceil(gamma),ceil(ad_sigma_max),k);
%         save( matname, 'R', 'f', 'estimated_bg', 'gamma', 'beta', 'alpha', 'ad_sigma_max', 'ad_sigma', 'k' );
%         figure,imshow(R,[]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Local majority voting-based distributing function 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         fPad=padarray(ones(size(f)),[xPad,yPad],0,'both');
        [MPad,NPad]=size(fPad);
%         fMask=padarray(ones(size(f)),[xPad,yPad],0,'both');
        
        % make each image dimension power to 2. (zero-padding)
        fPad = padarray(ones(size(f)),pre_pad_f,0,'pre') ;
        fPad = padarray(fPad,post_pad_f,0,'post') ;
        fMask = padarray(ones(size(f)),pre_pad_f,0,'pre') ;
        fMask = padarray(fMask,post_pad_f,0,'post') ;
        
        [X,Y]=ndgrid([0:floor(MPad/2),-floor(MPad/2)+1:-1],[0:floor(NPad/2),-floor(NPad/2)+1:-1]);

        %b can be an annealed just as a is annealed
        c=1/8; %c can also be adjusted if necessary (it is negligibly small)
        % gunhee
%         b=b_in*2^(-k-shrink);%64 for the original image
        b=b_in*2^(-k);%64 for the original image
        g=0.5*(1-erf((X.^2+Y.^2-b^2)/(b^2+c)));%smaller the "support", slower the computation
        ghat=fft2(g);
        gScale=ifft2(fft2(fMask).*ghat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Keep track of how many pixels are changing; the algorithm stops when none of the pixels (in psi) move
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        change=1;
        iter=0;
        while change~=0
            tic
            iter=iter+1 ;

            % chiPad=padarray(chi,[xPad,yPad],'replicate','both');
            chiPad = padarray(chi,pre_pad_f,'replicate','pre') ;
            chiPad = padarray(chiPad,post_pad_f,'replicate','post') ;

            for p = 1:P
                chiPad(:,:,p)=(ifft2(fft2(chiPad(:,:,p)).*ghat))./(gScale+epsilon);
                % chi(:,:,p)=chiPad(xPad+1:end-xPad,yPad+1:end-yPad,p);
                chi(:,:,p)=chiPad(pre_pad_f(1)+1:end-post_pad_f(1),pre_pad_f(2)+1:end-post_pad_f(2),p);
            end
            chi(:,:,1)=chi(:,:,1)+R;
            [maxval,inwhat]=max(chi,[],3);
            newPsi=zeros(size(psi));
            chi=zeros(M,N,P);
            newP=0;
            for p=1:P %as chi and psi evolve, defragment chi to get rid of space allotted to nonexisting dimensions
                temp=zeros(size(psi));
                temp(inwhat==p)=1;
                if(sum(sum(temp))~=0) %check if any dimension of chi has been wiped out
                    newP=newP+1;
                    %%%%%%%%Maintain the numbers assigned to the masks%%%%%
                    chi(:,:,newP)=temp; %defragmenting the embedding function chi
                    maskIndex(newP)=maskIndex(p); %as we get rid of non existing dimension, we need to reassign mask values to keep the hues from changing
                    newPsi=newPsi+maskIndex(p)*temp;
                end
            end
            change=sum(sum(((newPsi-psi)~=0))); %how many pixels in the mask have moved?            
            P=newP;
            psi=newPsi;
            disp(sprintf('k=%02d; a=%02d P=%03d; change=%05d; image %03d: Approximate time taken: %2.2f seconds',k,dil(j),P,change,iter,toc));
        end
    end
end

% save tmp_psi01.mat psi;
