% ComputeBlur smoothes the image and downsamples by a factor of 2
%
% USAGE     
%
%       smallImage=ComputeBlur(bigImage)
%
% INPUTS    
%
%       bigImage The image to be smoothed and downsampled
%
% COPYRIGHT
%
%       This file is part of the Matlab code provided for the following
%       reproducible paper:
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
%       Latest modifications: November 9, 2007.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function smallImage = ComputeBlur(bigImage)
bigImage = 0.25*(bigImage + horzcat(bigImage(:,2:end), bigImage(:,1)) + vertcat(bigImage(2:end,:),bigImage(1,:)) + [bigImage(2:end,2:end) bigImage(2:end,1); bigImage(1,2:end) bigImage(1,1)]);
smallImage = bigImage(1:2:end,1:2:end);