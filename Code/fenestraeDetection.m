function [fenestrae,statsFenestrae] = fenestraeDetection(img0, method, post_processing,thresSeg,sizeLOG)


if ~exist('method','var')
    method = 5;
end
if ~exist('post_processing','var')
    post_processing = 1;
end
if strcmp(method,'threshold')
    method=1 ;
elseif strcmp(method,'adaptive threshold')
    method=2;
elseif strcmp(method,'closing')
    method=3;
elseif strcmp(method,'LoG kernels')
    method=4 ;
end

% This is the pre-processing section as described in the report

switch method
    case 1
        
        med_img = medfilt2(img0, [20,20]);
        
        img = img0 - med_img;
        fenestrae = -img;
        
        T = 0.72;
        R = 3;
        
    case 2
        
        med_img = medfilt2(img0, [20,20]);
        
        img = img0 - med_img;
        n = 30;
        mean_img = medfilt2(img, [20,20]);
        sd_minus_img = sqrt(conv2(subplus(img - mean_img).^2, ...
            1/(n^2) * ones(n), 'same'));
        
        if sd_minus_img ~= 0
            fenestrae = subplus(mean_img - img) ./ (sd_minus_img);
        else
            fenestrae = subplus(mean_img - img);
        end
        
        T = 0.47;
        R = 3;
        
    case 3
        
        med_img = medfilt2(img0, [20,20]);
        
        img = img0 - med_img;
        fenestrae = (imclose(img, strel('disk',4,8)) - img);
        
        T = 0.54;
        R = 4;
        
    case 4
        
        med_img = medfilt2(img0, [20,20]);
        
        img = img0 - med_img;
        sigma_max = 7;
        kernel = mexicanKernel(sigma_max, 8*sigma_max);
        for sigma = 3:2:(sigma_max)
            kernel = kernel + mexicanKernel(sigma, 8*sigma_max);
        end
        kernel = 2 * kernel / sigma_max;
        
        fenestrae = subplus(-conv2(img, kernel, 'same'));
        
        T = 0.34;
        R = 2;
        
        
    case 5
        if ~exist('thresSeg');  thresSeg = 0.275;   end
        if ~exist('sizeLOG');   sizeLOG  = 13;      end
        
        qq1     =(imfilter(img0,fspecial('log',sizeLOG*[1 1],4),'replicate'));
        qq2     = bwlabel(qq1>(max(qq1(:))*thresSeg));
        qq3     = regionprops(qq2,'area','majoraxis','minoraxis','eccentricity');

        qq4     = ismember(qq2,find(([qq3.Area]>9 )& ([qq3.Eccentricity]<0.9 )));
        post_processing =0;
                        post_processing =1;T=0.5;R=2;
       fenestrae = qq4;
    
    case 6
        
        [rows,cols,levs]= size(img0);
        qq1=(imfilter(img0,fspecial('log',[9 9],5),'replicate'));
        qq2=qq1';
        
        [P1,locs1] = findpeaks(qq1(:),'minpeakdistance',5,'minpeakwidth',5,'maxpeakwidth',35,'minpeakheight',max(qq1(:))/3);
        [P2,locs2] = findpeaks(qq2(:),'minpeakdistance',5,'minpeakwidth',5,'maxpeakwidth',35,'minpeakheight',max(qq1(:))/3);
        [cc1,rr1] = ind2sub([rows cols],locs1);
        [rr2,cc2] = ind2sub([cols rows],locs2);
        
        qq1b=zeros(size(qq1));
        qq2b = qq1b';
        
        qq1b(locs1)=1;
        qq2b(locs2)=1;
        qq2c=qq2b';
        qq3 =(qq2c+qq1b)>1;
        [cc3,rr3]=find(qq3);
        
        qq4= imdilate(qq3,strel('disk',4));
        
       fenestrae = qq4;
                post_processing =1;T=0.5;R=3;
                
                
                
                        
    case 7
       
        qq1     =round(1000*(imfilter(img0,fspecial('log',11*[1 1],5),'replicate')));
        P = houghpeaks(qq1,1000,'NHoodsize',17*[1 1],'threshold',0.33*max(qq1(:)));
        
        
        
        locs = sub2ind(size(qq1),(P(:,1)),(P(:,2)));
        qq2=zeros(size(qq1));
        
        qq2(locs)=1;
        
        
        
        locs = sub2ind(size(qq1),(P(:,1)),(P(:,2)));
        qq2=zeros(size(qq1));
        
        qq2(locs)=1;
        
        qq3=bwlabel((qq1>(0.33*max(qq1(:)))));
        qq5 = qq3.*qq2;
        qq6 = qq5(:);
        qq7 = unique(qq6);
        qq8 = histcounts(qq6,max(qq7));
        qq9 = ismember(qq3,find(qq8(2:end)>1));
        qq2 = qq2.*(1-qq9);
        
        qq4= imdilate(qq2,strel('disk',5));
        
        %

       fenestrae = qq4;
        post_processing =0;
                                post_processing =1;T=0.5;R=1;
        
end
fenestrae = (fenestrae-min(fenestrae(:))) / (max(fenestrae(:))-min(fenestrae(:)));

% post processing
if post_processing
    fenestrae = imdilate(fenestrae > T, strel('disk', R, 8));
end

[fenestrae_L,numFenestrae]= bwlabel(fenestrae);
fenestrae_P     = regionprops(fenestrae_L,'Area');
statsFenestrae.numFenestrae   = numFenestrae;
statsFenestrae.sizeFenestrae  = mean([fenestrae_P.Area]);
statsFenestrae.areaFenestrae  = sum([fenestrae_P.Area]);

end

