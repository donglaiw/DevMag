% function to warp images with different dimensions
function [warpI2,mask]=U_warpImage(im,vx,vy)

[height2,width2,nchannels]=size(im);
[height1,width1]=size(vx);

[xx,yy]=meshgrid(1:width2,1:height2);
[XX,YY]=meshgrid(1:width1,1:height1);
XX=XX+vx;
YY=YY+vy;
mask=XX<1 | XX>width2 | YY<1 | YY>height2;
XX=min(max(XX,1),width2);
YY=min(max(YY,1),height2);
%{
% dilate boundary for nan image
if(nnz(isnan(im))~=0)
    ind = isnan(im);
    ind2 = ~isnan(im);
    pre_im=im;
    im(ind)=0;
    H=fspecial('gaussian',10,0.5);
    im =imfilter(im,H,'replicate');
    im(ind2)=pre_im(ind2);
end
%}
for i=1:nchannels
    foo=interp2(xx,yy,im(:,:,i),XX,YY,'bilinear');
    warpI2(:,:,i)=foo;
end

if max(warpI2(:))>5
    warpI2(warpI2>255)=255;
else
    warpI2(warpI2>1)=1;    
end
warpI2(warpI2<0)=0;
%mask=1-mask;
