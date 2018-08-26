%%
%do_elsd=1;
test_id=3;
if do_elsd
    %im = im2double(rgb2gray(imread('~/Desktop/Lines/syn3.jpg')));    
    switch test_id
        case 0
            % synthetic line
            sz=[400 400];center=[50 300 50 100];r=10;
            dist = sqrt(sum((center(1:2)-center(3:4)).^2));
            im = zeros(sz);[yy,xx]=meshgrid(1:sz(2),1:sz(1));
            pts = [ xx(:) yy(:)]; tmp1 = zeros(size(pts,1),1);
            len= cross([bsxfun(@minus,pts,center(1:2)) tmp1],[bsxfun(@minus,pts,center(3:4)) tmp1],2);
            im = reshape(sigmf(sqrt(sum(len.^2,2))/dist -r,[20 5]),sz);
            %imshow(im);return
        case 1
            % synthetic circle
            sz=[400 400];center=[200,100];r=[40 40];
            im = zeros(sz);[yy,xx]=meshgrid(1:sz(2),1:sz(1));
            im = sigmf(((xx-center(1))/r(1)).^2+((yy-center(2))/r(2)).^2-1,[20 5]);
        case 2
            % synthetic ellipse
            sz=[400 400];center=[200,100];r=[40 20];
            im = zeros(sz);[yy,xx]=meshgrid(1:sz(2),1:sz(1));
            im = sigmf(((xx-center(1))/r(1)).^2+((yy-center(2))/r(2)).^2-1,[20 5]);
        case 3
            %name= '/home/Stephen/Dropbox/neal_ellipse/HDR_Crescent_Moon_Earthshine.JPG';
            name= '/Users/Stephen/Dropbox (MIT)/DeviationMagnification/Data/Saturn/20150104/PIA18277.tif';
            im = im2double(imread(name));
    end
    %[a,b,c]=elsd_mex(255*im,1,0,'haha.txt');
    [a,b,c]=elsd_mex_circle(255*im,1,0,'haha.txt');
    disp([a b c])
    d=load('haha.txt');
end
cla
%imagesc(im),hold on
imshow(im),hold on
for id=1:size(d,1);
    geo_type = d(id,1);
    geo_param = d(id,2:end);
    U_plotGeo(geo_type,geo_param)
end
len=U_lenGeo(d(:,1),d(:,2:end))
axis tight,axis equal

%{


len=U_lenGeo(d(:,1),d(:,2:end));
[~,len_id]=sort(len(:)','descend');
for id = len_id
cla,imshow(im),hold on,U_plotGeo(d(id,1),d(id,2:end))
title(id)
d(id,:)
waitforbuttonpress
end

switch geo_type
    case 0
        U_plotLine(im,geo_param);
    case {1,2}
        U_plotCircle(im,geo_param);
    case 2
        U_plotEllipse(im,geo_param);
end

%}