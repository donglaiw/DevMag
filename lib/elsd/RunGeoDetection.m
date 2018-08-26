%%
%name= 'cameraman.tif';
%name = '~/Downloads/pia18308_full.jpg';
name='/Users/Stephen/Dropbox (MIT)/Proj/DeviationMagnification/Data/Benchmark/window_3.jpg';
% gray double image from 0-1
im = im2double(imread(name));
if numel(size(im))==3
    im = rgb2gray(im);
end

geo_id = 0;% all
%geo_id = 1;% line
%geo_id = 2;% circle
% angle_thres,density_thres,nfa_thres
% angle_thres (0-180): C/sin(angle_thres)
% density_thres (0-1): for refinement
% nfa > -log10(thres(2))
thres=[22.5,1,1e1];
%thres(1)=75;

do_dsp=0;
do_dsp=1;

switch geo_id
    case 0
        num = elsd_mex(255*im,thres,0.8,'p_geo.txt');
        param = load('p_geo.txt');
        fprintf('Found %d segments.\n',num);
    case 1
        num = elsd_mex_line(255*im,thres,0,'p_line.txt');
        param = load('p_line.txt');
        fprintf('Found %d line segments.\n',num);
    case 2
        [~,num]=elsd_mex_circle(255*im,thres,0,'p_circle.txt');
        param = load('p_circle.txt');
        fprintf('Found %d circular segments.\n',num);
end

if do_dsp
    arc_len=U_lenGeo(param(:,1),param(:,2:end));
    [~,len_id]=sort(arc_len(:),'descend');
    for id = len_id'
        cla,imshow(im),hold on,U_plotGeo(param(id,1),param(id,2:end))
        title(sprintf('id: %d',id))
        val = param(id,:);
        waitforbuttonpress
    end
end
