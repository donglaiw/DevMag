function varargout = gui_pipeline(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_pipeline_OpeningFcn, ...
    'gui_OutputFcn',  @gui_pipeline_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


function varargout = gui_pipeline_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function gui_pipeline_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

guidata(hObject, handles);
% init from varargin and input text
handles=init_handles(handles,varargin);
% init magnification  directory
guidata(hObject, handles);

function handles=init_handles(handles,var_in)
handles.opt = var_in{2};
switch handles.opt
    case 1
        % do label
        handles.fn='';
        handles.DD_in='./';
        param_label = var_in{3};
        handles.imgs = param_label{1};
        
        ind_name = find(handles.imgs=='/',1,'last');
        if ~isempty(ind_name)
            handles.DD_in = handles.imgs(1:ind_name);
            handles.imgs= {handles.imgs(1+ind_name:end)};
            handles.imgs_r = param_label{2};
        end
        handles.DD_out = param_label{3};
        handles.DD_out_label = param_label{4};
        handles.snap_thres = param_label{5};
        
        param_ana = var_in{4};
        handles.param = param_ana;
        handles=init_mag(handles);
        handles=init_label(handles);
        handles=init_dir(handles);
end


%% section 1: image directory
function popupmenu3_Callback(hObject, eventdata, handles)
handles=init_dir(handles);
guidata(hObject, handles);

function handles=init_dir(handles)
% init
if ~isfield(handles,'imgs')
    handles.folders = dir(handles.DD_in);
    handles.folders = handles.folders(cell2mat({handles.folders.isdir}));
    handles.folders(1:2) = [];
    set(handles.popupmenu3,'String',{handles.folders.name});
    handles.fn = [handles.folders(get(handles.popupmenu3,'Value')).name '/'];
    handles.imgs = U_getims([handles.DD_in handles.fn]);
    handles.imgs = {handles.imgs.name};
end
handles.img_id = 1;

set(handles.sel_id,'String',num2str(handles.img_id));
handles=update_img(handles);

function NextImage_btn_Callback(hObject, eventdata, handles)
if handles.img_id == numel(handles.imgs)
    handles.img_id= 1;
else
    handles.img_id =handles.img_id+1;
end
set(handles.sel_id,'String',num2str(handles.img_id));
handles=update_img(handles);
guidata(hObject, handles);

function PrevImage_btn_Callback(hObject, eventdata, handles)
if handles.img_id == 1
    handles.img_id = numel(handles.imgs);
else
    handles.img_id = handles.img_id-1;
end
set(handles.sel_id,'String',num2str(handles.img_id));
handles=update_img(handles);
guidata(hObject, handles);

function sel_id_Callback(hObject, eventdata, handles)
% sanity check
new_id = min(numel(handles.imgs),max(1,round(str2double(get(handles.sel_id,'String')))));
if new_id~= handles.img_id
    handles.img_id = new_id;
    set(handles.sel_id,'String',num2str(handles.img_id));
    handles=update_img(handles);
end
guidata(hObject, handles);

function display_img(handles)
cla(handles.axes_im)
axes(handles.axes_im);
imagesc(handles.img);
axis equal,axis off

cla(handles.axes_label_im)
axes(handles.axes_label_im);
imagesc(handles.img);
axis equal,axis off


function handles=update_img(handles)
handles.img_name = handles.imgs{handles.img_id};
handles.img = im2double(imread([handles.DD_in handles.fn handles.img_name]));
if handles.imgs_r(handles.img_id)~=1
    handles.img = imresize(handles.img,handles.imgs_r(handles.img_id));
    handles.img(handles.img>1) = 1;
    handles.img(handles.img<0) = 0;
end
handles.img_amp=[];
handles.img_sz = size(handles.img);
set(findobj('Tag','text_ImgId'),'String',sprintf('(          / %d)',numel(handles.imgs)));
set(findobj('Tag','text_ImgName'),'String',handles.imgs{handles.img_id});
display_img(handles)
undisplay_patches(handles)
undisplay_mag(handles)
handles = U_loadLabel(handles);
handles = U_loadAmp(handles);


%% section 2: geometry label
function handles=init_label(handles)
% init geometry
handles.geoLast = 0;
handles.geos_id = [];
handles.geos = {};
handles.dsp_geo = 0;
handles.geo_patches = {};
handles.autogeos = {};
handles.geoID = Tag2labelType(get(get(handles.panel_labelType,'SelectedObject'),'Tag'));
% init display patch size
handles.label_psz_c = str2double(get(handles.label_corner,'String'));
handles.label_psz_s = str2double(get(handles.label_sample,'String'));
handles.label_ran_c = -handles.label_psz_c:handles.label_psz_c;
handles.patches = cell(1,5);

function newpos = U_patch_pos(handles,pos)
newpos = round(pos);
newpos = max(newpos,1);
newpos(:,1) = min(newpos(:,1),handles.img_sz(1));
newpos(:,2) = min(newpos(:,2),handles.img_sz(2));

function [c,t]=U_pt2circ(pts)
% hard to do for ellipse
% pts: nx2
t = pts(2,:)-pts(1,:); u = pts(3,:)-pts(1,:); v = pts(3,:)-pts(2,:);
w = cross(t,u);
t2 = sum(t.^2); u2 = sum(u.^2); w2 = sum(w.^2);
c = p1+(t2*sum(u.*v)*u-u2*sum(t.*v)*t)/(2*w2);
r = 1/2*sqrt(t2*u2*sum(v.^2)/w2);

function dm_next_Callback(hObject, eventdata, handles)
ind = find(~isempty(handles.geos));
if numel(ind)~=0
    if handles.dsp_geo >= ind(end)
        handles.dsp_geo = ind(1);
    else
        next_ind = ind(ind>handles.dsp_geo);
        handles.dsp_geo =next_ind(1);
    end
    tmp_pos = handles.geos{handles.dsp_geo}.getPosition();
    handles=update_patches(handles,{round(tmp_pos(:,[2 1])),[]});
    guidata(hObject, handles);
end

function handles = update_patches(handles,pos,dsp,p_ind)
if ~exist('dsp','var');dsp=1;end
if ~exist('p_ind','var');p_ind=1:2;end
% callback function
% update the patch around the end points
switch handles.geoID
    case 1
        % line
        handles.patches_pt = U_patch_pos(handles,pos{1});
    case 2
        % circle
        handles.patches_pt = pos{2};
end
%{
handles.patches_pt = round(pos);
handles.patches_pt = max(handles.patches_pt,1);
handles.patches_pt(:,1) = min(handles.patches_pt(:,1),handles.img_sz(1));
handles.patches_pt(:,2) = min(handles.patches_pt(:,2),handles.img_sz(2));
%}
handles=update_patches_width(handles,3,p_ind);
if dsp
    display_patches(handles,3);
end
% remember all the strips
handles.geo_patches{handles.geoLast} = handles.patches;

function handles=update_patches_width(handles,opt_w,p_ind)
if(opt_w==1 || opt_w==3)
    % update the corner patch
    if handles.geoID==1
        for i=p_ind
            if (handles.patches_pt(i,1)>handles.label_psz_c ...
                    && handles.patches_pt(i,1)<=handles.img_sz(1)-handles.label_psz_c ...
                    && handles.patches_pt(i,2)>handles.label_psz_c ...
                    && handles.patches_pt(i,2)<=handles.img_sz(2)-handles.label_psz_c)
                handles.patches{i} = handles.img( ...
                    handles.patches_pt(i,1)+handles.label_ran_c, ...
                    handles.patches_pt(i,2)+handles.label_ran_c,:);
            else
                handles.patches{i} = zeros(handles.label_psz_c*2+1,handles.label_psz_c*2+1,3);
                
                st_x = max(1,handles.patches_pt(i,1)-handles.label_psz_c);
                lt_x = min(handles.img_sz(1),handles.patches_pt(i,1)+handles.label_psz_c);
                st_y = max(1,handles.patches_pt(i,2)-handles.label_psz_c);
                lt_y = min(handles.img_sz(2),handles.patches_pt(i,2)+handles.label_psz_c);
                
                pst_x = handles.label_psz_c-(handles.patches_pt(i,1)-st_x)+1;
                plt_x = handles.label_psz_c+(lt_x-handles.patches_pt(i,1))+1;
                pst_y = handles.label_psz_c-(handles.patches_pt(i,2)-st_y)+1;
                plt_y = handles.label_psz_c+(lt_y-handles.patches_pt(i,2))+1;
                
                handles.patches{i}(pst_x:plt_x,pst_y:plt_y,:) = handles.img(st_x:lt_x,st_y:lt_y,:);
            end
        end
    end
end
if(opt_w==2 || opt_w==3)
    % update the sampled patch
    switch handles.geoID
        case 1
            [handles.patches{3},handles.patches{5}] = U_sampleIm(handles.img,1,handles.label_psz_s,handles.mag_N,handles.patches_pt(1,[2 1]),handles.patches_pt(2,[2 1]));
        case 2
            radius = handles.patches_pt(3);
            center = handles.patches_pt([2 1]);
            th_ran = handles.patches_pt(5:6);
            [handles.patches{3},handles.patches{5}] = U_sampleIm(handles.img,2,handles.label_psz_s,handles.mag_N,th_ran(1),th_ran(2),radius,center);
    end
    if handles.mag_matt
        offset = 1e-3;
        Mask = handles.patches{3}+offset;
        step=4;
        st=1;lt=3;
        %lt=ceil(handles.label_psz_s/step);
        %lt=ceil(2*handles.label_psz_s/step);
        
        Mask(st:lt,:,:) = 1;
        Mask(end-(lt:-1:st)+1,:,:) = 0;
        %Mask =  CreatMaskLine_0427(handles.patches{3}, param.p1, param.p2);
        handles.patches{4} =runMattingfun(handles.patches{3}+offset,Mask);
    end
end

function undisplay_patches(handles)
cla(handles.axes_label_pt1),axes(handles.axes_label_pt1),axis off
cla(handles.axes_label_pt2),axes(handles.axes_label_pt2),axis off
cla(handles.axes_label_patch),axes(handles.axes_label_patch),axis off
cla(handles.axes_label_pix),axes(handles.axes_label_pix),axis off
cla(handles.axes_label_matt),axes(handles.axes_label_matt),axis off


function display_patches(handles,opt_w)
if(opt_w==1 || opt_w==3)&&handles.geoID==1
    % display corner patch
    cla(handles.axes_label_pt1)
    axes(handles.axes_label_pt1);
    imagesc(handles.patches{1});
    set(handles.axes_label_pt1,'NextPlot','add')
    % center pt
    plot(handles.label_psz_c+1 ,handles.label_psz_c+1,'rx','LineWidth',5);
    axis equal,axis on,axis tight
    
    cla(handles.axes_label_pt2)
    axes(handles.axes_label_pt2);
    imagesc(handles.patches{2});
    set(handles.axes_label_pt2,'NextPlot','add')
    % center pt
    plot(handles.label_psz_c+1,handles.label_psz_c+1,'rx','LineWidth',5);
    axis on,axis tight
    
end
if(opt_w==2 || opt_w==3)
    % display sample patch
    cla(handles.axes_label_patch)
    axes(handles.axes_label_patch);
    imagesc(handles.patches{3});
    set(handles.axes_label_patch,'NextPlot','add')
    center = handles.label_psz_s*handles.mag_N+1;
    % center line
    plot(1:size(handles.patches{3},2),ones(1,size(handles.patches{3},2))*center,'r-','LineWidth',0.5);
    axis on,axis tight
    
    
    cla(handles.axes_label_pix)
    axes(handles.axes_label_pix);
    plot(rgb2gray(handles.patches{3}(center,:,:)));
    axis on,axis tight
end
if handles.mag_matt
    % matting result
    cla(handles.axes_label_matt)
    axes(handles.axes_label_matt);
    imagesc(handles.patches{4});
    axis on,axis tight
    colormap gray
end


%{
% full circle
function out = U_imcircle(axe,pts)
if ~exist('pts','var')
    out = imellipse(axe,'PositionConstraintFcn',@(pos) [pos(1) pos(2) max(pos(3:4)) max(pos(3:4))]);
else
    out = imellipse(axe,pts,'PositionConstraintFcn',@(pos) [pos(1) pos(2) max(pos(3:4)) max(pos(3:4))]);
end
%}
function out = U_imcircle(axe,pts)
if ~exist('pts','var')
    out = impoly(axe,'Closed',false);
else
    out = impoly(axe,pts,'Closed',false);
end


function y=IsValidLabel(x)
% for label deletion
y = ~isempty(x)&&x.isvalid;

function label_btn_new_Callback(hObject, eventdata, handles)
handles.geoLast = handles.geoLast+1;
handles.geos_id = [handles.geos_id handles.geoID];


if exist([handles.DD_out handles.fn 'a_' handles.img_name(1:end-4) '_l.txt'],'file')
    % with elsd preprocession (snap and intersect)
    % default snapping effect
    handles.hand_select = 'snap';
else
    handles.hand_select = 'none';
end
switch handles.geoID
    case 1
        handles.geos{handles.geoLast} = imline(handles.axes_label_im);
    case 2
        %handles.geos{handles.geoLast} = imellipse(handles.axes_label_im);
        % need to change to freehand
        handles.geos{handles.geoLast} = U_imcircle(handles.axes_label_im);
end
set(handles.geos{1},'UserData',handles.geoLast);
if handles.geoID==1
    handles.geos{handles.geoLast}.setColor([1 0 0]);
end
%handles.geos{handles.geoLast}.Deletable = true;
% for new line
cur_pos = handles.geos{handles.geoLast}.getPosition();
if handles.geoID==2
    % check the last pt
    if size(cur_pos,1)>1
        if sum((cur_pos(end,:)-cur_pos(1,:)).^2)<100
            cur_pos(end,:)=cur_pos(1,:);
        end
    end
end

if ~isempty(handles.autogeos)&&~isempty(handles.autogeos{handles.geoID})
    [snap_pos,param] = U_snap(handles.geoID,handles.autogeos{handles.geoID},cur_pos);
    setPosition(handles.geos{handles.geoLast},snap_pos);
    set(handles.geos{handles.geoLast},'UserData',param);
else
    snap_pos = cur_pos;
    % for circle: has to snap ?
    param=[];
end
handles = update_patches(handles,{round(snap_pos(:,[2 1])),param});



% for later line changes
% doesn't work as the button is still being clicked
%addNewPositionCallback(handles.geos{handles.geoLast},@(pos) ...
%    setPosition(handles.geos{handles.geoLast},U_snap(handles.geoID,handles.autogeos{handles.geoID},pos)));
addNewPositionCallback(handles.geos{handles.geoLast},@(pos) ...
    U_label_callback(handles,pos));

handles.geos{handles.geoLast}.setColor([1,0,0]);
guidata(hObject, handles);

function [gid,cur_pt] = U_geoID(handles,pos)
% get the selected geo
tmp_pt = zeros(numel(handles.geos),4);
Isvalid = cellfun(@(x) IsValidLabel(x),handles.geos);
for i=find(Isvalid)
    tmp_pt(i,:) = reshape(handles.geos{i}.getPosition(),1,[]);
end
dist = sum(bsxfun(@minus, tmp_pt,reshape(pos,1,[])).^2,2);
[~,gid] = min(dist);
cur_pt = reshape(tmp_pt(gid,:),[2 2]);

function U_label_callback(handles,pos)
% shortcut key callback
handles = update_handles_callback(handles);
set(handles.figure1,'KeyPressFcn',@(x,y)U_label_key_callback(y,handles,pos));

switch handles.geoID
    case 1
        % in case no key press
        [~,tmp_pos] = U_geoID(handles,pos);
        handles.patches_pt = U_patch_pos(handles,tmp_pos(:,[2 1]));
        update_patches(handles,{round(pos(:,[2 1])),[]});
    case 2
        % bug: no dynamic estimation yet
        %update_patches(handles,pos);
end



function handles = update_handles_callback(handles)
% update parameter
handles.mag_amp = str2double(get(handles.dm_amp,'String'));
handles.label_psz_c = round(str2double(get(handles.label_corner,'String')));
handles.label_psz_s = round(str2double(get(handles.label_sample,'String')));
handles.label_ran_c = -handles.label_psz_c:handles.label_psz_c;
%handles.img = handles.img;
%handles.img_sz = size(handles.img);

function U_label_key_callback(e,handles,pos)
% in the callback function, need to update the parameter
cur_id = U_geoID(handles,pos);
switch e.Key
    case 'q'
        if ~isempty(handles.autogeos{handles.geoID})
            switch handles.hand_select
                case 'snap'
                    [pos,param] = U_snap(handles.geoID,handles.autogeos{handles.geoID},pos);
                    setPosition(handles.geos{handles.geoLast}, pos);
                case 'none'
            end
            % corner display
            update_patches(handles,{round(pos(:,[2 1])),param});
        end
    case {'s','e','d','f','j','i','k','l'}
        % first pt
        delta = zeros(2,2);
        update=-1;
        switch e.Key
            case 's'; delta(1,1) = delta(1,1)-1;update=1;
            case 'f'; delta(1,1) = delta(1,1)+1;update=1;
            case 'e'; delta(1,2) = delta(1,2)-1;update=1;
            case 'd'; delta(1,2) = delta(1,2)+1;update=1;
            case 'j'; delta(2,1) = delta(2,1)-1;update=2;
            case 'l'; delta(2,1) = delta(2,1)+1;update=2;
            case 'i'; delta(2,2) = delta(2,2)-1;update=2;
            case 'k'; delta(2,2) = delta(2,2)+1;update=2;
        end
        %[pos delta]
        %handles.geos{handles.geoLast}
        pos = pos +delta;
        
        % bug: no dynamic update yet
        % second pt
        update_patches(handles,{round(pos(:,[2 1])),[]});
        setPosition(handles.geos{cur_id}, pos);
    case {'z'}
        % delete
        handles.geos{cur_id}.delete();
        undisplay_patches(handles);
    otherwise
end

function [new_pos,param] = U_snap(geoID,autogeos,pos)
[dist, proj_pt,param] = U_dist_geo(geoID,autogeos,pos(:,[2 1])');
%{
figure(2),imshow(zeros([602 900])),hold on, U_plotGeo(zeros(size(autogeos,1),1),autogeos)
%}
%[~,snap_id]= sort(dist,'ascend');
[~,snap_id]= min(dist);
switch geoID
    case 1
        % 1. display the closest
        %pos_snap = handles.autogeos{handles.geoID}(snap_id(1),[2 1 4 3]);
        % 2. interp the drawn pos
        pos_snap = proj_pt([2 1 4 3],snap_id(1));
        %pos_snap = proj_pt([1 2 3 4],snap_id(1));
        new_pos= reshape(pos_snap,[2 2])';
        param=[];
        %disp([pos reshape(pos_snap,[2 2])'])
    case 2
        new_pos = reshape(proj_pt(snap_id,:),2,[])';
        param = [autogeos(snap_id,1:4) param(snap_id,[1 end])];
        % full circle
        if param(5)>=param(end)
            param(end) = param(end)+2*pi;
        end
end

function [dist,out,param] = U_dist_geo(geo_id,geos,geo_ref)
dist = zeros(1,size(geos,1));
param = [];
switch geo_id
    case 1
        out = zeros(4,size(geos,1));
        for i=1:size(geos,1)
            [tmp1,out(1:2,i)] = U_pt2lineseg(geo_ref(1:2),geos(i,:));
            [tmp2,out(3:4,i)] = U_pt2lineseg(geo_ref(3:4),geos(i,:));
            dist(i) = tmp1+tmp2;
        end
    case 2
        % distance to arc
        out = zeros([size(geos,1) 2*size(geo_ref,2)]);
        param = zeros([size(geos,1),size(geo_ref,2)]);
        for i=1:size(geo_ref,2)
            [tmp1,out(:,(i-1)*2+(1:2)),param(:,i)] = U_pt2arcseg(geo_ref(:,i)',geos);
            dist = dist+tmp1;
        end
end
function [dist,pt,theta] = U_pt2arcseg(p,arc)
% arc: x,y,rx,ry,t,th1,th2
dist = zeros(1,size(arc,1));
theta = zeros(size(arc,1),1);
rot_id = find(arc(:,5)==0);
tmp_pos = bsxfun(@rdivide,bsxfun(@minus,p,arc(:,1:2)),arc(rot_id,3:4));
dist(rot_id) = abs(1-sum(tmp_pos(rot_id,:).^2,2)).*arc(rot_id,3);
theta(rot_id) = mod(atan2(tmp_pos(rot_id,2),tmp_pos(rot_id,1)),2*pi);

% need to consider the rotaion of the ellipse
rot_id = find(arc(:,5)~=0);
for rot=rot_id'
    A= [[cos(arc(rot,5)) -sin(arc(rot,5))];[sin(arc(rot,5)) cos(arc(rot,5))]];
    tmp_pos2 = A*tmp_pos(rot,:);
    dist(rot) = abs(1-sum(bsxfun(@rdivide,tmp_pos2,arc(rot,3:4)).^2,2)).*arc(rot,3);
end
theta(rot_id) = mod(atan2(tmp_pos(rot_id,2),tmp_pos(rot_id,1))-arc(rot_id,5),2*pi);
% bug: not consider the finite segment
pt = zeros(size(arc,1),2);

pt(:,2) = arc(:,1) + arc(:,3).*cos(theta).*cos(arc(:,5)) - arc(:,4).*sin(theta).*sin(arc(:,5));
pt(:,1) =  arc(:,2) + arc(:,3).*cos(theta).*sin(arc(:,5)) + arc(:,4).*sin(theta).*cos(arc(:,5));
% figure(2),hold on,plot(pt(:,1),pt(:,2),'bx'),plot(p(1),p(2),'ro')


function [dist,pt] = U_pt2lineseg(p,line)
% given similar orientation, prefer long segments
distance = @(x,y) sum((x-y).^2);
v = line(1:2);
w = line(3:4);
% Return minimum distance between line segment vw and point p
l2 = distance(w,v);
if (l2 == 0)
    dist = distance(w,p);
    pt = w;
else
    % Consider the line extending the segment, parameterized as v + t (w - v).
    % We find projection of point p onto the line.
    % It falls where t = [(p-v) . (w-v)] / |w-v|^2
    t = dot(p - v, w - v) / l2;
    pt = v + t * (w - v);  % Projection may fall out the segment
    if (t < 0.0)
        dist= distance(p, v);       % Beyond the 'v' end of the segment
    elseif (t > 1.0)
        dist= distance(p, w);    % Beyond the 'w' end of the segment
    else
        dist =  distance(p, pt);
    end
end
dist = sqrt(dist);

function label_btn_save_Callback(hObject, eventdata, handles)
pts = zeros(handles.geoLast,5);
Isvalid = cellfun(@(x) IsValidLabel(x),handles.geos);
for i=find(Isvalid)
    tmp_pt = handles.geos{i}.getPosition();
    switch handles.geos_id(i)
        case 1
            % line
            pts(i,:) = [handles.geos_id(i) reshape(tmp_pt(:,[1 2])',1,[])];
        case 2
            % ellipse
            pts(i,:) = [handles.geos_id(i) tmp_pt(1:2) tmp_pt(1:2)+tmp_pt(3:4)-1];
    end
end
pts(Isvalid==0,:) = [];
if ~exist([handles.DD_out handles.fn],'dir')
    mkdir([handles.DD_out handles.fn])
end
pts(:,1)=pts(:,1)-1;
pts(:,2:5)=pts(:,[3 2 5 4])-1;
dlmwrite([handles.DD_out handles.fn handles.DD_out_label],pts)
% reload

%display_img(handles); % remove previous labels
%handles = U_loadLabel(handles);
guidata(hObject, handles);


function handles = U_loadLabel(handles)
if ~exist([handles.DD_out handles.fn],'dir')
    mkdir([handles.DD_out handles.fn])
end
handles.geos =[];
handles.geo_patches =[];
if exist([handles.DD_out handles.DD_out_label],'file')
    % c-index from 0
    pts = load([handles.DD_out handles.DD_out_label]);    
    if numel(pts)>0
        pts(:,2:5) = pts(:,[3 2 5 4]);
        handles.geoLast = size(pts,1);
        handles.geos_id = 1+pts(:,1)';
        handles.geos = cell(1,handles.geoLast);
        for i=1:handles.geoLast
            switch handles.geos_id(i)
                case 1
                    handles.geos{i} = imline(handles.axes_label_im,pts(i,[2 4]),pts(i,[3 5]));
                case 2
                    %handles.geos{i} = imellipse(handles.axes_label_im,[pts(i,2:3),pts(i,4:5)-pts(i,2:3)+1]);
                    handles.geos{i} = U_imcircle(handles.axes_label_im,[pts(i,2:3),pts(i,4:5)]);
            end
            set(handles.geos{1},'UserData',i);
            %handles.geos{handles.geoLast}.Deletable = true;
            handles.geos{i}.setColor([1 0 0]);
            handles.geoID = handles.geos_id(i);
            cur_pos = handles.geos{i}.getPosition();
            
            % only display the patch for the last geo
            switch handles.geoID
                case 1
                    handles = update_patches(handles,{round(cur_pos(:,[2 1])),[]},i==handles.geoLast);
                case 2
                    % bug: need to generate discrete control pt
                    [x,y] = U_samplePts(pts(i,:),1);
                    max_pt =5;
                    if size(x,1)>max_pt
                        ind=unique(round(linspace(1,size(x,1),max_pt)));
                        x=x(ind);y=y(ind);
                    end
                    handles = update_patches(handles,{[],[x y]},i==handles.geoLast);
            end
            
            handles.geo_patches{i} = handles.patches;
        end
        for i=1:handles.geoLast
            addNewPositionCallback(handles.geos{i},@(pos) ...
                U_label_callback(handles,pos));
        end
        handles.dsp_geo = handles.geoLast;
    end
end
handles.autogeos =[];
%[handles.DD_out handles.fn 'a_' handles.img_name(1:end-4) '_l.txt']
if exist([handles.DD_out handles.fn 'l1_' handles.img_name(1:end-3) 'txt'],'file')
    tmp = load([handles.DD_out handles.fn 'l1_' handles.img_name(1:end-3) 'txt']);
    tmp(:,2:5) = tmp(:,[3 2 5 4]);
    arc_len=U_lenGeo(tmp(:,1),tmp(:,2:end));
    tmp(arc_len<handles.snap_thres,:)=[];
    for i=1:3
        % save each kind of geometry separately
        handles.autogeos{i} = tmp(tmp(:,1)==i-1,2:end);
        %{
        if i==1
            % descending sort the endpt of line by x
            [~,pt_id]=sort(handles.autogeos{i}(:,[1 3]),'descend');
            for j=find(pt_id==2)
                handles.autogeos{i}(j,:) = handles.autogeos{i}(j,[3 4 1 2]);
            end
        end
        %}
    end
end


% section 3: magnification parameter
function handles=init_mag(handles)
handles.mag_l = str2double(get(handles.dm_l,'String'));
handles.mag_h = str2double(get(handles.dm_h,'String'));
handles.mag_amp = str2double(get(handles.dm_amp,'String'));
handles.mag_N = str2double(get(handles.dm_N,'String'));
handles.mag_type = Tag2ampType(get(get(handles.panel_ampType,'SelectedObject'),'Tag'));
handles.mag_anti = get(handles.dm_anti,'Value');
handles.mag_pre = get(handles.dm_pre,'Value');
handles.mag_syn = get(handles.dm_syn,'Value');
handles.mag_matt = get(handles.dm_matt,'Value');
handles.raw_edge = [];
handles.filter_edge = [];

%handles.mag_lowwave = str2double(get(handles.dm_lowwave,'String'));
%handles.mag_highwave = str2double(get(handles.dm_highwave,'String'));


function undisplay_mag(handles)
cla(handles.axes_dm_raw),axes(handles.axes_dm_raw),axis off
cla(handles.axes_dm_filter),axes(handles.axes_dm_filter),axis off
cla(handles.axes_dm_syn),axes(handles.axes_dm_syn),axis off



function id = Tag2ampType(str)
id = -1;
switch str
    case 'btn_mag_pos'
        id = 0;
    case 'btn_mag_curve'
        id = 1;
    case 'btn_mag_both'
        id = 2;
end

function display_amp(handles)
cla(handles.axes_dm_syn)
axes(handles.axes_dm_syn);
if ~isempty(handles.img_amp)
    imagesc(handles.img_amp);
end
axis equal,axis off

cla(handles.axes_dm_raw)
axes(handles.axes_dm_raw);
if isfield(handles,'raw_edge')&&~isempty(handles.raw_edge)
    plot(handles.raw_edge);
end
axis tight
cla(handles.axes_dm_filter)
axes(handles.axes_dm_filter);
if isfield(handles,'filter_edge')&&~isempty(handles.filter_edge)
    plot(handles.filter_edge);
end
axis tight

function param=U_handle2param(handles)
param.alpha = handles.mag_amp;
param.N = handles.mag_N; % samples per pixel
param.lambda_l = handles.mag_l;
param.lambda_h = handles.mag_h;
param.preprocess = handles.mag_pre; % perprocessing of the input image.
param.anti_aliasing = handles.mag_anti; % perprocessing of the input image.
param.curve_mag = handles.mag_type;
param.user_syn = handles.mag_syn;
param.method = 'edge';
param.DISPLAYFIGURE = false;
disp([handles.mag_type handles.mag_pre handles.mag_anti])

function btn_ana_Callback(hObject, eventdata, handles)
% analyze the last one

if ~isempty(handles.geos(handles.dsp_geo))
    tmp_pt = handles.geos{handles.dsp_geo}.getPosition();
    param.shape_type = handles.geos_id(handles.dsp_geo);
    switch handles.geos_id(handles.dsp_geo)
        case 1
            % line
            % label
            param.shape_type = 1;
            param.p1 = tmp_pt(1,:); % in case of a circle p1 and p2 are the start and end angles.
            param.p2 = tmp_pt(2,:);
        case 2
            % circle
            %center = (tmp_pt([2 1])+tmp_pt([4 3]))*0.5;
            param.p1 = 0;
            param.p2 = 2*pi;
            tmp_param = get(handles.geos{i},'UserData');
            radius = tmp_param(3);
            center = tmp_param(1:2);
            
            param.shape_type = 2;
            param.shape_param = [radius center];
    end
    [raw_edge, filter_edge] = U_analysis(handles.patches{4},handles.param);
    handles.raw_edge = raw_edge;
    handles.filter_edge = filter_edge;
    display_amp(handles)
    guidata(hObject, handles);
end


function btn_syn_Callback(hObject, eventdata, handles)
Isvalid = cellfun(@(x) IsValidLabel(x),handles.geos);
im_sz = size(handles.img);
warp_field = zeros([im_sz(1:2) 3]);
if im_sz(3)==1
    warp_field(:,:,1) = handles.img;
else
    warp_field(:,:,1) = rgb2gray(handles.img);
end

%param=U_handle2param(handles);
param=handles.param;
raw_edge=[];
filter_edge=[];
for i=find(Isvalid)
    tmp_pt = handles.geos{i}.getPosition();
    param.shape_type = handles.geos_id(i);
    switch  handles.geos_id(i)
        case 1
            % line
            % label
            param.shape_type = 1;
            param.p1 = tmp_pt(1,:); % in case of a circle p1 and p2 are the start and end angles.
            param.p2 = tmp_pt(2,:);
        case 2
            % circle
            %center = (tmp_pt([2 1])+tmp_pt([4 3]))*0.5;
            param.p1 = 0;
            param.p2 = 2*pi;
            tmp_param = get(handles.geos{i},'UserData');
            radius = tmp_param(3);
            center = tmp_param(1:2);
            
            param.shape_type = 2;
            param.shape_param = [radius center];
    end
    %[img_amp, raw_edge, filter_edge, curve_edge] = LineCircMagnifyDev(handles.img, param,img_amp);
    [raw_edge, filter_edge] = U_analysis(handles.geo_patches{i}{4},param);
    warp_field(:,:,2) = warp_field(:,:,2) + U_pointsToGrid(handles.geo_patches{i}{5}(:,1:2), filter_edge(:).*handles.geo_patches{i}{5}(:,3), im_sz(1:2),'nearest');
    warp_field(:,:,3) = warp_field(:,:,3) + U_pointsToGrid(handles.geo_patches{i}{5}(:,1:2), filter_edge(:).*handles.geo_patches{i}{5}(:,4), im_sz(1:2),'nearest');
    %[raw_edge, filter_edge, img_amp] = Mag_gui0120(handles.patches{3}, param,handles.img);
    %img_amp(img_amp>1)=1;img_amp(img_amp<0)=0;
end
% synthesis
% gaussian hack:
%im= handles.img;save Tali_weight im warp_field
%{
figure(2),imagesc(warp_field(:,:,2));
im= handles.img;figure(2),imshow(im);hold on,plot(handles.geo_patches{i}{5}([1 end],1),handles.geo_patches{i}{5}([1 end],2),'b-')
%}
%{
dev_ran = handles.mag_amp;
tmp_f = fspecial('gaussian',30,10);
warp_field = imfilter(warp_field,tmp_f);
warp_field = warp_field/max(abs(warp_field(:)))*dev_ran;
img_amp = U_warpImage(handles.img,warp_field(:,:,1),warp_field(:,:,2));
%}
handles.warp_field = warp_field;
warp_field(:,:,2:3)=warp_field(:,:,2:3)*param.alpha;
[ux, uy] = U_SynFlow(warp_field(:,:,2)~=0,warp_field, param.Dist_th);
handles.img_amp = U_warpImage(handles.img,ux,uy);
handles.raw_edge = raw_edge;
handles.filter_edge = filter_edge;

display_amp(handles)
guidata(hObject, handles);


function btn_saveamp_Callback(hObject, eventdata, handles)
if ~exist([handles.DD_out handles.fn],'dir')
    mkdir([handles.DD_out handles.fn])
end
save_pre = [handles.DD_out handles.fn 'a_' handles.imgs{handles.img_id}(1:end-4)];
num_im = numel(dir([save_pre '_*.png']));
imwrite(handles.img_amp,sprintf([save_pre '_%d.png'],num_im+1))
% save deviation
mag_save.param = handles.param;
mag_save.raw_edge = handles.raw_edge;
mag_save.filter_edge = handles.filter_edge;
save('-v7.3',sprintf([save_pre '_%d.mat'],num_im+1),'mag_save')

function handles = U_loadAmp(handles)
if ~exist([handles.DD_out handles.fn],'dir')
    mkdir([handles.DD_out handles.fn])
end
amp_name = [handles.DD_out handles.fn 'a_' handles.img_name(1:end-4) '_1.png'];
if exist(amp_name,'file')
    handles.img_amp = imread(amp_name);
    display_amp(handles);
end



function labelType_btn_Callback(hObject, eventdata, handles)
handles.geoID = Tag2labelType(get(hObject,'Tag'));
%disp(handles.geoID )
guidata(hObject, handles);

function id = Tag2selectType(str)
id = -1;
switch str
    case 'none'
        id = 0;
    case 'snap'
        id = 1;
    case 'intersect'
        id = 2;
end

function id = Tag2labelType(str)
id = -1;
switch str
    case 'label_btn_line'
        id = 1;
    case 'label_btn_ellipse'
        id = 2;
end

function label_corner_Callback(hObject, eventdata, handles)
handles.label_psz_c = str2double(get(handles.label_corner,'String'));
handles.label_ran_c = -handles.label_psz_c:handles.label_psz_c;
% update the point position (may be updated in callback)
if ~isempty(handles.geos{handles.geoLast})
    tmp_pos = getPosition(handles.geos{handles.geoLast});
    handles.patches_pt = U_patch_pos(handles,tmp_pos(:,[2 1]));
end

handles = update_patches_width(handles,1);
display_patches(handles,1);
guidata(hObject, handles);

function label_sample_Callback(hObject, eventdata, handles)
handles.label_psz_s = str2double(get(handles.label_sample,'String'));
% update the point position (may be updated in callback)
if ~isempty(handles.geos{handles.geoLast})
    tmp_pos = getPosition(handles.geos{handles.geoLast});
    switch handles.geoID
        case 1
            handles.patches_pt = U_patch_pos(handles,tmp_pos(:,[2 1]));
        case 2
            % center radius theta
            handles.patches_pt = U_pt2circ(tmp_pos(:,[2 1]));
    end
end
handles=update_patches_width(handles,2);
display_patches(handles,2);
guidata(hObject, handles);




function dm_matt_Callback(hObject, eventdata, handles)
handles.mag_matt = get(handles.dm_matt,'Value');
guidata(hObject, handles);



function dm_anti_Callback(hObject, eventdata, handles)
handles.mag_anti = get(handles.dm_anti,'Value');
guidata(hObject, handles);

function dm_pre_Callback(hObject, eventdata, handles)
handles.mag_pre = get(handles.dm_pre,'Value');
guidata(hObject, handles);

function dm_l_Callback(hObject, eventdata, handles)
handles.mag_l = str2double(get(handles.dm_l,'String'));
guidata(hObject, handles);

function dm_h_Callback(hObject, eventdata, handles)
handles.mag_h = str2double(get(handles.dm_h,'String'));
guidata(hObject, handles);

function dm_N_Callback(hObject, eventdata, handles)
handles.mag_N = str2double(get(handles.dm_N,'String'));
guidata(hObject, handles);

function dm_amp_Callback(hObject, eventdata, handles)
handles.mag_amp = str2double(get(handles.dm_amp,'String'));

if ~isempty(handles.warp_field)
    % redo warping
    warp_field=handles.warp_field;
    warp_field(:,:,2:3)=handles.warp_field(:,:,2:3)*handles.mag_amp;
    [ux, uy] = U_SynFlow(warp_field(:,:,2)~=0,warp_field, handles.param.Dist_th);
    handles.img_amp = U_warpImage(handles.img,ux,uy);
    display_amp(handles)
end

%handles.mag_type = Tag2ampType(get(get(handles.panel_ampType,'SelectedObject'),'Tag'));
guidata(hObject, handles);

function dm_ampType_Callback(hObject, eventdata, handles)
handles.mag_type = Tag2ampType(get(hObject,'Tag'));
%handles.mag_type
guidata(hObject, handles);

function dm_syn_Callback(hObject, eventdata, handles)
handles.mag_syn = get(handles.dm_syn,'Value');
guidata(hObject, handles);















function dm_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dm_amp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dm_l_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dm_N_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dm_lowwave_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dm_highwave_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function label_corner_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function label_sample_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sel_id_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end