function handles=U_pipeline_main(handles,opt,varargin)
global geoKey cur_geo_id param_default
switch opt
    case 'ana'
        if numel(handles.geos)>=cur_geo_id(1) && ~isempty(handles.geo_pts{cur_geo_id(1)}{1})
            param = U_pipeline_init(handles,'param_dm');
            if handles.geo_pts{cur_geo_id(1)}{1}(8)-handles.geo_pts{cur_geo_id(1)}{1}(7)>2*pi-0.01
                param.shape_type = 2;
            end
            [handles.raw_edges{cur_geo_id(1)}, handles.filter_edges{cur_geo_id(1)}] = U_analysis(handles.geo_patches{cur_geo_id(1)}{4},param);
            handles = U_pipeline_display(handles,'amp');
        end
    case 'syn'
        Isvalid = U_pipeline_util('geo_valid',handles.geos);
        im_sz = handles.img_sz;
        warp_field = zeros([im_sz(1:2) 3]);
        if im_sz(3)==1
            warp_field(:,:,1) = handles.img;
        else
            warp_field(:,:,1) = rgb2gray(handles.img);
        end
        param = U_pipeline_init(handles,'param_dm');
        for i=find(Isvalid)
            param.shape_type = 1;
            if handles.geo_pts{i}{1}(8)-handles.geo_pts{i}{1}(7)>2*pi-0.01
                param.shape_type = 2;
            end
            %[img_amp, raw_edge, filter_edge, curve_edge] = LineCircMagnifyDev(handles.img, param,img_amp);
            [handles.raw_edges{i}, handles.filter_edges{i}] = U_analysis(handles.geo_patches{i}{4},param);
            warp_field(:,:,2) = warp_field(:,:,2) + U_pointsToGrid(handles.geo_patches{i}{5}(:,1:2), handles.filter_edges{i}(:).*handles.geo_patches{i}{5}(:,3), im_sz(1:2),'nearest');
            warp_field(:,:,3) = warp_field(:,:,3) + U_pointsToGrid(handles.geo_patches{i}{5}(:,1:2), handles.filter_edges{i}(:).*handles.geo_patches{i}{5}(:,4), im_sz(1:2),'nearest');
        end
        handles.warp_field = warp_field;
        warp_field(:,:,2:3)=warp_field(:,:,2:3)*param.alpha;
        [ux, uy] = U_SynFlow(warp_field(:,:,2)~=0,warp_field, param.Dist_th);
        handles.img_amp = U_warpImage(handles.img,ux,uy);
        handles = U_pipeline_display(handles,'amp');
    case 'label_new'
        new_id = numel(handles.geo_ids)+1;
        handles.geo_ids = [handles.geo_ids handles.geoID];
        if ~isempty(handles.geos)
            % with elsd preprocession (snap and intersect)
            % default snapping effect
            handles.hand_select = 'snap';
        else
            handles.hand_select = 'none';
        end
        handles.geos{new_id} = imfreehand(handles.axes_label_im,'Closed',false);
        
        % for new line
        tmp_pos = handles.geos{new_id}.getPosition();
        % check closure:
        if size(tmp_pos,1)>1
            if sum((tmp_pos(end,:)-tmp_pos(1,:)).^2)<100
                tmp_pos(end,:)=tmp_pos(1,:);
            end
        end
        
        no_snap=0;
        snap_thres = (max(handles.img_sz)*0.02)^2;
        if ~isempty(handles.autogeos)&&~isempty(handles.autogeos{1+handles.geoID})
            [snap_pos,param] = U_snap(handles.geoID,handles.autogeos{1+handles.geoID},tmp_pos,snap_thres);
            param = [param new_id];
            if isempty(snap_pos)
                no_snap=1;
            end
        end
        if no_snap
            % no snap
            snap_pos = tmp_pos;
            param =[U_geofit(tmp_pos(:,1),tmp_pos(:,2),handles.geoID) new_id];
        end
        handles.geos{new_id}.delete;
        switch handles.geoID
            case 0
                handles.geos{new_id} = imline(handles.axes_label_im,snap_pos(:,2),snap_pos(:,1));
            case {1,2}
                handles.geos{new_id} = U_imcircle(handles.axes_label_im,snap_pos(:,[2 1]),param);
        end
        
        %cur_pos = snap_pos(:,[2 1]);
        if cur_geo_id(1)>0 && ~isempty(handles.geos{cur_geo_id(1)}) && handles.geos{cur_geo_id(1)}.isvalid
            handles.geos{cur_geo_id(1)}.setColor(handles.geoC{3});
        end
        
        handles.geos{new_id}.setColor(handles.geoC{1});
        handles.geo_pts{new_id}{1} = param;
        handles.geo_pts{new_id}{2} = snap_pos;
        handles.geo_patches{new_id} = cell(1,5);
        cur_geo_id = [new_id 2];
        handles = U_pipeline_update(handles,'label_patch',3,1:2,1);
        
    case 'label_drag'
        pos = varargin{1};
        if numel(geoKey)==0
            % mouse dragging
            % recalculate geo and control pt
            [cur_id,pid] = U_geoFind(handles,pos);
            geo_type = handles.geo_ids(cur_id);
            % fit new parameter/update patch
            param = U_geofit(pos(:,2),pos(:,1),geo_type);
            handles.geo_pts{cur_id}{1} = param;
            handles.geo_pts{cur_id}{2} = pos;
            setPosition(handles.geos{cur_id}, pos);
            cur_geo_id = [cur_id,pid];
        else
            % keyboard change
            geoKey ='';
        end
    case 'label_elsd'
        tmp_im = imresize(im2double(snapshot(handles.cam)),param_default(5));
        elsd_mex_all(255*tmp_im,param_default(1:4),handles.elsd_file);
        handles = U_video_io(handles,'label_load_elsd');
        % correction for resolution
        for i=1:3
            handles.autogeos{1}{1}(:,2:4) = handles.autogeos{1}{1}(:,2:4)/param_default(5);
            if i>1
                handles.autogeos{i}{2} = handles.autogeos{i}{2}/param_default(5);
            end
        end
    case 'label_key'
        if ~isempty(handles.geos{cur_geo_id(1)})
            geo_type = handles.geo_ids(cur_geo_id(1));
            pos = handles.geo_pts{cur_geo_id(1)}{2};
            switch geoKey
                case 'q'
                    if ~isempty(handles.autogeos{geo_type+1})
                        [pos,param]= U_snap(handles.geoID,handles.autogeos{geo_type+1},pos);
                        handles.geo_pts{cur_geo_id(1)}{1}(1:end-1) = param;
                        handles.geo_pts{cur_geo_id(1)}{2} = pos;
                        
                    end
                    % corner display
                    handles=U_pipeline_update(handles,'label_patch',3,1:2,1);
                    setPosition(handles.geos{cur_geo_id(1)}, pos);
                case 'p'
                    % extension
                    if ~isempty(handles.autogeos{geo_type+1})
                        [pos,param]= U_snap(handles.geoID,handles.autogeos{geo_type+1},pos);
                        handles.geo_pts{cur_geo_id(1)}{1}(1:end-1) = param;
                        handles.geo_pts{cur_geo_id(1)}{2} = pos;
                        
                    end
                    % corner display
                    handles=U_pipeline_update(handles,'label_patch',3,1:2,1);
                    setPosition(handles.geos{cur_geo_id(1)}, pos);
                case {'s','e','d','f','j','i','k','l'}
                    % first pt
                    delta = zeros(2,2);
                    up_pind = [1 cur_geo_id(2)];
                    switch geoKey
                        % global motion
                        case 's'; delta(:,2) = delta(:,2)-1;
                        case 'f'; delta(:,2) = delta(:,2)+1;
                        case 'e'; delta(:,1) = delta(:,1)-1;
                        case 'd'; delta(:,1) = delta(:,1)+1;
                            % local motion
                        case 'j'; delta(2,2) = delta(2,2)-1;up_pind(1)=[];
                        case 'l'; delta(2,2) = delta(2,2)+1;up_pind(1)=[];
                        case 'i'; delta(2,1) = delta(2,1)-1;up_pind(1)=[];
                        case 'k'; delta(2,1) = delta(2,1)+1;up_pind(1)=[];
                    end
                    pos([1 cur_geo_id(2)],:) = pos([1 cur_geo_id(2)],:) + delta;
                    
                    % dynamic update
                    param = U_geofit(pos(:,1),pos(:,2),geo_type);
                    handles.geo_pts{cur_geo_id(1)}{1}(1:end-1) = param;
                    handles.geo_pts{cur_geo_id(1)}{2} = pos;
                    setPosition(handles.geos{cur_geo_id(1)}, pos(:,[2 1]));
                    handles=U_pipeline_update(handles,'label_patch',3,up_pind,1);
                case {'z'}
                    % delete
                    handles.geos{cur_geo_id}.delete();
                    handles.geo_pts{cur_geo_id} = [];
                    U_pipeline_display(handles,'unpatch')
                otherwise
            end
        end
end