function handles=U_pipeline_update(handles,opt,varargin)
global cur_pos cur_geo_id
% callback functions
switch opt
    case 'img'
        set(handles.img_sel,'String',num2str(handles.img_id));
        set(handles.img_num,'String',sprintf('(     / %d)',numel(handles.img_names)));
        set(handles.img_name,'String',handles.img_names{handles.img_id});
        set(handles.img_label,'String',sprintf('%d',handles.img_label_id));
        set(handles.img_amp,'String',sprintf('%d',handles.img_amp_id));
        
        handles.img_amp=[];
        handles.img = im2double(imread([handles.DD_in handles.img_names{handles.img_id}]));
        
        im_name = handles.img_names{handles.img_id}(1:end-4);
        handles.elsd_file = [handles.DD_out 'pre1_' im_name '_pt2.txt'];
        handles.label_file = sprintf([handles.DD_out 'pre1_' im_name '_l%d.txt'],handles.img_label_id);
        handles.amp_file = sprintf([handles.DD_out 'pre1_' im_name '_l%d_a%d.mat'],handles.img_label_id,handles.img_amp_id);
        handles = U_pipeline_io(handles,'dm_load');
        handles = U_pipeline_io(handles,'label_load_elsd');
        handles=U_pipeline_display(handles,'img');
        handles=U_pipeline_display(handles,'unpatch');
        handles=U_pipeline_display(handles,'unmag');
        handles = U_pipeline_io(handles,'label_load_manual');        
    case 'img_label'
        handles.img_amp_id = str2double(get(handles.img_label,'String'));
        label_name = sprintf([handles.elsd_file(1:end-6) '_l%d.txt'],handles.img_label_id);
        if exist(label_name,'file')
            
        end
    case 'img_amp'
        handles.img_amp_id = str2double(get(handles.img_amp,'String'));
        amp_name = sprintf([handles.elsd_file(1:end-6) '_l%d_a%d.txt'],handles.img_label_id,handles.img_amp_id);
        if exist(amp_name,'file')
            
        end
    case 'img_next'
        if handles.img_id == numel(handles.img_names)
            handles.img_id= 1;
        else
            handles.img_id =handles.img_id+1;
        end
        handles=U_pipeline_update(handles,'img');
    case 'img_prev'
        if handles.img_id == 1
            handles.img_id = numel(handles.img_names);
        else
            handles.img_id = handles.img_id-1;
        end
        handles=U_pipeline_update(handles,'img');
    case 'img_sel'
        new_id = min(numel(handles.img_names),max(1,round(str2double(get(handles.sel_id,'String')))));
        if new_id~= handles.img_id
            handles.img_id = new_id;
            handles=U_pipeline_update(handles,'img');
        end
    case 'label_corner'
        handles.label_psz_c = str2double(get(handles.label_corner,'String'));
        handles.label_ran_c = -handles.label_psz_c:handles.label_psz_c;
        % update the point position (may be updated in callback)
        if ~isempty(handles.geos{cur_geo_id(1)})
            handles = U_pipeline_update(handles,'label_corner_patch',1,1:2);
            U_pipeline_display(handles,'patch',1)
        end
    case 'label_height'
        handles.label_psz_s = str2double(get(handles.label_height,'String'));
        % update the point position (may be updated in callback)
        if ~isempty(handles.geos{cur_geo_id(1)})
            handles=U_pipeline_update(handles,'patch_width',2);
            U_pipeline_display(handles,'patch',2)
        end
    case 'label_N'
        handles.mag_N = str2double(get(handles.label_N,'String'));
        if ~isempty(handles.geos{cur_geo_id(1)})
            handles=U_pipeline_update(handles,'label_patch');
            % update analysis
            handles=U_pipeline_main(handles,'ana');
        end
    case 'label_patch'
        opt_w = varargin{1};
        patch_id = varargin{2};
        dsp = varargin{3};
        
        if(opt_w==1 || opt_w==3)
            % update the stripe corner            
            pts = handles.geo_pts{cur_geo_id(1)}{2}(patch_id,:);
            pts = U_pipeline_util('pos_thres',pts,handles.img_sz);            
            for i=1:numel(patch_id)
                pid = i+2-numel(patch_id);
                if (pts(i,1)>handles.label_psz_c ...
                        && pts(i,1)<=handles.img_sz(1)-handles.label_psz_c ...
                        && pts(i,2)>handles.label_psz_c ...
                        && pts(i,2)<=handles.img_sz(2)-handles.label_psz_c)
                    handles.geo_patches{cur_geo_id(1)}{pid} = handles.img( ...
                        pts(i,1)+handles.label_ran_c, ...
                        pts(i,2)+handles.label_ran_c,:);
                else
                    handles.geo_patches{cur_geo_id(1)}{pid} = zeros(handles.label_psz_c*2+1,handles.label_psz_c*2+1,3);
                    st_x = max(1,pts(i,1)-handles.label_psz_c);
                    lt_x = min(handles.img_sz(1),pts(i,1)+handles.label_psz_c);
                    st_y = max(1,pts(i,2)-handles.label_psz_c);
                    lt_y = min(handles.img_sz(2),pts(i,2)+handles.label_psz_c);
                    
                    pst_x = handles.label_psz_c-(pts(i,1)-st_x)+1;
                    plt_x = handles.label_psz_c+(lt_x-pts(i,1))+1;
                    pst_y = handles.label_psz_c-(pts(i,2)-st_y)+1;
                    plt_y = handles.label_psz_c+(lt_y-pts(i,2))+1;
                    handles.geo_patches{cur_geo_id(1)}{pid}(pst_x:plt_x,pst_y:plt_y,:) = handles.img(st_x:lt_x,st_y:lt_y,:);
                end
            end
            if patch_id(end)==1
                % special case: need to update both patches
                handles.geo_patches{cur_geo_id(1)}{1} =  handles.geo_patches{cur_geo_id(1)}{2};
            end
        end
        if(opt_w==2 || opt_w==3)
            % update the stripe height
            tmp_param = handles.geo_pts{cur_geo_id(1)}{1};            
            num_pt = U_lenGeo(tmp_param(1),tmp_param(2:end))*handles.mag_N;
            tmp_pts = U_shape2pts(tmp_param(1),tmp_param(2:end),num_pt);
            [handles.geo_patches{cur_geo_id(1)}{3},handles.geo_patches{cur_geo_id(1)}{5}] = U_SampleStripe(tmp_pts{1}(:,2)', tmp_pts{1}(:,1)', handles.img, handles.label_psz_s,handles.mag_N);
            % figure(1),clf,imagesc(handles.patches{3})
            if handles.mag_matt
                offset = 1e-3;
                Mask = handles.geo_patches{cur_geo_id(1)}{3}+offset;
                st=1;lt=3;
                %lt=ceil(handles.label_psz_s/step);
                %lt=ceil(2*handles.label_psz_s/step);
                Mask(st:lt,:,:) = 1;
                Mask(end-(lt:-1:st)+1,:,:) = 0;
                %Mask =  CreatMaskLine_0427(handles.patches{3}, param.p1, param.p2);
                handles.geo_patches{cur_geo_id(1)}{4} =runMattingfun(handles.geo_patches{cur_geo_id(1)}{3}+offset,Mask);
            end
            % figure(1),clf,imagesc(handles.patches{4})
        end
        if dsp
            handles=U_pipeline_display(handles,'patch',opt_w);
        end
    case 'label_geo'
        pos = varargin{1};
        up_pind = varargin{2};        
        error('undone')
        switch handles.geoID
            case {2,3}
                pos0=pos;
                handles.geo_pts{cur_geo_id(1)}{1} = U_geofit(pos(:,2),pos(:,1),handles.geoID);
                pos = U_shape2pts(handles.geoID,tmp_param,handles.Ns(handles.geoID));
                pos = pos{1}(:,[2 1]);
                if norm(pos(1,:)- pos0(1,:))>norm(pos(1,:)- pos0(end,:))
                    pos =pos(end:-1:1,:);
                end
                if up_pind(1)~=1;up_pind=[1 up_pind];end
            case 4
                tmp_param = cell(size(pos));
                for i=1:numel(pos)
                    tmp_param{i} = U_geofit(pos{i}(:,1),pos{i}(:,2),handles.geoID);
                end
        end
        % update patch
        handles=U_pipeline_update(handles,'patch',{pos(:,[2 1]),[handles.geoID tmp_param]},3,up_pind);
    case 'dm_next_pt'
        cur_geo_id(2) = cur_geo_id(2)+1;
        if  cur_geo_id(2) > handles.Ns(handles.geo_ids(cur_geo_id(1))+1)
            cur_geo_id(2) = 1;
        end
        % update corner 
        handles = U_pipeline_update(handles,'label_patch',1,cur_geo_id(2),0);
        handles=U_pipeline_display(handles,'patch_vertex');
    case 'dm_next_shape'
        if handles.geos{cur_geo_id(1)}.isvalid
            handles.geos{cur_geo_id(1)}.setColor(handles.geoC{3});
        end
        ind = U_pipeline_util('geo_valid',handles.geos);
        if numel(ind)~=0
            if cur_geo_id(1) >= ind(end)
                cur_geo_id(1) = ind(1);
            else
                next_ind = ind(ind>cur_geo_id(1));
                cur_geo_id(1) =next_ind(1);
            end           
            handles.geos{cur_geo_id(1)}.setColor(handles.geoC{1});        
            cur_geo_id(2) = 2;            
            handles=U_pipeline_display(handles,'patch',3);
        end
    case 'dm_amp'
        handles.mag_amp = str2double(get(handles.dm_amp,'String'));
        % update synthesis
        if ~isempty(handles.warp_field)
            warp_field=handles.warp_field;
            warp_field(:,:,2:3)=handles.warp_field(:,:,2:3)*handles.mag_amp;
            [ux, uy] = U_SynFlow(warp_field(:,:,2)~=0,warp_field, handles.param.Dist_th);
            handles.img_amp = U_warpImage(handles.img,ux,uy);
            U_pipeline_display(handles,'amp')
        end
        %handles.mag_type = Tag2ampType(get(get(handles.panel_ampType,'SelectedObject'),'Tag'));
    case 'dm_matt'
        tmp_matt =  get(handles.label_matt,'Value');
        if tmp_matt ~= handles.mag_matt
            handles.mag_matt= tmp_matt;
            switch tmp_matt
                case 0
                    % display raw
                case 1
            end
        end
end