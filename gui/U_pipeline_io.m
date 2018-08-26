function handles = U_pipeline_io(handles,opt,varargin)
global cur_geo_id param_p_default
% constraint 1: smaller gui file
% constraint 2: fewer .m files
% constraint 3: matlab only search the first function in .m file
switch opt
    case 'label_load_manual'
        if exist(handles.label_file,'file')
            param = load(handles.label_file);
            if numel(param)>0
                % saved label from 1
                num_geo = size(param,1);
                handles.geo_ids = param(:,1)';
                handles.geos = cell(1,num_geo);
                handles.geo_pts = cell(1,num_geo);
                for i=1:num_geo
                    switch handles.geo_ids(i)
                        case 0
                            handles.geos{i} = imline(handles.axes_label_im,param(i,[3 5]),param(i,[2 4]));
                            %tmp_pts = reshape(param(i,[3 2 5 4]),[2,2])';
                            tmp_pts = reshape(param(i,2:5),[2,2])';
                        case {1,2}
                            % at least five points to estimate circle/ellipse
                            %tmp_pts = U_shape2pts(handles.geo_ids(i),param(i,[3 2 5 4 6:end]),handles.Ns(1+handles.geo_ids(i)));
                            tmp_pts = U_shape2pts(handles.geo_ids(i),param(i,[2:end]),handles.Ns(1+handles.geo_ids(i)));
                            tmp_pts = tmp_pts{1};
                            handles.geos{i} = U_imcircle(handles.axes_label_im,tmp_pts(:,[2 1]),param(i,:));
                    end
                    handles.geos{i}.setColor(handles.geoC{3-2*(i==num_geo)});
                    
                    handles.geo_pts{i} = cell(1,2);
                    handles.geo_pts{i}{1} = param(i,:);
                    %handles.geo_pts{i}{1}(2:5) = param(i,[3 2 5 4]);
                    handles.geo_pts{i}{2} = tmp_pts;
                    cur_geo_id = [i 2];
                    % only display the patch for the last geo
                    handles = U_pipeline_update(handles,'label_patch',3,1:2,i==num_geo);
                end
            end
        else
            handles=U_pipeline_init(handles,'init_label',param_p_default);
        end
    case 'label_load_elsd'
        if ~exist(handles.elsd_file,'file')
            gui_label_nowindow({handles.DD_in,handles.DD_out,{handles.img_names{handles.img_id}}});
        end
        elsd_param = load([handles.elsd_file(1:end-6) '.txt']);
        if elsd_param(5)~=1
            handles.img = imresize(handles.img,elsd_param(5));
            handles.img(handles.img>1) = 1;
            handles.img(handles.img<0) = 0;
        end
        handles.img_sz = size(handles.img);
        
        handles.autogeos =[];
        tmp = load(handles.elsd_file);
        %tmp(:,2:5) = tmp(:,[3 2 5 4]);
        for i=1:3
            % save each kind of geometry separately
            % c: index from 0
            if i==1
                handles.autogeos{i} = cell(1);
                handles.autogeos{i}{1} = tmp(tmp(:,1)==i-1,2:end);
                handles.autogeos{i}{1}(:,2:5) = handles.autogeos{i}{1}(:,2:5)+1;
            else
                % param and pts representation
                handles.autogeos{i} = cell(1,2);
                lid = find(tmp(:,1)==i-1);
                handles.autogeos{i}{1} = tmp(lid,2:end);
                % shape center
                handles.autogeos{i}{1}(:,1:2) = handles.autogeos{i}{1}(:,1:2)+1;
                % id=2;V_label(handles.img,handles.autogeos{i}{1}(id,:))
                %aa = U_shape2pts(i*ones(numel(lid),1),handles.autogeos{i}{1}(:,1:end),N(i));
                tmp_pt = U_shape2pts((i-1)*ones(numel(lid),1),handles.autogeos{i}{1}(:,1:end),handles.Ns(i),handles.img_sz(1:2));
                bid = [arrayfun(@(x) size(tmp_pt{x},1)~=handles.Ns(i),1:numel(tmp_pt))];
                tmp_pt(bid)=[];
                % param representation
                handles.autogeos{i}{1}(bid,:)=[];
                % pts representation
                handles.autogeos{i}{2} = reshape(cell2mat(tmp_pt)',2*handles.Ns(i),[])';
            end
        end
        % reference for snapping
        % figure(1),imshow(handles.img);hold on,U_plotGeo(tmp(:,1),tmp(:,2:end));
        % title('Pre-processed shapes for snapping')
        
    case 'label_save'
        pts = ones(numel(handles.geo_pts),10);
        Isvalid = U_pipeline_util('geo_valid',handles.geos);
        for i=find(Isvalid)
            pts(i,:) = handles.geo_pts{i}{1};
        end
        pts(Isvalid==0,:) = [];
        % saved label from 1
        dlmwrite(handles.label_file,pts)
    case 'dm_save'
        mag_save=[];
        % save parameter
        
        % save image
        mag_save.img_amp = handles.img_amp;
        % save deviation
        mag_save.raw_edge = handles.raw_edge;
        mag_save.filter_edge = handles.filter_edge;
        save('-v7.3',handles.amp_file,'mag_save')
    case 'dm_load'
        % load amp parameters
        if exist(handles.amp_file,'file')
            load(handles.amp_file);
            % display amp parameter
            handles=U_pipeline_init(handles,'init_mag',mag_save.param);
            % display amp image
            handles.img_amp = mag_save.img_amp;
            U_pipeline_display(handles,'amp');
            % store analysis result
        else
            handles=U_pipeline_init(handles,'init_mag',param_p_default);
            handles=U_pipeline_init(handles,'param_dsp',param_p_default);
        end
        
end
