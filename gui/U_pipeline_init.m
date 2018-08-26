function handles = P_pipeline_init(handles,opt,varargin)
% constraint 1: smaller gui file
% constraint 2: fewer .m files
% constraint 3: matlab only search the first function in .m file
global geoKey cur_geo_id
switch opt
    case 'init_const'
        handles.geoC = {[1,0,0],[0,1,0],[0,0,1]};
        handles.Ns = [2,3,5];
        tmp_type = get(handles.panel_labelType,'SelectedObject');
        handles.geoID = U_pipeline_util('tag_label',get(tmp_type,'Tag'));
        handles.geo_vertex = cell(1,2);
        % init display patch size
        handles.label_psz_c = str2double(get(handles.label_corner,'String'));
        handles.label_psz_s = str2double(get(handles.label_height,'String'));
        handles.label_ran_c = -handles.label_psz_c:handles.label_psz_c;
    case 'init_img'
        handles.img_id = 1;
        handles.img_label_id = 1;
        handles.img_amp_id = 1;
        geoKey = '';
        cur_geo_id = [0,0];% geo id and control pt id                
        handles = U_pipeline_update(handles,'img');
    case 'init_mag'
        param =varargin{1};
        handles.mag_l = param(1);
        handles.mag_h = param(2);
        handles.mag_amp = param(3);
        handles.mag_N = param(4);
        handles.mag_type = param(5);
        handles.mag_anti = param(6);
        handles.mag_pre = param(7);
        handles.mag_syn = param(8);
        handles.mag_syn_box = param(9:10);
        handles.mag_matt = param(11);
        handles.raw_edges = [];
        handles.filter_edges = [];
    case 'init_label'
        % init geometry
        handles.geos = {};
        handles.geo_ids = [];
        handles.geo_pts = {};
        handles.geo_patches = {};
        
    case 'param_save' % for save
        handles =[handles.mag_l handles.mag_h handles.mag_amp handles.mag_N ...
            handles.mag_type handles.mag_anti handles.mag_pre ...
            handles.mag_syn handles.mag_syn_box handles.mag_matt];
    case 'param_dm'
        param.alpha = handles.mag_amp;
        param.N = handles.mag_N; % samples per pixel
        param.lambda_l = handles.mag_l;
        param.lambda_h = handles.mag_h;
        param.preprocess = handles.mag_pre; % perprocessing of the input image.
        param.anti_aliasing = handles.mag_anti; % perprocessing of the input image.
        param.curve_mag = handles.mag_type;
        param.user_syn = handles.mag_syn;
        param.shape_type = 1;
        param.Dist_th = 20;
        handles = param;
    case 'param_dsp'
        param = varargin{1};
        % magnification panel
        set(handles.dm_l,'String',sprintf('%d',param(1)));
        set(handles.dm_h,'String',sprintf('%d',param(2)));
        set(handles.dm_amp,'String',sprintf('%d',param(3)));
        set(handles.label_N,'String',sprintf('%d',param(4)));
        %set(handles.panel_ampType,'SelectedObject');
        set(handles.dm_anti,'Value',param(6));
        set(handles.dm_pre,'Value',param(7));
        set(handles.dm_syn,'Value',param(8));
        set(handles.label_matt,'Value',param(11));
end
