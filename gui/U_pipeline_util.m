function out = U_pipeline_util(opt,varargin)
switch opt
    case 'geo_valid'
        out = find(cellfun(@(x) ~isempty(x)&&x.isvalid ,varargin{1}));    
    case 'pos_thres'
        pos = varargin{1};
        im_sz = varargin{2};
        out = round(pos);
        out = max(out,1);
        out(:,1) = min(out(:,1),im_sz(1));
        out(:,2) = min(out(:,2),im_sz(2));
    case 'tag_snap'
        out = -1;
        switch varargin{1}
            case 'none';out = 0;
            case 'snap';out = 1;
            case 'intersect';out = 2;
        end
    case 'tag_amp'        
        switch varargin{1}
            case 'btn_mag_pos';out = 0;
            case 'btn_mag_curve';out = 1;
            case 'btn_mag_both';out = 2;
        end
    case 'tag_label'
        out = -1;
        switch varargin{1}
            case 'label_btn_line';out = 0;
            case 'label_btn_circle';out = 1;
            case 'label_btn_ellipse';out = 2;
            case 'label_btn_best';out = 3;
        end
end
