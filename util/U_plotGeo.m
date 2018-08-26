function U_plotGeo(geo_type,geo_params,linwidth,cc)
if ~exist('linwidth','var');linwidth = 3;end

% geo_param output from elsd
for i=1:numel(geo_type)
    geo_param = geo_params(i,:);
    if exist('cc','var') 
        if size(cc,1)>i
            cc_s= cc(i,:);
        else
            cc_s= cc(1,:);
        end
    end
    switch geo_type(i)
        case 0
            % line
            if ~exist('cc','var')
                plot(geo_param([2 4]),geo_param([1 3]),'r-','LineWidth',linwidth);
            else
                
                plot(geo_param([2 4]),geo_param([1 3]),'color',cc_s,'LineWidth',linwidth);
            end
            %line(geo_param([1 3]),geo_param([2 4]),'g');
        case 1
            % circle
            %theta =  (geo_param(end-1):(pi/180):geo_param(end));
            %xx = geo_param(1) + geo_param(3) * sin(theta);
            %yy = geo_param(2) + geo_param(4) * cos(theta);
            [xx,yy] = U_arc(geo_param(1),geo_param(2),geo_param(3),geo_param(4),0,100,geo_param(6:7));
            if ~exist('cc','var')
                plot( yy, xx,'g-','LineWidth',linwidth)
            else
                plot( yy, xx,'color',cc_s,'LineWidth',linwidth)
            end
        case 2
            % ellipse
            [xx,yy] = U_arc(geo_param(1),geo_param(2),geo_param(3),geo_param(4),geo_param(5),100,geo_param(6:7));
            if ~exist('cc','var')
                plot( yy, xx, 'b-','LineWidth',linwidth)
            else
                plot( yy, xx, 'color',cc_s,'LineWidth',linwidth)
            end
    end
end

