%% Step 1. Read Image
Dir_in = 'data/';
Dir_out = 'result/';
im_name='sand.jpg';


im_name_p = im_name(1:find(im_name=='.',1,'last')-1);
out_name = ['a1_' im_name_p]; % save name
im = im2double(imread([Dir_in im_name]));
% possible resize
im =imresize(im,im_ratio);
im(im>1)=1;im(im<0)=0;
im_sz = size(im);

% need to convert image to pgm if needed
% imwrite(imread([Dir_in im_name]),[Dir_in im_name_p '.pgm'])

%% Step 2. setup parameters
% 1.1 rough line detection
% default for elsd parameters
thres=[22.5,0.7,1e0,1];
im_ratio=1;
snap_thres=5;
% 1.2 sub-pixel accurate line estimation
param.N = 2;% sampled points per pixel
param.delta =10; %half width of the stripe
param.matt_offset = 1e-3;
param.matt_ran = ceil([2 (param.delta/3)]);
param.matt_high = 1e-3;
% 1.3 deformation analysis
param.lambda_l = 10*param.N;
param.lambda_h = 300*param.N;
param.curve_mag = 0;
param.alpha = 10;
param.method = 'edge';
param.contour_close = 5;
% 1.4 amplification
param.Dist_th = 50;

%% Step 3. Find Shapes (line/ellipse)
%% 3.1 automatic shapes (elsd)
if ~exist([Dir_out 'l1_' im_name_p '.txt'],'file')
    % all possible shapes without voting comparsion
    % useful for finding all instances for one specific shape
    % e.g. manual snapping
    sn = [Dir_out 'l1_' im_name_p '.txt'];
    % option 1: binary file
    cmd_p = sprintf(' %.2f %.2f %.2f %d %s',thres(1),thres(2),thres(3),thres(4),sn); 
    system(['lib/elsd/elsd_io ' Dir_in im_name_p '.pgm' cmd_p]);
    % option 2: mex file
    %elsd_mex_all(255*im,thres,sn);
end

%% 3.2 label shapes
if ~exist([Dir_out out_name '.txt'],'file')
    % manual
    gui_pipeline([],1,{[Dir_in im_name],im_ratio,Dir_out,[out_name '.txt'],snap_thres},param);
else
    % each cell is a contour
    % each contour is matrix of parameters
    % assume one contour
    tmp = load([Dir_out out_name '.txt']);
    m_shapes = {tmp};
    
end

%% 3.3 convert parametric shapes to points
m_pts = cell(1,numel(m_shapes));
%m_shape_type =[];
for i=1:numel(m_shapes)
    tmp_pt =[];
    for j = size(m_shapes{i},1)
        % for each shape
        [x,y] = U_samplePts(m_shapes{i}(j,:),param.N);
        if ~isempty(tmp_pt)
            % connection of shapes, assume in order
            tmp_param = [1,tmp_pt(end,:) x(1) y(1)];
            [con_x,con_y] = U_samplePts(tmp_param,N);
            tmp_pt =[tmp_pt; [con_x,con_y]; [x,y]];
            %tmp_shape_type = [size(tmp_pt,1)+[1,numel(con_x)],1;size(tmp_pt,1)+numel(con_x)+[1,numel(x)],1+m_shapes{i}(j,1)];
            %m_shape_type =[m_shape_type; tmp_shape_type];
            
        else
            tmp_pt =[x,y];
            %m_shape_type =[m_shape_type; [1,numel(x),1+m_shapes{i}(j,1)]];
        end
    end
    m_pts{i} =tmp_pt;
end


%% Step 4. Compute Deformation
% 4.1 initial 2D warping field
warp_field = zeros([im_sz(1:2) 3]);
if im_sz(3)==1
    warp_field(:,:,1) = im;
else
    warp_field(:,:,1) = rgb2gray(im);
end

for i=1:numel(m_pts)
    % pts to strips
    [Stripe,vxn,vyn] = p_sampleStripe(m_pts{i}(:,2)', m_pts{i}(:,1)', im, param.delta, param.N);
    % matting on strips
    % avoid numerical 0
    Mask = Stripe+param.matt_offset;
    Mask2 = Mask;
    % mask positive and negative
    Mask2(param.matt_ran(1):param.matt_ran(2),:,:) = 1;
    Mask2(end-(param.matt_ran(2):-1:param.matt_ran(1))+1,:,:) = 0;
    Stripe_m =runMattingfun(Mask,Mask2);
    
    % strips to deviation
    param.shape_type=1;
    if sqrt((m_pts{i}(1,2)-m_pts{i}(end,2))^2+(m_pts{i}(1,1)-m_pts{i}(end,1))^2)<param.contour_close
        % closed shape
        param.shape_type=2;
    end
    [raw_edge, filter_edge] = P_analysis(Stripe_m, param);
    %filter_edge=raw_edge;
    % deviation to 2D field
    warp_field(:,:,2) = warp_field(:,:,2) + U_pointsToGrid(m_pts{i}(:,[2,1]), filter_edge(:).*vxn(:), im_sz(1:2),'nearest');
    warp_field(:,:,3) = warp_field(:,:,3) + U_pointsToGrid(m_pts{i}(:,[2,1]), filter_edge(:).*vyn(:), im_sz(1:2),'nearest');
end
warp_field(:,:,2:3)=warp_field(:,:,2:3)*param.alpha;
%%
% 4.2 diffused 2D warping field
[ux, uy] = P_synFlow(warp_field(:,:,2)~=0,warp_field, param.Dist_th);

%% Step 5. Amplify Deformation
% 2.3 amplification result
im_amp = P_warpImage(im,ux,uy);

%%
% 3. display
figure(2),
subplot(321),title('labeled geometry')
imshow(im),hold on
for sid=1:numel(m_shapes)
    U_plotGeo(m_shapes{sid}(:,1),m_shapes{sid}(:,2:end))
end


subplot(322),title('deviation signal: raw v.s. filtered')
plot([raw_edge; filter_edge]'),axis tight
legend('raw','filtered')

subplot(312),title('warping field: [x,y]')
imagesc([warp_field(:,:,2) warp_field(:,:,3)])
colormap jet; colorbar


subplot(313),title('image: input v.s. magnified')
imagesc([im zeros(im_sz(1),20,im_sz(3)) im_amp])



