% this is an implemnatation that is based on Anat's colorization paper.
% Input: FlowMap is binary image indicating where are the "scribbles" i.e.,
% hard constraints on the flow.
% ImFlow is a 3 channel image. First channel: intensity image. 2 channel -
% flow in y, 3 channel - flow in x.

% donglai: 
% 5x speed up from vectorization
% ~2x speed up from cutting down matrix size (only connected pixels)
function [ux, uy]=U_SynFlow(FlowMap,ImFlow, Dist_th)
if(isempty(Dist_th))
    Dist_th = 50;
end
wd=1;
% pad image for boundary condition
FlowMap = padarray(FlowMap,[wd wd],'replicate');
ImFlow = padarray(ImFlow,[wd wd],'replicate');




n=size(ImFlow,1); m=size(ImFlow,2);
imgSize=n*m;
u = ImFlow(:,:,2);
v = ImFlow(:,:,3);
LocalDirection = zeros(n,m);
LocalDirection(FlowMap~=0) = atan2(v(FlowMap~=0), u(FlowMap~=0));
[~, IND] = bwdist(LocalDirection);
if(find(IND))
    LocalDirection = LocalDirection(IND);
end



D = double(bwdist(FlowMap));
% remove boundary
D([1:wd (end-wd+1):end],:)=Dist_th+1;
D(:,[1:wd (end-wd+1):end])=Dist_th+1;

% D = D-min(D(:));
% D = D./max(D(:));



num_w = (2*wd+1)^2;
patch_ind = reshape(bsxfun(@plus, (-wd:wd)'*n,(-wd:wd)),1,[]);
patch_cen = ceil(num_w/2);
patch_other = setdiff(1:num_w,patch_cen);

% [pts_x,pts_y]=ind2sub([n,m],pts);

[patch_dy,patch_dx] = meshgrid(-wd:wd,-wd:wd);
patch_dir = [patch_dy(:) patch_dx(:)];
% will be nan
patch_dir(patch_cen,:) =[];
patch_dir = bsxfun(@rdivide,patch_dir, sqrt(sum(patch_dir.^2,2)));

% area without the scribble
pts = reshape(find(D<Dist_th & ~FlowMap),[],1);
lblInds=find(FlowMap);


num_pt = numel(pts);

row_inds=repmat((1:num_pt)',[1,num_w]);
col_inds=bsxfun(@plus, pts,patch_ind);

% many unused
%pts_all = reshape(1:imgSize,[n,m]);
pts_all = unique(col_inds);
lblInds = intersect(lblInds,pts_all);
% pts_other = setdiff(pts_all((wd+1):end-wd,(wd+1):end-wd),pts);
pts_other = setdiff(pts_all,pts);
pts_other = [lblInds;setdiff(pts_other,lblInds)];
%{
% intersect(pts,lblInds))
% not empty as flow field has width
[xx,yy]=ind2sub([n,m],setdiff(lblInds,pts_all));
[xx2,yy2]=ind2sub([n,m],lblInds);
imshow(ImFlow(:,:,1));hold on;
plot(yy2,xx2,'r-')
plot(yy,xx,'bx')
imagesc(ImFlow(:,:,2))
%}
cosDir = cos(LocalDirection(pts));
sinDir = sin(LocalDirection(pts));
% N: number of pts index
% P: number of patch index
% Nx(P-1)
dir_val = squeeze(abs(sum(bsxfun(@times, patch_dir, reshape([sinDir';cosDir'],[1,2,num_pt])),2)))';



gvals = ImFlow(col_inds);
t_val = ImFlow(pts);
c_var = mean(bsxfun(@minus,gvals,mean(gvals,2)).^2,2);
csig=c_var*10;
mgv=min(bsxfun(@minus,gvals(:,patch_other),t_val).^2,[],2);

csig(csig<-mgv/log(0.01)) = -mgv(csig<-mgv/log(0.01))/log(0.01);
csig(csig<0.000002)=0.000002;
gvals=exp(-bsxfun(@rdivide,bsxfun(@minus,gvals(:,patch_other),t_val).^2,csig)).*dir_val;
gvals=bsxfun(@rdivide,gvals,sum(gvals,2));

vals(:,patch_other) = -bsxfun(@times,gvals,exp(-D(pts)/(Dist_th/2)^2));                
vals(:,patch_cen)=1;


num_pt_other = numel(pts_other);
num_pt_all = num_pt_other+num_pt;
% need to relabel pts
relabel = zeros(1,imgSize);
relabel([pts;pts_other])= (1:num_pt_all)';
A=sparse([row_inds(:); ((num_pt+1):num_pt_all)'],relabel([col_inds(:); pts_other]),[vals(:); ones(num_pt_other,1)],num_pt_all,num_pt_all);
b=zeros(size(A,1),1);

nI = ImFlow;
for t=2:3
    curIm=ImFlow(:,:,t);
    b(num_pt+(1:numel(lblInds)))=curIm(lblInds);
    new_vals=A\b;
    outIm = nI(:,:,t);
    outIm([pts;pts_other]) = new_vals;    
    nI(:,:,t)=outIm;
end


ux = nI((wd+1):end-wd,(wd+1):end-wd,2);
uy = nI((wd+1):end-wd,(wd+1):end-wd,3);

%{
im0=ones(20);
aa= zeros(20);
bb= zeros(20);
aa(5:15,9:11)=1;
[ux, uy] = U_SynFlow(aa~=0,cat(3,im0,aa,bb), 3);
[ux2, uy2] = U_SynFlow_old(aa~=0,cat(3,im0,aa,bb), 3);
imagesc([[ux uy];[ux2 uy2]])
%}
