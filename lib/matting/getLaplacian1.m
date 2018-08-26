function [A,A1]=getLaplacian1(I,consts,epsilon,win_size)
  
  if (~exist('epsilon','var'))
    epsilon=0.0000001;
  end  
  if (isempty(epsilon))
    epsilon=0.0000001;
  end
  if (~exist('win_size','var'))
    win_size=1;
  end     
  if (isempty(win_size))
    win_size=1;
  end     

  neb_size=(win_size*2+1)^2;
  [h,w,c]=size(I);
  n=h; m=w;
  img_size=w*h;
  consts=imerode(consts,ones(win_size*2+1));
  
  indsM=reshape([1:img_size],h,w);

  tlen=sum(sum(1-consts(win_size+1:end-win_size,win_size+1:end-win_size)))*(neb_size^2);

  row_inds=zeros(tlen ,1);
  col_inds=zeros(tlen,1);
  vals=zeros(tlen,1);
  len=0;
  % The following is  faster vsion of the original in matlab vectorization
  % style
  kernel = ones(2*win_size+1)/neb_size;
  win_mu_all_orig = imfilter(I, kernel);
  win_mu_all = permute(win_mu_all_orig, [ 3 1 2]);
  CMatrix = zeros(3,3,size(I,1), size(I,2));
  for c1 = 1:3
      for c2 = 1:3
        CMatrix(c1,c2,:,:) = imfilter(I(:,:,c1).*I(:,:,c2),kernel) - win_mu_all_orig(:,:,c1).*win_mu_all_orig(:,:,c2)+epsilon/neb_size*(c1==c2);
        
      end
  end
  winI_pre = zeros(neb_size,3,size(I,1),size(I,2));
  win_inds_pre = zeros(neb_size,size(I,1),size(I,2));
  count = 1';
  for k = -1:1
      for j = -1:1
            winI_pre(count,:,:,:) = permute(circshift(I, -[j,k,0]), [3 1 2])-win_mu_all;            
            win_inds_pre(count,:,:) = circshift(indsM,-[j,k]);
            count = count +1;
      end
  end
  killIndices = repmat(permute(~consts(2:end-1,2:end-1),[3 4 1 2]), [9,9,1,1]);
  row_inds = repmat(permute(win_inds_pre(:,2:end-1,2:end-1),[1 4 2 3]),[1 neb_size 1 1]);
  row_inds = row_inds(killIndices);
  col_inds = repmat(permute(win_inds_pre(:,2:end-1,2:end-1),[4 1 2 3]),[neb_size 1 1 1]);
  col_inds = col_inds(killIndices);
  
  sz = size(CMatrix);
  CMatrix2 = reshape(CMatrix,[sz(1),sz(2),sz(3)*sz(4)]);
  
  RHS = permute(MultiSolver(CMatrix2, eye(3)),[1 3 2]);
  CMatrix2 = reshape(RHS, sz);
  
  tvals2 = zeros(neb_size,neb_size,size(I,1), size(I,2));
  for k = 1:neb_size
      for j = k:neb_size
            tvals2(k,j,:,:) = (1+sum(sum(repmat(winI_pre(k,:,:,:),[3 1 1 1]).*repmat(permute(winI_pre(j,:,:,:),[2 1 3 4]),[1 3 1 1]).*CMatrix2,1),2))/neb_size;
            tvals2(j,k,:,:) = tvals2(k,j,:,:);
      end
  end
  vals2 = tvals2(:,:,2:end-1,2:end-1);
  vals2 = vals2(:);
  vals = vals2(killIndices);
  %{
  for j=1+win_size:w-win_size
    for i=win_size+1:h-win_size
      if (consts(i,j))
        continue
      end  
      %win_inds=indsM(i-win_size:i+win_size,j-win_size:j+win_size);
      %win_inds=win_inds(:);
      %win_inds = win_inds_pre(:,i,j);
      %assert(mean(abs(win_inds(:)-win_inds_orig(:)))<1e-8, sprintf('Failed: %d, %d', i,j));
      %winI=I(i-win_size:i+win_size,j-win_size:j+win_size,:);
      %winI=reshape(winI,neb_size,c);
      %win_mu = (win_mu_all(:,i,j));
      %win_mu=mean(winI,1)';
      %assert(mean(abs((win_mu- (win_mu_all(:,i,j)))))< 1e-8,sprintf('Failed: %d, %d', i,j));
      %win_var_orig=inv(winI'*winI/neb_size-win_mu*win_mu' +epsilon/neb_size*eye(c));
      %win_var = inv(CMatrix(:,:,i,j));
      win_var=CMatrix2(:,:,i,j);
      %assert(mean(abs(win_var(:)-win_var_orig(:)))/mean(abs(win_var_orig(:)))<1e-8, sprintf('Failed: %d, %d', i,j));
      
      %winI=winI-repmat(win_mu',neb_size,1);
      winI = winI_pre(:,:,i,j);
      %assert(mean(mean(abs(winI-winI_pre(:,:,i,j))))<1e-8, sprintf('Failed: %d, %d', i,j));
      %winI = bsxfun(@minus, winI, win_mu');
      %tvals=(1+winI*win_var*winI')/neb_size;
      %assert(mean(mean(abs(tvals2(:,:,i,j)-tvals)))<1e-8,  sprintf('Failed: %d %d',i,j));
      %tvals=(1+winI*(win_var\winI'))/neb_size;
      % assert(mean(abs(tvals_orig(:)-tvals(:)))<1e-8, 'Failed: %d %d',i,j);
      tvals = tvals2(:,:,i,j);
      %row_inds(1+len:neb_size^2+len)=reshape(repmat(win_inds,1,neb_size), neb_size^2,1);
      %col_inds(1+len:neb_size^2+len)=reshape(repmat(win_inds',neb_size,1), neb_size^2,1);
      vals(1+len:neb_size^2+len)=tvals(:);
      len=len+neb_size^2;
    end
  end  
  %}
  %vals=vals(1:len);
  %row_inds=row_inds(1:len);
  %col_inds=col_inds(1:len);
  A=sparse(row_inds,col_inds,vals,img_size,img_size);
  
  sumA=sum(A,2);
  A=spdiags(sumA(:),0,img_size,img_size)-A;
  
return


