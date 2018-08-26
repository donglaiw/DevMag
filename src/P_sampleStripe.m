% the function get a list of 2D points and sample in the local normal
% direction the image Im

function [Stripe,vxn,vyn] = U_SampleStripe(xPts, yPts, Im, delta, samplesPerPixel)

Np = length(xPts);
[w, h,c] = size(Im);

%comput the local gradient
[vxn, vyn] = computeNormalDirection(xPts, yPts);
k=1;

out = zeros([1 numel(vxn) c]);
ind = 1/samplesPerPixel:1/samplesPerPixel:delta;
Stripe = zeros([1+2*numel(ind) numel(vxn) c]);


xPtsj = [];
yPtsj = [];
totalRows = 0;
for j=[-ind(end:-1:1) 0 ind]
    xPtsj = [xPtsj xPts + vxn*j];
    yPtsj = [yPtsj yPts + vyn*j];
    totalRows = totalRows + 1;
end

for i=1:c
    Stripe(:,:,i) = reshape(U_samplePoint(Im(:,:,i), [xPtsj;yPtsj]),[numel(xPts), totalRows])';
end
%{
    Stripe(k,:,:) = out;
    k = k+1;    
end
%}
%}
%{

for j=[-ind(end:-1:1) 0 ind]
    xPtsj = xPts + vxn*j;
    yPtsj = yPts + vyn*j;
    for i=1:c
        out(:,:,i) = U_samplePoint(Im(:,:,i), [xPtsj;yPtsj]);
    end
    Stripe(k,:,:) = out;
    k = k+1;    
end
%}
% need to handle out of boundary
Stripe=real(Stripe);
Stripe(isnan(Stripe))=0;
