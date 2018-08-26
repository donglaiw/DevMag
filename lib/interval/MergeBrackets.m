function [lower,upper] = MergeBrackets(Left, Right)

% Detect when right < left (empty Ii), and later remove it (line #29, 30) 
Right = Right(:); 
Left = Left(:);

IdxKeep = Right>=Left; 
Right = Right(IdxKeep); 
Left = Left(IdxKeep); 
% sort the rest by left bound 
[Left,iorder] = sort(Left); 
Right = Right(iorder);

% Allocate, as we don't know yet the size, we assume the largest case 
lower = zeros(size(Left)); 
upper = zeros(size(Right));

% Nothing to do 
if isempty(lower) 
return 
end

% Initialize 
l = Left(1); 
u = Right(1); 
k = 0; 
% Loop on brakets 
for i=1:length(Left) 
if Left(i) > u % new Jk detected 
% Stack the old one 
k = k+1; 
lower(k) = l; 
upper(k) = u; 
% Reset l and u 
l = Left(i); 
u = Right(i); 
else 
u = max(u, Right(i)); 
end 
end 
% Stack the last one 
k = k+1; 
lower(k) = l; 
upper(k) = u; 
IdxKeep = true(k,1); 
% Remove the tails 
lower = lower(IdxKeep); 
upper = upper(IdxKeep);