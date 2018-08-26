function [L,Lan] = U_elliptic(a,b,t1,t2)
%
% ELLIPSEARC computes the arc length of an ellipse of equation:
%
%     x(t) = a.cos(t)
%     y(t) = b.sin(t)
%
% between angles t1 and t2 (in radians). The ellipse arc length is
% computed numerically by dividing the arc in small straight segments.
%
% If angles are omitted, the perimeter of the complete ellipse is given.
% The second output is the approximated perimeter of the complete ellipse
% using the Ramanujan formula (see http://en.wikipedia.org/wiki/Ellipse): 
%
%           L = pi.( 3.(a+b) - ((3.a+b).(a+3.b))^(1/2) )
%
%       See http://en.wikipedia.org/wiki/Ellipse
%
% Examples:
%
%   >> ellipsearc(10,5)
%
%   ans =
%
%      48.4422
%
%   >> ellipsearc(10,5,0,pi/3)
%
%   ans =
%
%       7.0496
%
%   >> Req=6378.137;                     % earth radius at equator
%   >> Rpo=6356.752;                     % earth radius at pole
%   >> ellipsearc(Req,Rpo,-pi/2,pi/2)    % meridian length
%
%   ans =
%
%     2.0004e+004
%
%   Luc Masset, luc.masset@ulg.ac.be (March 2010) 

%arguments
if nargin ~= 2 & nargin ~= 4,
 error('Requires two or four inputs.')
 return
end
if nargin == 2,
 t1=0;
 t2=pi/2;
end

%approximated analytical perimeter (Ramanujan)
Lan=pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)));

%tests
if t2 <= t1,
 L=0;
 return
end

%number of straight segments (complete ellipse)
nseg=100000;    % you may adjust this parameter

%number of straight segments (approx. same precision for any ellipse arc)
nseg=nseg*(t2-t1)/2/pi;
nseg=max(nseg,2);

%numerical computation of arc length (division in small straight segments)
tt=(t1:(t2-t1)/nseg:t2);
xx=a*cos(tt);
yy=b*sin(tt);
dx=diff(xx);
dy=diff(yy);
L=sum(sqrt(dx.^2+dy.^2));

%complete ellipse
if nargin == 2,
 L=4*L;
end

return