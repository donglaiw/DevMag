% Takes in a list of coordinates and a value for each coordinate and
% returns a grid

% Coordinates Nx2 matrix, each row is x, y
function [ output ] = U_pointsToGrid( coordinates, values, imageSize, interpolationMethod)
    if (nargin < 4)
       interpolationMethod = 'nearest'; 
    end

    assert(numel(imageSize) == 2, 'Expects imageSize to have only a width and a height.')
    assert(size(coordinates,2) == 2, 'Coordinates must be an Nx2 matrix.');
    assert(size(coordinates,1) == size(values,1), 'Coordinates and values must have the same number of rows.');
    
    % Forward warping splatting
    output = zeros(imageSize);
    weights = zeros(imageSize);
    kernel = ones(3);
    for k = 1:size(coordinates,1)
       x = coordinates(k,1);
       y = coordinates(k,2);
       if (strcmp(interpolationMethod, 'nearest'))
           xp = round(x);
           yp = round(y);
           xBounds = and(xp >= 2, xp <= imageSize(2)-1);
           yBounds = and(yp >= 2, yp <= imageSize(1)-1);
           if (and(xBounds,yBounds))
               output(yp+(-1:1), xp+(-1:1)) = values(k)*kernel + output(yp+(-1:1), xp+(-1:1)) ;
               weights(yp+(-1:1), xp+(-1:1)) = kernel + weights(yp+(-1:1), xp+(-1:1));
           end
       elseif(strcmp(interpolationMethod, 'linear'))
           error('Linear interpolation doesnt work.');
           x1 = floor(x);
           x2 = x1 + 1;
           y1 = floor(y);
           y2 = y1 + 1;
           weightx = x-x1;
           weighty = y-y1;
           
           output(y1, x1) = output(y1,x1) + values(k).*(1-weightx).*(1-weighty);
           output(y1, x2) = output(y1,x2) + values(k).*(weightx).*(1-weighty);
           output(y2, x1) = output(y2,x1) + values(k).*(1-weightx).*(weighty);
           output(y2, x2) = output(y2,x2) + values(k).*(weightx).*(weighty);
           
       end
       
    end
    output = output./(eps+weights);

end

