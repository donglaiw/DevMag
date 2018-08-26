function [rawDeviationSignal, filteredDeviationSignal, ver_c] = U_analysis(inputStrip, param, m_pts)
% for debug
DISPLAYFIGURE =0;

shape_type = param.shape_type;
% # of samples per pixel
if(~isfield(param,'N'))
    N = 1;
else
    N = param.N;
end
% low cut off  (in pixels)
if(~isfield(param, 'lambda_l') )
    lambda_l = 10*N;
else
    lambda_l = param.lambda_l;
end
% high cut off (in pixels)
if(~isfield(param,'lambda_h'))
    lambda_h = 300*N;
else
    lambda_h = param.lambda_h;
end

if(~isfield(param,'mag_curve'))
    mag_curve = 2;
else
    mag_curve = param.mag_curve;
end
% W = fspecial('gaussian', [size(inputStrip,1) 1], size(inputStrip,1)/10);
% W = repmat(W, [1 size(inputStrip,2)]);
W = ones(size(inputStrip));
%inputStrip = inputStrip.*win;
epsilon = 1e-3;

[Gx, Gy] = imgradientxy(inputStrip,'centralDifference');
% inputStrip = Gy;
% [Gx, Gy] = imgradientxy(inputStrip,'centralDifference');
Rm = median(inputStrip,2);
Rmx = median(Gy,2);
rawDeviationSignal = zeros(1,size(inputStrip,2));
for i=1:size(inputStrip,2)
    if(DISPLAYFIGURE)
        if(i==2)
            figure;
        end
        if(mod(i,10)==0)
            plot(inputStrip(:,i)); hold on
        end
    end
  
    rawDeviationSignal(i) = sum((inputStrip(:,i)-Rm).*W(:,i).*Rmx)./sum(Rmx.*Rmx.*W(:,i));
end


if (param.detrend);
   rawDeviationSignal =detrend(rawDeviationSignal); 
end
if (param.removeMean)
    rawDeviationSignal =rawDeviationSignal-mean(rawDeviationSignal); 
end
if (shape_type==1)
    x = 1:length(rawDeviationSignal);

    curve_fun = [];
    ver_c = zeros(size(x));
    NoCSignal = rawDeviationSignal-ver_c;
    if(DISPLAYFIGURE)
        f1 = figure;
        subplot(1,2,1); plot(ver_c); hold on; plot(rawDeviationSignal,'r');
    end
    
    %bandpassing
    frequency_highCutoff = 2/lambda_l; % 2 pixels per cycle maps to 1 (the normalized Nyquist frequency)
    frequency_lowCutoff = 2/lambda_h;
    signalLength = size(rawDeviationSignal,2);
    B = fir1(signalLength, [frequency_lowCutoff, frequency_highCutoff]);
    % bandpassing
    signalIn = [2*NoCSignal(:,1)-NoCSignal(:,end-1:-1:2) NoCSignal 2*NoCSignal(:,end)-NoCSignal(:,end-1:-1:2);];
    %signalIn = [repmat(NoCSignal(:,1), [1 numel(NoCSignal)]) NoCSignal repmat(NoCSignal(:,end), [1 numel(NoCSignal)]);];
    signalFiltered = convn(signalIn,B,'same');
    if (param.antialiasing)
        slope = m_pts(1,:) - m_pts(end,:);
        orientation =  atan2(slope(2), slope(1));
        signalFiltered = removeAliasingFrequencies(signalFiltered, param.N, orientation);
    end
    signalFiltered = signalFiltered(:,signalLength+1:2*signalLength);
    %signalFiltered = NoCSignal;
    if(DISPLAYFIGURE)
        figure(f1);
        subplot(1,2,2);
        plot(NoCSignal,'r');
        plot(signalFiltered);
        hold on;
    end
elseif (shape_type==2)
    pts2Sample = numel(rawDeviationSignal);
    if (param.removeGlobalMotionCircleEllipse)
       lambda_h = min(lambda_h, 0.95*pts2Sample); 
       
    end
    frequency_highCutoff = 2/lambda_l; % 2 pixels per cycle maps to 1 (the normalized Nyquist frequency)
    frequency_lowCutoff = 2/lambda_h;
    if (frequency_lowCutoff < frequency_highCutoff)
        B = fir1(pts2Sample, [frequency_lowCutoff, frequency_highCutoff],kaiser(pts2Sample+1,0.5));
        %B = fir1(pts2Sample, [frequency_lowCutoff, frequency_highCutoff]);
        % bandpassing
        signalIn = [rawDeviationSignal(:); rawDeviationSignal(:); rawDeviationSignal(:);];
        signalFiltered = conv(signalIn,B,'same');
        if (param.antialiasing)
            slope = m_pts(1,:) - m_pts(end,:);
            orientation =  atan2(slope(2), slope(1));
            signalFiltered = removeAliasingFrequencies(signalFiltered, param.N, orientation);
        end
        signalFiltered = signalFiltered(pts2Sample+1:2*pts2Sample);
    else
       signalFiltered = zeros(size(rawDeviationSignal)); 
    end
    
    filteredDeviationSignal = signalFiltered';
    ver_c = 0;
    curve_fun = 0;
    if(DISPLAYFIGURE)
        figure(3);
        plot(rawDeviationSignal); hold on;
        plot(signalFiltered,'r');
    end
end

if (shape_type==1)
    if(mag_curve==1)
        filteredDeviationSignal = ver_c;
    else
        %warpingSignal =signalFiltered;
        filteredDeviationSignal = signalFiltered;
    end
    
    if(mag_curve==2)
        filteredDeviationSignal = filteredDeviationSignal + ver_c;
    end
end


if isfield(param,'signal_0') && param.signal_0 ~= 0
    switch param.signal_0
        case 1
            rawDeviationSignal = rawDeviationSignal-rawDeviationSignal(1);
            filteredDeviationSignal=filteredDeviationSignal-filteredDeviationSignal(1);
        case 2
            rawDeviationSignal = rawDeviationSignal-rawDeviationSignal(end);
            filteredDeviationSignal=filteredDeviationSignal-filteredDeviationSignal(end);
    end
end

if (param.detrend);
   filteredDeviationSignal = filteredDeviationSignal -linspace(filteredDeviationSignal(1), filteredDeviationSignal(end), numel(filteredDeviationSignal));
end

