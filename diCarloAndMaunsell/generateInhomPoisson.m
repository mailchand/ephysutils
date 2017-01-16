function spikes = generateInhomPoisson(time, rate)
% spikes = generateInhomPoisson(time, rate)
%
% generates spike times from a time-dependent rate function
% as an inhomogeneous Poisson process
% rate function is linearly interpolated between the samples
%
% inputs:
%   time    -   array (N,1), time stamps for the rate samples
%   rate    -   array (N,1), rate samples
% output:
%   spikes  -   spike times
%
%
%
% (c) Tatiana Engel, Boahen and Shenoy labs, 2016
% check dimensions
if (size(time,2)~=1)
    time = time';
    if (size(time,2)~=1)
        error('Time should be an Nx1 array.');
    end;
end;
if (size(rate,2)~=1)
    rate = rate';
    if (size(rate,2)~=1)
        error('Rate should be an Nx1 array.');
    end;
end;
if (all(size(time)~=size(rate)))
    error('Dimensions of time and rate vectors should be the same');
end;


% calculate cumulative rate
deltaT = time(2:end)-time(1:end-1);
r = [0; cumsum( rate(1:end-1).*deltaT )];
deltaR = r(2:end)-r(1:end-1);

% generate 1.5 as many spikes as expected on average for exponential distribution with rate 1
numX = ceil( 1.5*r(end) );

% generate exponential distributed spikes with the average rate 1
notEnough = true;
x = [];
xend = 0;
while (notEnough)
    x = [ x; xend+cumsum( exprnd(1.0, numX, 1) )];
    % check that we generated enough spikes
    xend = x(end);
    notEnough = xend<r(end);
end
% trim extra spikes
x(x>r(end)) = [];

% for each x find index of the last rate which is smaller than x
indJ = arrayfun(@(y) find(r<y,1,'last'), x);

% compute resclaed spike times
try
    spikes = [];
    if ~isempty(indJ)
        spikes = time(indJ) + (x-r(indJ)).*deltaT(indJ)./deltaR(indJ);
    end
catch
    
end

end