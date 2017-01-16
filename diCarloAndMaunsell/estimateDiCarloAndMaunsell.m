function [lambda, betas] = estimateDiCarloAndMaunsell(binned, RT, tBinned,varargin)
% estimateDiCarloAndMaunsell(dN, leftLim, tBinned, varargin)
%
%   binned - binned spike trains at 1 ms
%   RT - RTs
%   tBinned - time for the binned spike trains
%
% CC, Shenoylab, 2017
nIterations = 2;
useSlow = false;
assignopts(who, varargin);
Latency = [-0.4:0.02:1];
Betas = [-0.4:0.1:1.5];

if size(binned,1)~=length(tBinned)
   fprintf('\n I found that the firing rate array has the wrong dimensions. I am transposing');
   binned = binned';
end

if max(RT) > 10000
    fprintf('\n The RTs seem suspiciously large, I am dividing by 1000. Ensure your RTs are in the normal range');
    RT = RT./1000;
    
end
    

for iter=1:nIterations
    fprintf('\n %d',iter);
    allErr = [];
    rates = [];
    for nLat = 1:length(Latency)
        if mod(nLat,10)==0
            fprintf('\n Testing lambda =%3.2f, ',Latency(nLat));
        end
        for nBeta = 1:length(Betas)
            preSpikes = [];
            postSpikes = [];
            Spre = [];
            Spost = [];
            if useSlow
                for trialId = 1:size(binned,2)
                    tCurr = Latency(nLat)*[nanmean(RT)./1000] + Betas(nBeta).*([RT(trialId) - nanmean(RT)]./1000);
                    %                 keyboard
                    %
                    preSpikes = [preSpikes binned(tBinned < tCurr,trialId)'];
                    postSpikes = [postSpikes binned(tBinned >= tCurr & tBinned < (RT(trialId)+100)./1000,trialId)'];
                end
            else
                [A,preSpikesFast,postSpikesFast] = fastCalc(RT, binned', Latency(nLat), Betas(nBeta), tBinned);
                
                 preSpikes = preSpikesFast(1:A(1));
                 postSpikes = postSpikesFast(1:A(2));
            end
            Spre = [preSpikes - mean(preSpikes)].^2;
            Spost = [postSpikes - mean(postSpikes)].^2;
            allErr(nLat, nBeta) = sum(Spre(:)) + sum(Spost);
            % rates(nLat,nBeta,:) = [nanmean(Spre) nanmean(Spost)];
        end
    end
    
    ErrSurf{iter} = allErr;
    
    O = allErr;
    [IX,IY] = find(O < prctile(O(:),1));
    
    [ifx,ify] = find(O == min(O(:)));
    
    try
        if ~isempty(ifx)
            lambda(iter) = Latency(ifx(1));
        else
            lambda(iter) = NaN;
        end
        if ~isempty(ify)
            betas(iter) = Betas(ify(1));
        else
            betas(iter) = NaN;
        end
    catch
        lambda(iter) = NaN;
        betas(iter) = NaN;
    end
    
    tmpLat = Latency(IX);
    tmpBeta = Betas(IY);
    
    
    Latency = [min(tmpLat):0.0025:max(tmpLat)];
    Betas = [min(tmpBeta):.0025:max(tmpBeta)];
    
end



