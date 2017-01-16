%
nTrials = 100;
t = [-0.4:0.01:1.5]';
zeroTime = 0.0;
RT = 0.305+gamrnd(0.7,0.4,1,nTrials);
RT(RT > 1.5) = 1;
figure(1);
hist(RT*1000,30);
hold on;
set(gca,'visible','off')
getAxesP([200 1500],[200:200:1500],'RT (ms)',-1,-1,[0 30],[0 30],'Trials',200,1,[1 1]);

figure;
NeuronalLatency = 0.25;
simulatedBeta = 1;
simulatedBeta = 0.1;
simulatedLambda = NeuronalLatency./mean(RT);

[RT,i] = sort(RT);

figure(2);
spksRT = [];
spksRamp = [];
FR = [];
for k = 1:nTrials
    FR = 5*ones(1,length(t));
    % NeuronalLatency = NeuronalLatency + 1*[RT(k)-nanmean(RT)];
    FR(t > zeroTime + NeuronalLatency + simulatedBeta*(RT(k)-mean(RT)) & t < zeroTime + RT(k)+0.05) = 40;
    % FR(t > zeroTime + RT(k)-0.1 & t < zeroTime + RT(k)) = 50;
    spksRT(k).times = generateInhomPoisson(t(1:length(t)), FR );
    X = spksRT(k).times';
    line([X; X],[k+0.2*ones(1,length(spksRT(k).times)); k+0.8*ones(1,length(spksRT(k).times))],'color','k');
    
    hold on;
    line([RT(k); RT(k)],[k+0.2; k+0.8],'color','m');
end
axis tight;
getAxesP([t(1) t(end)+zeroTime]-zeroTime,[t(1):0.2:t(end)+zeroTime]-zeroTime,8,-2,'Time (ms)',[t(1) t(end)],[t(1):0.2:t(end)],-1,-10,'Time (ms)',[1 0]);
% drawLines(0);
set(gca,'visible','off');

drawLines(NeuronalLatency);
drawLines(0);
data = [];

%%




%%
lims = [min(t) max(t)]*1000;
[dN,tBinned] = binspikes(spksRT,1000,[lims]./1000);
leftLim = floor(RT*1000);
    

% Now choose an offset for each value of beta and neuronal latency.
Latency = [-0.4:0.02:1];
Betas = [-0.4:0.1:1.5];
tic;
for iter=1:2
    fprintf('\n %d',iter);
    allErr = [];
    rates = [];
    for nLat = 1:length(Latency)
        if mod(nLat,10)==0
            fprintf('..%3.2f, ',Latency(nLat));
        end
        for nBeta = 1:length(Betas)
            preSpikes = [];
            postSpikes = [];
            Spre = [];
            Spost = [];
%             for trialId = 1:size(dN,2)
%                 tCurr = Latency(nLat)*[nanmean(leftLim)./1000] + Betas(nBeta).*([leftLim(trialId) - nanmean(leftLim)]./1000);
%                 %                 keyboard
%                 %
%                 preSpikes = [preSpikes dN(tBinned < tCurr,trialId)'];
%                 postSpikes = [postSpikes dN(tBinned >= tCurr & tBinned < (leftLim(trialId)+100)./1000,trialId)'];
%             end
      
            
            [A,preSpikesFast,postSpikesFast] = fastCalc(leftLim, dN', Latency(nLat), Betas(nBeta), tBinned);
            
            preSpikesFast = preSpikesFast(1:A(1));
            postSpikesFast = postSpikesFast(1:A(2));
            
            preSpikes = preSpikesFast;
            postSpikes = postSpikesFast;
            
            Spre = [preSpikes - nanmean(preSpikes)].^2;
            Spost = [postSpikes - nanmean(postSpikes)].^2;
            allErr(nLat, nBeta) = nansum(Spre(:)) + nansum(Spost);
            rates(nLat,nBeta,:) = [nanmean(Spre) nanmean(Spost)];
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
    Betas = [min(tmpBeta):.005:max(tmpBeta)];
    
end
toc;

%%
% figure;

for trialId = 1:nTrials
    subplot(231);
    tCurr = lambda(2)*[nanmean(leftLim)./1000] + betas(2).*([leftLim(trialId) - nanmean(leftLim)]./1000);
    % line([tCurr; tCurr],[trialId+0.2; trialId+0.8],'color','b');
    data(trialId) = tCurr;
    
    hold on;
    line([RT(trialId); RT(trialId)],[trialId+0.2; trialId+0.8],'color','m');
end
drawLines(lambda(iter)*mean(leftLim)./1000);
plot(data,1:100,'r-');

text(-0.2, 101, sprintf('lambda=%3.2f, beta=%3.2f',lambda(2),betas(2)));
text(-0.2, 101, sprintf('simlambda=%3.2f, simbeta=%3.2f',simulatedLambda,simulatedBeta));


%%
iter = 2;
lambda(2) = simulatedLambda
leftLim = RT
betas = [0 0];
% figure;
for trialId = 1:nTrials
    subplot(231);
    tCurr = lambda(2)*[nanmean(leftLim)./1000] + betas(2).*([leftLim(trialId) - nanmean(leftLim)]./1000);
    % line([tCurr; tCurr],[trialId+0.2; trialId+0.8],'color','b');
    data(trialId) = tCurr;
    
    hold on;
    line([RT(trialId); RT(trialId)],[trialId+0.2; trialId+0.8],'color','m');
end
drawLines(lambda(iter)*mean(leftLim));
plot(data,1:nTrials,'r-');

text(-0.2, 101, sprintf('lambda=%3.2f, beta=%3.2f',lambda(2),betas(2)));
text(-0.2, 101, sprintf('simlambda=%3.2f, simbeta=%3.2f',simulatedLambda,simulatedBeta));