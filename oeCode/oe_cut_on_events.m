function [D,uAB,dims] = oe_cut_on_events(data,events,win,Fs)
% [D,uAB,dims] = oe_cut_on_events(data,events,win,Fs)
%
% data:     samples x channel
% events:   [timestamp valuesA valuesB]
% win:      [onset offset]
% Fs:       sampling rate
% 
% D:        time x channel x rep x A x B
% uAB:      unique event values corresponding to {A, B}
% dims:     structure with fields named for each dimension
%


if nargin < 4 || isempty(Fs), Fs = 1; end


swin = round(win .* Fs);
svec = swin(1):swin(2)-1;

ons = round(events(:,1) .* Fs);

evA = events(:,2); uevA = unique(evA);
evB = events(:,3); uevB = unique(evB);

c = sum(uevA(1) == evA & uevB(1) == evB);
D = zeros(length(svec),size(data,2),c,length(uevA),length(uevB),'like',data);
for i = 1:length(uevA)
    for j = 1:length(uevB)
        ind = uevA(i) == evA & uevB(j) == evB;
                
        sidx = ons(ind)+svec;
        for k = 1:size(sidx,1)
            tmp(:,:,k) = data(sidx(k,:),:);
        end
        
        D(:,:,:,i,j) = tmp;
    end
end

uAB = {uevA, uevB};

dims.Time        = 1;
dims.Channels    = 2;
dims.Repetition  = 3;
dims.Avalues     = 4;
dims.Bvalues     = 5;