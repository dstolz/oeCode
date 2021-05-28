%%  Create neighborhood 

 load('G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG_Data\ECoG_channel_map_layout.mat')
 p = layout.pos(1:end-2,:);
 n = [];
for i = 1:length(p)
    n(i).label = layout.label{i};
    
    x = p(i,1);
    y = p(i,2);
    
    z = sqrt((x-p(:,1)).^2+(y-p(:,2)).^2);
    
    ind = z > 0 & z < 1.5; % <= sqrt(2)
    
    n(i).neighblabel = layout.label(ind);
end
