function Ystate = ash_panel(pimt,HH,T,rsd);
% cooper version for just aggregate shocks.
% =========================================================
% NOW SET UP THE RANDOM VARIABLES USED FOR THE SIMULATION
% =========================================================
ne = length(pimt); % shocks
rand('state',rsd);

rnd2 = rand(HH,T); % so draw a panel of shocks

np=HH;
% create cumulative pimt
cpimt = pimt*triu(ones(ne));

cp70 = pimt^70*triu(ones(ne));
[id1x,id1y] = max(((ones(np,1)*cp70(1,:) - rnd2(1:np,1)*ones(1,ne)).^-1)');
Statei = id1y';
%Statei = ones(np,1)*round((ne+1)/2);

for t=2:T
 id1 = cpimt(Statei(:,t-1)',:);
 [id1x,id1y] = max(((id1 - rnd2(1:np,t)*ones(1,ne)).^-1)');
 Statei(:,t) = id1y';
end

Ystate = Statei;

