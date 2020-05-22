clc
clear;
close all;

%% Parameterization
aalpha = 0.3; % capital share
ttheta = 0.8; % decreasing return technology
nnu = 1 - (1-aalpha).*ttheta; % technology exponent

ddelta = 0.15; % depreciation
ppsi = 0.15; % adjustment cost parameter
kksi = 0.5; % collateral constraint
llambda = 0.05; % external financing cost
r = 0.04; % risk-free interest rate

% Technology Process
rhoLogZ = 0.7; % correlation of technology AR1 
sigmaLogZ = 0.2; % standard deviation of technology AR1
nZ = 3;
muZ = 1;
[vLogZ,mTranZ] = tauchen(nZ,log(muZ),rhoLogZ,sigmaLogZ,3);
vZ = exp(vLogZ)';
zMax = max(vZ);
zMin = min(vZ);

% Tax Process
taoHigh = 0.3;
taoLow = 0.1;
nTao = 2;
vTao = [taoHigh; taoLow];

%% To calibrate
rhoTao = 0.8; % s.t. it changes about once every 5 years
mTranTao = [rhoTao,     1 - rhoTao;...
                    1 - rhoTao,     rhoTao];

cchi = 2; % s.t. workers work about 1/3 of the time, i.e., L=1/3; 
              % some numbers between 1.5 and 2.5 should work;

%% Combine the two shocks into one shock
mTran_z_tao = kron(mTranZ,mTranTao);
[A1,A2] = meshgrid(vTao,vZ);
temp=cat(2,A2',A1');
mGrid_z_tao=reshape(temp,[],2);
N = nTao * nZ;

%% Functions
% syms z nu K alpha L theta w
% f_product_sym = z.^nu.*(K.^alpha .* L.^(1-alpha)).^theta;
% f_wage_sym = simplify(diff(f_product_sym,L));
% f_wage_sym_subs=subs(f_wage_sym,[alpha,theta,nu],[aalpha,ttheta,nnu]);
% simplify(f_wage_sym_subs)
% solve(f_wage_sym_subs==w,L)

f_product =@(aalpha,k,labor,nnu,ttheta,z)z.^nnu.*(k.^aalpha .* labor.^(1-aalpha)).^ttheta;
f_labor = @(aalpha,k,nnu,ttheta,wage,z)(wage./(z.^(nnu).*ttheta.*(1 - aalpha).*k.^(aalpha .* ttheta))).^(1./(ttheta - 1 - aalpha.*ttheta));

f_adjust = @(ddelta,k,kPrime,ppsi)ppsi.*(kPrime - (1-ddelta)*k).^2./(2*k);

f_bond_limit_upper = @(ddelta,kksi,kPrime)kksi.*(1 - ddelta).*kPrime;
bond_limit_lower = 0;

%% Grid for State Variables
nK = 20; % number of capital points
nP = 15; % number of bond points

% kSteadyState = (aalpha .* muZ ./(r + ddelta)).^(1 ./ (1-aalpha))
kSteadyState = (aalpha .*ttheta .* muZ ./(r + ddelta)).^(1 ./ (1 - aalpha .*ttheta));
kSpread = 0.5;

% kMax = (zMax ./ ddelta) .^(1./(1 - aalpha.*ttheta ));
% kMin = 0.0001;
% kMax = min(2*kSteadyState, kMax); % Tighten the grid

kMax = kSteadyState*(1 + kSpread);
kMin = kSteadyState *(1 - kSpread);

vK = curvspace(kMin,kMax,nK,2)'; % I use curved grid to enhance accuracy

vP_max = f_bond_limit_upper(ddelta,kksi,vK);
pMin = 0;

mP = zeros(nP,nK); % each column is bond grid for each k state
for iK = 1:nK
    mP(:,iK) = curvspace(pMin,vP_max(iK),nP,1/2)';
end

temp = zeros(nZ*nTao*nP,4,nK);
mGrid = zeros(nZ*nTao*nP*nK,4);

for iK = 1:nK
    temp(:,1:(end-1),iK) = gridmake(vZ,vTao,mP(:,iK));
    temp(:,end,iK)=vK(iK);
    mGrid(((iK-1)*nZ*nTao*nP+1):(iK)*nZ*nTao*nP,:) = temp(:,:,iK);
end


