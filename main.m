clc
clear;
close all;

% Difference from the toymodel.pdf
% \rhoLogZ is 0.9, instead of 0.7 in the document. Otherwise the collateral constraint is always binding

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
rhoLogZ = 0.9; % correlation of technology AR1 
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
rhoTao = 0.7; % s.t. it changes about once every 5 years
mTranTao = [rhoTao,     1 - rhoTao;...
                    1 - rhoTao,     rhoTao];
% Ystate = ash_panel(mTranTao,1,100,nTao) % for simulations

cchi = 2; % s.t. workers work about 1/3 of the time, i.e., L=1/3; 
              % some numbers between 1.5 and 2.5 should work;



%% Combine the two shocks into one shock
mTran_z_tao = kron(mTranZ, mTranTao);
[A1,A2] = meshgrid(vTao,vZ);
temp = cat(2,A2',A1');
mGrid_z_tao=reshape(temp,[],2);
nShock = nTao * nZ;

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
f_invest = @(ddelta,k,kPrime)kPrime - (1-ddelta).*k;

f_divident = @(profitAfterTax,invest,costAdjust,p,pPrime,r,ttao)...
                            (profitAfterTax - invest - costAdjust + pPrime - p.*(1 + r.*(1-ttao)) )...
                             .*(1 + llambda.*(profitAfterTax - invest - costAdjust + pPrime - p.*(1 + r.*(1-ttao))<0))   ;                                                                    

f_bond_limit_upper = @(ddelta,kksi,kPrime)kksi.*(1 - ddelta).*kPrime;

%% Grid for State Variables
nK = 20; % number of capital points
nP = 10; % number of bond points

% kSteadyState = (aalpha .* muZ ./(r + ddelta)).^(1 ./ (1-aalpha))
kSteadyState = (aalpha .*ttheta .* muZ ./(r + ddelta)).^(1 ./ (1 - aalpha .*ttheta));
kSpread = 0.9;

% kMax = (zMax ./ ddelta) .^(1./(1 - aalpha.*ttheta ));
% kMin = 0.0001;
% kMax = min(2*kSteadyState, kMax); % Tighten the grid

kMax = kSteadyState*(1 + kSpread);
kMin = kSteadyState *(1 - kSpread);

% vK = curvspace(kMin,kMax,nK,2)'; % I use curved grid to enhance accuracy
vK = linspace(kMin,kMax,nK)'; % I use curved grid to enhance accuracy

vP_max = f_bond_limit_upper(ddelta,kksi,vK);
pMin = 0;

mP = zeros(nP,nK); % each column is the bond grid for each k state
for iK = 1:nK
%     mP(:,iK) = curvspace(pMin,vP_max(iK),nP,1/2)';
    mP(:,iK) = linspace(pMin,vP_max(iK),nP)';
end

% temp = zeros(nZ*nTao*nP,4,nK);
% mGrid = zeros(nZ*nTao*nP*nK,4);
% 
% for iK = 1:nK
%     temp(:,1:(end-1),iK) = gridmake(vTao,vZ,mP(:,iK));
%     temp(:,end,iK)=vK(iK);
%     mGrid(((iK-1)*nZ*nTao*nP+1):(iK)*nZ*nTao*nP,:) = temp(:,:,iK);
% end

% mTran_z_tao = kron(mTranZ,mTranTao);



%% To guess Wage to get equilibrium


% Strategy
% According to the utility function of the consumer, 
% marginal disutility of labor is cchi,
% marginal utility of consumption is wage/(p*(1+r) + wage*L - p' + E(.)),
% where p islast period's saving for this period, p' is this period's
% saving for next period, and E(.) is divident.
% I solve the firm's problem and computed the p, p', labor, E(.) on the
% aggregate level, and plug it into the MU of consumption wage/(p*(1+r) +
% wage*L - p' + E(.))
% Then I play around with wage to get the marginal (dis)utilities equal.

% I plotted the error (marginal utility minus marginal labor) against wage and
% see the former goes up with the latter.

% So I use bisection method to find wage.

% wage = 1.6;
% wage_min = 1.59;
% wage_max = 1.61;

wage_min = 0.6;
wage_max = 1;

wage_left = wage_min;
wage_right = wage_max;

tol = 1e-3;
error = 10;
maxIteration = 100;
ITER = 0;

tic
while abs(error) > tol && ITER < maxIteration
    
    ITER = ITER + 1;
    wage = 0.5 *(wage_left + wage_right);
    display(['wage = ',num2str(wage), '.    ITER = ',num2str(ITER)]);
    
% 
% nn= 10;
% vWage = linspace(1.59,1.61,nn)';
% vError = zeros(nn,1);

% for iWage = 1:nn
%     wage = vWage(iWage);

    %% Precompute some values to retrieve
    mLabor = f_labor(aalpha,vK,nnu,ttheta,wage,mGrid_z_tao(:,1)'); % nK by nShock
    mY = f_product(aalpha,vK,mLabor,nnu,ttheta,mGrid_z_tao(:,1)'); % nK by nShock
    mProfitAfterTax = (1-mGrid_z_tao(:,2)').*(mY - wage.*mLabor); % nK by nShock

    mCostAdjust = f_adjust(ddelta,vK,vK',ppsi); % nK by nK_prime
    mInvest = f_invest(ddelta,vK,vK'); % nK by nK_prime

    %% Value Function Iteration
    iter = 0;
    maxIter = 1000;
    distance = 10;
    tolerance = 1e-5;
    
    if ITER == 1
        mValue0 = repmat(mProfitAfterTax * (1 + 1/r),1,1,nP);
        mValue0 = permute(mValue0,[1,3,2]);
        % mValue0 = zeros(nK,nP,nShock);
        % for iP = 1:nP
        %     mValue0(:,iP,:)= mProfitAfterTax * (1 + 1/r);
        % end

        mValue1 = mValue0;

        mPolicyK = zeros(nK,nP,nShock);
        mPolicyP = zeros(nK,nP,nShock);
        mPolicyIndexK = zeros(nK,nP,nShock);
        mPolicyIndexP = zeros(nK,nP,nShock);
    end
    
%     tic
    while distance > tolerance && iter < maxIter
        iter = iter + 1;

        for iShock = 1:nShock % shock        
            ttao = mGrid_z_tao(iShock,2);

            for iK = 1:nK % k
                k = vK(iK);
                vP = mP(:,iK);
                profitAfterTax = mProfitAfterTax(iK,iShock);
                vInvest = mInvest(iK,:)'; % kPrime by 1
                vCostAdjust = mCostAdjust(iK,:)';% kPrime by 1

                for iP = 1:nP % p
                    p = vP(iP);

                    % kPrime by pPrime
                    mDivident = f_divident(profitAfterTax,vInvest,vCostAdjust,p,mP',r,ttao); % kPrime by pPrime

                    mValuePrime = zeros(nK,nP,nShock); % kPrime by pPrime by shockPrime
                    for iShockPrime = 1:nShock
                        mValuePrime(:,:,iShockPrime) = mTran_z_tao(iShock,iShockPrime) * mValue0(:,:,iShockPrime);
                    end
                    mValuePrimeExpected = sum(mValuePrime,3); % sum by the 3rd dimension to get a kPrime by pPrime matrix

                    x = mDivident + 1/(1+r)*mValuePrimeExpected;

                    [rows,cols]=find(x==max(x,[],'all'));
                    iKPrime = min(rows);
                    iPPrime = min(cols);
                    mPolicyIndexK(iK,iP,iShock) = iKPrime;
                    mPolicyIndexP(iK,iP,iShock) = iPPrime;

                    mPolicyK(iK,iP,iShock) = vK(iKPrime);
                    mPolicyP(iK,iP,iShock) = mP(iPPrime,iKPrime);

                    mValue1(iK,iP,iShock) = max(x,[],'all');

                end % p
            end % k
        end % shock

        distance = max(abs(mValue1 - mValue0),[],'all');
        mValue0 = mValue1;

        if mod(iter,100)==0
            display("Iteration = " + iter + ".  Difference = " + distance)
        end
    end% while
%     toc

    if distance <= tolerance
        display("Iteration = " + iter + ".  Difference = " + distance + ".   Converged.")
    else
        display("Iteration = " + iter + ".  Difference = " + distance + ".   Convergence failed.")
    end

    %% Compute Stationary Distribution
    mDist0 = (1/(nK*nP*nShock))*ones(nK,nP,nShock);
    distance=100;
    tolerance=1e-5;
    iteration=0;
    maxIter = 1000;

    while distance>tolerance && iteration < maxIter
        mDist1 = zeros(nK,nP,nShock);
        for iShock = 1:nShock
    %     for iTao=1:nTao
    %         for iZ=1:nZ
    %             iShock = (iZ - 1)*nTao + iTao;

            for iK=1:nK
    %                 vP = mP(:,iK);

                for iP=1:nP

                    iKPrime = mPolicyIndexK(iK,iP,iShock);          
                    iPPrime = mPolicyIndexP(iK,iP,iShock); 

                    prob = mDist0(iK,iP,iShock);

                    for iShockPrime=1:nShock
%                         iZPrime = ceil(iShockPrime/nTao);
%                         iTaoPrime = (mod(iShockPrime,nTao)==0)*nTao + (mod(iShockPrime,nTao)~=0)*mod(iShockPrime,nTao);
                        prob_shockPrime = prob*mTran_z_tao(iShock,iShockPrime);

                        mDist1(iKPrime,iPPrime,iShockPrime) = mDist1(iKPrime,iPPrime,iShockPrime) + prob_shockPrime;
                    end

                end % p 
            end % k
        end % shock

        distance=sum(abs(mDist0-mDist1),'all');
        mDist0 = mDist1;
        iteration = iteration + 1;
    end

    mDist1_4D= reshape(mDist1,nK,nP,nTao,nZ);

%% Aggregate Labor Demand
    % mDist1_3D = reshape(mDist1,nK,nP,nShock);

    mDist1_kz = sum(mDist1,[2]);
    mDist1_kz= reshape(mDist1_kz,nK,nShock);
    laborAgg = sum(mLabor.*mDist1_kz,'all');

    % mLaborDemand = f_labor(aalpha,vK,nnu,ttheta,wage,vZ'); % nK by nZ
    % mDist1_kz = sum(mDist1,[2 3]);
    % mDist1_kz= reshape(mDist1_kz,nK,nZ);
    % laborAgg = sum(mLaborDemand.*mDist1_kz,'all');

%% Aggregate Bond State
    mDist1_kp = sum(mDist1,[3]);
    pAgg = sum(mP'.*mDist1_kp,'all');

    % mDist1_kp = sum(mDist1,[3 4]);
    % pAgg = sum(mP'.*mDist1_kp,'all');

%% Aggregate Bond Demand
    pPrimeAgg = sum(mPolicyP.*mDist1,'all');

    assert(abs(pPrimeAgg - pAgg)<1e-4); % bond market clear

%% Aggregate Divident
    mProfitAfterTax_3D = permute(repmat(mProfitAfterTax,1,1,nP),[1 3 2]);

    mK_3D = repmat(vK,1,nP,nShock);
    mP_3D = repmat(mP',1,1,nShock);
    % mLabor = f_labor(aalpha,vK,nnu,ttheta,wage,mGrid_z_tao(:,1)'); % nK by nShock
    % mY = f_product(aalpha,vK,mLabor,nnu,ttheta,mGrid_z_tao(:,1)'); % nK by nShock
    % mProfitAfterTax = (1-mGrid_z_tao(:,2)').*(mY - wage.*mLabor); % nK by nShock

    % mCostAdjust = f_adjust(ddelta,vK,vK',ppsi); % nK by nK_prime
    % mInvest = f_invest(ddelta,vK,vK'); % nK by nK_prime

    mInvest_3D = f_invest(ddelta,mK_3D,mPolicyK);
    mCostAdjust_3D = f_adjust(ddelta,mK_3D,mPolicyK,ppsi);
    mDivident_3D = f_divident(mProfitAfterTax_3D,mInvest_3D,mCostAdjust_3D,mP_3D,mPolicyP,r,ttao); % kPrime by pPrime
    dividentAgg = sum(mDivident_3D.*mDist1,'all');

    consumption = r.*pPrimeAgg + wage.*laborAgg + dividentAgg;
%     consumption_2 = wage/cchi;

%     error = consumption_1 - consumption_2;
    % error_labor = laborAgg - 1/3;
    MarU_consumption = wage./consumption;
    MarU_labor = cchi;
    error = MarU_consumption - MarU_labor;
%     vError(iWage)= error;
    
    display("wage = "+num2str(wage))
    display("marginal utility of consumption = "+num2str(MarU_consumption))
    display("marginal utility of labor = "+num2str(MarU_labor))
    display("diff = "+num2str(error))
    
    if error > 0  % decrease wage
        wage_right = wage;       
    elseif error < 0  % increase wage
        wage_left = wage; %  
    end
end

% figure
% plot(vWage,vError)
% xlabel('wage')
% ylabel('MarU_c - MarU_labor')

if abs(error)<=tol 
    display("Equilibrium wage found. Wage = "+num2str(wage)+". Error = "+ num2str(error))
    display("Labor Demand/Supply = "+num2str(laborAgg));
else
    display("Equilibrium wage not found.")
end

toc

%% Plot figures
% make grid for k and p
mK = repmat(vK,1,nP)';

% Value Function
figure
subplot(1,3,1)
h1=surf(mK,mP,mValue1(:,:,1)','EdgeColor', 'black', 'FaceColor', [255,100,0]/255, 'FaceAlpha', .5, 'Marker', '.' );
hold on
h2=surf(mK,mP,mValue1(:,:,2)','EdgeColor', 'black', 'FaceColor', [1,255,200]/255, 'FaceAlpha', .9, 'Marker', '.');
% legend([h1, h2], {'high tax', 'low tax'});
title('low productivity (z)')
xlabel('k')
ylabel('p')
zlabel('value')
xlim([min(vK),max(vK)])
ylim([0,max(mP,[],'all')])
% zlim([min(mValue1,[],'all'),max(mValue1,[],'all')])

subplot(1,3,2)
h1=mesh(mK,mP,mValue1(:,:,3)','EdgeColor', 'black', 'FaceColor', [255,100,0]/255, 'FaceAlpha', .5, 'Marker', '.');
hold on
h2=mesh(mK,mP,mValue1(:,:,4)','EdgeColor', 'black', 'FaceColor', [1,255,200]/255, 'FaceAlpha', .9, 'Marker', '.');
legend([h1, h2], {'high tax', 'low tax'},'Location','best');
title('medium productivity (z)')
xlabel('k')
ylabel('p')
zlabel('value')
xlim([min(vK),max(vK)])
ylim([0,max(mP,[],'all')])
% zlim([min(mValue1,[],'all'),max(mValue1,[],'all')])

subplot(1,3,3)
h1=mesh(mK,mP,mValue1(:,:,5)','EdgeColor', 'black', 'FaceColor', [255,100,0]/255, 'FaceAlpha', .5, 'Marker', '.');
hold on
h2=mesh(mK,mP,mValue1(:,:,6)','EdgeColor', 'black', 'FaceColor', [1,255,200]/255, 'FaceAlpha', .9, 'Marker', '.');
% legend([h1, h2], {'high tax', 'low tax'});
title('high productivity (z)')
xlabel('k')
ylabel('p')
zlabel('value')
xlim([min(vK),max(vK)])
ylim([0,max(mP,[],'all')])
% zlim([min(mValue1,[],'all'),max(mValue1,[],'all')])

savefig('fig_value')

% Policy Function - k
figure
subplot(1,3,1)
h1=surf(mK,mP,mPolicyK(:,:,1)','EdgeColor', 'black', 'FaceColor', [255,100,0]/255, 'FaceAlpha', .5, 'Marker', '.' );
hold on
h2=surf(mK,mP,mPolicyK(:,:,2)','EdgeColor', 'black', 'FaceColor', [1,255,200]/255, 'FaceAlpha', .9, 'Marker', '.');
% legend([h1, h2], {'high tax', 'low tax'});
title('low productivity (z)')
xlabel('k')
ylabel('p')
zlabel('k policy')
xlim([min(vK),max(vK)])
ylim([0,max(mP,[],'all')])
% zlim([min(mPolicyK,[],'all'),max(mPolicyK,[],'all')/2])

subplot(1,3,2)
h1=mesh(mK,mP,mPolicyK(:,:,3)','EdgeColor', 'black', 'FaceColor', [255,100,0]/255, 'FaceAlpha', .5, 'Marker', '.');
hold on
h2=mesh(mK,mP,mPolicyK(:,:,4)','EdgeColor', 'black', 'FaceColor', [1,255,200]/255, 'FaceAlpha', .9, 'Marker', '.');
legend([h1, h2], {'high tax', 'low tax'},'Location','best');
title('medium productivity (z)')
xlabel('k')
ylabel('p')
zlabel('k policy')
xlim([min(vK),max(vK)])
ylim([0,max(mP,[],'all')])
% zlim([min(mPolicyK,[],'all'),max(mPolicyK,[],'all')/2])

subplot(1,3,3)
h1=mesh(mK,mP,mPolicyK(:,:,5)','EdgeColor', 'black', 'FaceColor', [255,100,0]/255, 'FaceAlpha', .5, 'Marker', '.');
hold on
h2=mesh(mK,mP,mPolicyK(:,:,6)','EdgeColor', 'black', 'FaceColor', [1,255,200]/255, 'FaceAlpha', .9, 'Marker', '.');
% legend([h1, h2], {'high tax', 'low tax'});
title('high productivity (z)')
xlabel('k')
ylabel('p')
zlabel('k policy')
xlim([min(vK),max(vK)])
ylim([0,max(mP,[],'all')])
% zlim([min(mPolicyK,[],'all'),max(mPolicyK,[],'all')])

savefig('fig_policy_k')


% Policy Function - p
figure
subplot(1,3,1)
h1=surf(mK,mP,mPolicyP(:,:,1)','EdgeColor', 'black', 'FaceColor', [255,100,0]/255, 'FaceAlpha', .5, 'Marker', '.' );
hold on
h2=surf(mK,mP,mPolicyP(:,:,2)','EdgeColor', 'black', 'FaceColor', [1,255,200]/255, 'FaceAlpha', .9, 'Marker', '.');
% legend([h1, h2], {'high tax', 'low tax'});
title('low productivity (z)')
xlabel('k')
ylabel('p')
zlabel('p policy')
xlim([min(vK),max(vK)])
ylim([0,max(mP,[],'all')])
% zlim([min(mPolicyP,[],'all'),max(mPolicyP,[],'all')/2])

subplot(1,3,2)
h1=mesh(mK,mP,mPolicyP(:,:,3)','EdgeColor', 'black', 'FaceColor', [255,100,0]/255, 'FaceAlpha', .5, 'Marker', '.');
hold on
h2=mesh(mK,mP,mPolicyP(:,:,4)','EdgeColor', 'black', 'FaceColor', [1,255,200]/255, 'FaceAlpha', .9, 'Marker', '.');
legend([h1, h2], {'high tax', 'low tax'},'Location','best');
title('medium productivity (z)')
xlabel('k')
ylabel('p')
zlabel('p policy')
xlim([min(vK),max(vK)])
ylim([0,max(mP,[],'all')])
% zlim([min(mPolicyP,[],'all'),max(mPolicyP,[],'all')/2])

subplot(1,3,3)
h1=mesh(mK,mP,mPolicyP(:,:,5)','EdgeColor', 'black', 'FaceColor', [255,100,0]/255, 'FaceAlpha', .5, 'Marker', '.');
hold on
h2=mesh(mK,mP,mPolicyP(:,:,6)','EdgeColor', 'black', 'FaceColor', [1,255,200]/255, 'FaceAlpha', .9, 'Marker', '.');
% legend([h1, h2], {'high tax', 'low tax'});
title('high productivity (z)')
xlabel('k')
ylabel('p')
zlabel('p policy')
xlim([min(vK),max(vK)])
ylim([0,max(mP,[],'all')])
% zlim([min(mPolicyP,[],'all'),max(mPolicyP,[],'all')])

savefig('fig_policy_p')

%% Plot the distribution

% make grid for k and p
mK = repmat(vK,1,nP)';

figure
subplot(1,3,1)
h1=surf(mK,mP,mDist1_4D(:,:,1,1)','EdgeColor', 'black', 'FaceColor', [255,100,0]/255, 'FaceAlpha', .5, 'Marker', '.' );
hold on
h2=surf(mK,mP,mDist1_4D(:,:,2,1)','EdgeColor', 'black', 'FaceColor', [1,255,200]/255, 'FaceAlpha', .9, 'Marker', '.');
% legend([h1, h2], {'high tax', 'low tax'});
title('low productivity (z)')
xlabel('k')
ylabel('p')
zlabel('prob')

subplot(1,3,2)
h1=surf(mK,mP,mDist1_4D(:,:,1,2)','EdgeColor', 'black', 'FaceColor', [255,100,0]/255, 'FaceAlpha', .5, 'Marker', '.' );
hold on
h2=surf(mK,mP,mDist1_4D(:,:,2,2)','EdgeColor', 'black', 'FaceColor', [1,255,200]/255, 'FaceAlpha', .9, 'Marker', '.');
legend([h1, h2], {'high tax', 'low tax'},'Location','best');
title('medium productivity (z)')
xlabel('k')
ylabel('p')
zlabel('prob')

subplot(1,3,3)
h1=surf(mK,mP,mDist1_4D(:,:,1,3)','EdgeColor', 'black', 'FaceColor', [255,100,0]/255, 'FaceAlpha', .5, 'Marker', '.' );
hold on
h2=surf(mK,mP,mDist1_4D(:,:,2,3)','EdgeColor', 'black', 'FaceColor', [1,255,200]/255, 'FaceAlpha', .9, 'Marker', '.');
% legend([h1, h2], {'high tax', 'low tax'});
title('high productivity (z)')
xlabel('k')
ylabel('p')
zlabel('prob')

savefig('fig_distribution')

save result

