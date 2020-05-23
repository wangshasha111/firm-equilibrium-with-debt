


% make grid for k and p
mK = repmat(vK,1,nP)';
vKP = [mK(:),mP(:)];

% Value Function
figure
for iShock = 1:nShock
    temp = mValue1(:,:,iShock)';
    plot3(vKP(:,1),vKP(:,2),temp(:));
    hold on
end
legend('1','2','3','4','5','6')
hold off;
title('Value Function')
xlabel('k - capital')
ylabel('p - debt')
zlabel('value')
savefig('fig_value')

% Policy Function - k
figure
for iShock = 1:nShock
    temp = mPolicyK(:,:,iShock)';
    plot3(vKP(:,1),vKP(:,2),temp(:));
    hold on
end
legend('1','2','3','4','5','6')
hold off;
title('Policy Function - k^\prime','Interpreter','latex')
xlabel('k - capital')
ylabel('p - debt')
zlabel('k^\prime','Interpreter','latex')
savefig('fig_policy_k')

% Policy Function - k
figure
for iShock = 1:nShock
    temp = mPolicyP(:,:,iShock)';
    plot3(vKP(:,1),vKP(:,2),temp(:));
    hold on
end
legend('1','2','3','4','5','6')
hold off;
title('Policy Function - p^\prime','Interpreter','latex')
xlabel('k - capital')
ylabel('p - debt')
zlabel('p^\prime','Interpreter','latex')
savefig('fig_policy_p')




% Value Function
figure
for iShock = 1:nShock
    for iK = 1:nK
    plot3(repmat(vK(iK),nP,iShock),mP(:,iK),mValue1(iK,:,iShock)');
    hold on
    end
    hold on
end
hold off;
title('Value Function')
xlabel('k - capital')
ylabel('p - debt')
zlabel('value')

% Policy Function for k - capital
figure
for iShock = 1:nShock
    for iK = 1:nK
    plot3(repmat(vK(iK),nP,iShock),mP(:,iK),mPolicyK(iK,:,iShock)');
    hold on
    end
    hold on
end
hold off;
title('Policy Function - k (capital)')
xlabel('k - capital')
ylabel('p - debt')
zlabel('k^\prime','Interpreter','latex')

% Policy Function for p - debt
figure
for iShock = 1:nShock
    for iK = 1:nK
    plot3(repmat(vK(iK),nP,iShock),mP(:,iK),mPolicyP(iK,:,iShock)');
    hold on
    end
    hold on
end
hold off;
title('Policy Function - p (debt)')
xlabel('k - capital')
ylabel('p - debt')
zlabel('p^\prime','Interpreter','latex')