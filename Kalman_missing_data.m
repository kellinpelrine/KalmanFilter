clear;
clf;
close all;
rng(3);

n = 50;
a = .7;
tau = 20; %missing entries tau1 to tau2
tau1 = 20;
tau2 = 29; 
yplot = [0:51];
%generate errors in X
v = randn(n,1);
Q = 1; %error variance. NOTE: need to modify above function if changed 

vecX(:,1) = v(1);
for i = 2:n
    vecX(:,i) = a * vecX(:,i - 1) + v(i);
end

for i = 1:n
    H(:, i) = 1;
end
for i = tau1:tau2
    H(:, i) = 0;
end
H = permute(H, [2 1]);

for i = 1:n
    vecY(i) = H(i, 1)' * vecX(:, i);
end

% Xhat(:,:,1) = 0';
% P(:,:,1) = 1; %note: error variance = 1 here

% for i = 2:n + 1
%     if i ~= tau + 1
%         partK = P(:, :, i - 1)*H(:,:,i - 1)*inv( H(:,:,i - 1)' * P(:, :, i - 1) *H(:,:,i - 1));
%         Xhatt(:,:,i - 1) = Xhat(:, :, i - 1) + F(:, :, i - 1)* partK * (vecY(i - 1) - H(:, :, i - 1)' * Xhat(:,:, i - 1));
%         Xhat(:, :,i) = F(:, :, i - 1)*Xhatt(:, :, i - 1);
%         Pt(:,:,i - 1) = P(:,:,i - 1) - partK * H(:,:,i - 1)' *P(:,:,i - 1);
%         P(:,:,i) = F(:,:,i - 1)* Pt(:,:,i - 1) * F(:, :, i - 1)' + [1 0; 0 0];           
%     else
%         Xhatt(:,:,i - 1) = Xhat(:, :, i - 1);
%         Xhat(:, :,i) = F(:, :, i - 1)*Xhat(:, :, i - 1);
%         Pt(:,:, i - 1) = P(:, :, i - 1);
%         P(:,:,i) = F(:,:,i - 1)*  P(:,:,i - 1) * F(:,:, i - 1)' + [1 0; 0 0];
%     end
% end
Xhat = vecY;
for i = tau1:tau2
    Xhat(i) = a * Xhat(i - 1);
end


%plot(Xhat, 'color','black');
figure('Name','Blue: state. Red: missing data. Green: filter. Black: smoother.')
hold on
plot(vecX(1:tau1 - 1),'color','blue')
plot(yplot(tau2 + 2:n + 1),vecX(tau2 + 1:n),'color','blue')
plot(yplot(tau1 + 1:tau2 + 1),vecX(1,tau1:tau2),'color','red');
plot(yplot(tau1 + 1:tau2 + 1),Xhat(1, tau1:tau2),'color','green')

for i = 1:tau1 - 1
    P(i) = Q;
    Pt(i) = 0;
end
for i = tau1:tau2
    P(i) = a^2 * Pt(i - 1) + Q;
    Pt(i) = P(i);
end
for i = tau2 + 1:n    
    P(i) = Q;
    Pt(i) = 0;
end
P(tau2 + 1) = a^2 * Pt(tau2) + Q; %note: above loop just fills a 0 for P(tau2 + 1). This is correct value.

for i = 1:n
    J(i) = 0;
end
for i = tau1:tau2
    J(i) = a*Pt(i)*(1/P(i + 1));
end

Xhat2 = Xhat;
Xhat2(tau2 + 1) = a*Xhat(tau2); %prediction of next available point from last missing point

Xs = Xhat;
for i = tau2:-1:tau1
    Xs(i) = Xhat2(i) + J(i)*(Xs(i + 1) - Xhat2(i + 1));
end

plot(yplot(tau1 + 1:tau2 + 1),Xs(tau1:tau2),'color','black');

% for i = tau:n
%     J(:,:,i) = Pt(:,:,i) *F(:,:,i)' * inv(P(:,:,i + 1));
% end
% 
% S(:,:,n + 1) = Xhat(:,:,n + 1);
% for i = 1:n 
%     S(:,:,n + 1 - i) = Xhatt(:,:,n + 1 - i) + J(:,:,n + 1 - i) * (S(:,:,n + 2 - i) - Xhat(:,:,n + 2 - i));
% end
% 
% %plot(S(1,:));
% 
% 
% r(:,n + 1) = [0 0]';
% for i = n:-1:tau + 1
%     r(:,i) = H(:,:,i)*inv(H(:,:,i)' * P(:,:,i + 1) * H(:,:,i)) ...
%           + (F(:,:,i) - F(:,:,i)*P(:,:,i + 1)*H(:,:,i)*inv(H(:,:,i)' * P(:,:,i + 1) * H(:,:,i))*H(:,:,i)')*r(:,i + 1);
% end
% r(:,tau) = F(:,:,tau) * r(:,tau + 1);
% for i = tau - 1:-1:1
%     r(:,i) = H(:,:,i)*inv(H(:,:,i)' * P(:,:,i + 1) * H(:,:,i)) ...
%           + (F(:,:,i) - F(:,:,i)*P(:,:,i + 1)*H(:,:,i)*inv(H(:,:,i)' * P(:,:,i + 1) * H(:,:,i))*H(:,:,i)')*r(:,i + 1);
% end
% for i = 2:n + 1
%     S2(:,i) = Xhat(i) + P(:,:,i)*r(:,i - 1);
% end
% 
% %plot(S2(1,:));
% %vecY = [0 vecY'];
% %plot(vecY);



%%%%%%%%%%%%%%%%OLD CODE%%%%%%%%%%%%%%%%%

%set F to vary correctly with time
% for i = 1:tau1 - 1
%     F(:, :, i) = [a 0; 0 1];
% end
% for i = tau1:tau2
%     F(:, :, i) = [a 0; 1 0];
% end
% for i = tau2+1:n
%     F(:, :, i) = [a 0; 0 1];
% end
% 
% for i = 1:n
%     vecv(:,:,i) = [v(i) 0];
% end
% vecv = permute(vecv, [2 1 3]);
% 
% vecX(:,:,1) = [0 0]';
% for i = 1:n
%     vecX(:,:,i + 1) = F(:,:,i) * vecX(:, :, i) + vecv(:, :, i);
% end
% 
% for i = 1:n
%     H(:, :, i) = [1 0];
% end
% H(:, :, tau) = [0 0];
% H = permute(H, [2 1 3]);
% 
% for i = 1:n
%     vecY(i) = H(:, :, i)' * vecX(:, :, i + 1);
% end
% vecY = vecY';
% 
% Xhat(:,:,1) = [0 0]';
% P(:,:,1) = [1 0; 0 0]; %note: error variance = 1 here
% 
% for i = 2:n + 1
%     if i ~= tau + 1
%         partK = P(:, :, i - 1)*H(:,:,i - 1)*inv( H(:,:,i - 1)' * P(:, :, i - 1) *H(:,:,i - 1));
%         Xhatt(:,:,i - 1) = Xhat(:, :, i - 1) + F(:, :, i - 1)* partK * (vecY(i - 1) - H(:, :, i - 1)' * Xhat(:,:, i - 1));
%         Xhat(:, :,i) = F(:, :, i - 1)*Xhatt(:, :, i - 1);
%         Pt(:,:,i - 1) = P(:,:,i - 1) - partK * H(:,:,i - 1)' *P(:,:,i - 1);
%         P(:,:,i) = F(:,:,i - 1)* Pt(:,:,i - 1) * F(:, :, i - 1)' + [1 0; 0 0];           
%     else
%         Xhatt(:,:,i - 1) = Xhat(:, :, i - 1);
%         Xhat(:, :,i) = F(:, :, i - 1)*Xhat(:, :, i - 1);
%         Pt(:,:, i - 1) = P(:, :, i - 1);
%         P(:,:,i) = F(:,:,i - 1)*  P(:,:,i - 1) * F(:,:, i - 1)' + [1 0; 0 0];
%     end
% end
% 
% plot(vecX(1,:),'color','blue')
% hold on
% plot(yplot(tau1 + 2:tau2 + 2),vecX(1,tau1 + 1:tau2 + 1),'color','red');
% %plot (Xhat(1, :))
% 
% for i = tau:n
%     J(:,:,i) = Pt(:,:,i) *F(:,:,i)' * inv(P(:,:,i + 1));
% end
% 
% S(:,:,n + 1) = Xhat(:,:,n + 1);
% for i = 1:n 
%     S(:,:,n + 1 - i) = Xhatt(:,:,n + 1 - i) + J(:,:,n + 1 - i) * (S(:,:,n + 2 - i) - Xhat(:,:,n + 2 - i));
% end
% 
% %plot(S(1,:));
% 
% 
% r(:,n + 1) = [0 0]';
% for i = n:-1:tau + 1
%     r(:,i) = H(:,:,i)*inv(H(:,:,i)' * P(:,:,i + 1) * H(:,:,i)) ...
%           + (F(:,:,i) - F(:,:,i)*P(:,:,i + 1)*H(:,:,i)*inv(H(:,:,i)' * P(:,:,i + 1) * H(:,:,i))*H(:,:,i)')*r(:,i + 1);
% end
% r(:,tau) = F(:,:,tau) * r(:,tau + 1);
% for i = tau - 1:-1:1
%     r(:,i) = H(:,:,i)*inv(H(:,:,i)' * P(:,:,i + 1) * H(:,:,i)) ...
%           + (F(:,:,i) - F(:,:,i)*P(:,:,i + 1)*H(:,:,i)*inv(H(:,:,i)' * P(:,:,i + 1) * H(:,:,i))*H(:,:,i)')*r(:,i + 1);
% end
% for i = 2:n + 1
%     S2(:,i) = Xhat(i) + P(:,:,i)*r(:,i - 1);
% end
% 
% %plot(S2(1,:));
% %vecY = [0 vecY'];
% %plot(vecY);