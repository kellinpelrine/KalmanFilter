clear;
clf;
close all;

rng(1);
obserror = 0;
filtererror = 0;
filtererrors = 0;
smootherror = 0;
maerror3 = 0;
maerror5 = 0;
iter = 10000;
for i = 1:iter

    n = 50;
    a = .5;
    %generate errors in X. Note: error variance = 1 here
    v = randn(1, n);
    w = randn(1, n);

    vecX(1) = v(1);
    for i = 2:n
        vecX(i) = a * vecX(i - 1) + v(i);
    end

    for i = 1:n
        vecY(i) =  vecX(i) + w(i);
    end

    Xhat(1) = 0;
    P(1) = (1 - a^2)^(-1);   
    Xhatt(1) = P(1) * (P(1) + 1)^(-1) * (vecY(1));  
    Pt(1) = P(1) - P(1) * (P(1) + 1)^(-1) * P(1);

    for i = 2:n
        Xhat(i) = a*Xhatt(i - 1);
        P(i) = a* Pt(i - 1) * a + 1;
        Xhatt(i) = Xhat(i) + P(i) *(P(i) + 1)^(-1) * (vecY(i) - Xhat(i));
        Pt(i) = P(i) - P(i)*(P(i) + 1)^(-1) * P(i);
    end
    
    obserror = obserror + sum((vecX - vecY).^2)/n;
    filtererror = filtererror + sum((vecX - Xhat).^2)/n;
    filtererrors = filtererrors + sum((vecX(1:n - 1) - Xhat(2:n)).^2)/n;
    
    for i = 1:(n - 1)
        J(i) = Pt(i) * a * P(i + 1) ^ (-1);
    end
    Xs = zeros(1, n);
    Xs(n) = Xhatt(n);
    for i = n - 1:-1:1
        Xs(i) = Xhatt(i) + J(i) * (Xs(i + 1) - Xhat(i + 1));
    end
    smootherror = smootherror + sum((vecX - Xs).^2)/n;
    maerror3 = maerror3 + sum((movmean(vecY,3) - vecX).^2) /n;
    maerror5 = maerror5 + sum((movmean(vecY,5) - vecX).^2) /n;

end

obserror = obserror/iter;
filtererror = filtererror/iter;
filtererrors = filtererrors/iter;
smootherror = smootherror/iter;
maerror3 = maerror3/iter;
maerror5 = maerror5/iter;

% % for i = 2:n + 1
% %     if i ~= tau + 1
% %         partK = P(:, :, i - 1)*H(:,:,i - 1)*inv( H(:,:,i - 1)' * P(:, :, i - 1) *H(:,:,i - 1));
% %         Xhatt(:,:,i - 1) = Xhat(:, :, i - 1) + F(:, :, i - 1)* partK * (vecY(i - 1) - H(:, :, i - 1)' * Xhat(:,:, i - 1));
% %         Xhat(:, :,i) = F(:, :, i - 1)*Xhatt(:, :, i - 1);
% %         Pt(:,:,i - 1) = P(:,:,i - 1) - partK * H(:,:,i - 1)' *P(:,:,i - 1);
% %         P(:,:,i) = F(:,:,i - 1)* Pt(:,:,i - 1) * F(:, :, i - 1)' + [1 0; 0 0];           
% %     else
% %         Xhatt(:,:,i - 1) = Xhat(:, :, i - 1);
% %         Xhat(:, :,i) = F(:, :, i - 1)*Xhat(:, :, i - 1);
% %         Pt(:,:, i - 1) = P(:, :, i - 1);
% %         P(:,:,i) = F(:,:,i - 1)*  P(:,:,i - 1) * F(:,:, i - 1)' + [1 0; 0 0];
% %     end
% % end
% Xhat = vecY;
% for i = tau1:tau2
%     Xhat(i) = a * Xhat(i - 1);
% end

figure('Name','Blue: state. Red: observations. Green: filter.')

%plot(Xhat, 'color','black');
hold on
plot(vecX(1:n),'color','blue')
plot(vecY(1:n), 'color','red')
plot(Xhat(1:n), 'color','green')

figure('Name','Blue: state. Red: observations. Black: smoother.')
hold on
plot(vecX(1:n),'color','blue')
plot(vecY(1:n), 'color','red')
plot(Xs(1:n), 'color', 'black')

% 
% for i = 1:tau1 - 1
%     P(i) = Q;
%     Pt(i) = 0;
% end
% for i = tau1:tau2
%     P(i) = a^2 * Pt(i - 1) + Q;
%     Pt(i) = P(i);
% end
% for i = tau2 + 1:n    
%     P(i) = Q;
%     Pt(i) = 0;
% end
% P(tau2 + 1) = a^2 * Pt(tau2) + Q; %note: above loop just fills a 0 for P(tau2 + 1). This is correct value.
% 
% for i = 1:n
%     J(i) = 0;
% end
% for i = tau1:tau2
%     J(i) = a*Pt(i)*(1/P(i + 1));
% end
% 
% Xhat2 = Xhat;
% Xhat2(tau2 + 1) = a*Xhat(tau2); %prediction of next available point from last missing point
% 
% Xs = Xhat;
% for i = tau2:-1:tau1
%     Xs(i) = Xhat2(i) + J(i)*(Xs(i + 1) - Xhat2(i + 1));
% end
% 
% plot(yplot(tau1 + 1:tau2 + 1),Xs(tau1:tau2),'color','green');

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