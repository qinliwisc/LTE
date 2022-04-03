% function pGradF=VanLeerFlux(U,dt,dx,a)
% N=length(U);
% theta=zeros(N,1);
% phi=zeros(N,1);
% 
% pGradF=zeros(N,1);
% if a>0
%     mu=a*dt/dx;
%     theta(2:N-1)=max(-1000,min(1000,(U(2:N-1)-U(1:N-2))./(U(3:N)-U(2:N-1))));
%     phi(2:N-1)=(abs(theta(2:N-1))+theta(2:N-1))./(1+abs(theta(2:N-1)));
%     Flux(2:N-1)=a*U(2:N-1)+0.5*(1-mu)*a*phi(2:N-1).*(U(3:N)-U(2:N-1));
%     pGradF(3:N-1)=(Flux(3:N-1)-Flux(2:N-2))/dx;
% %     sigma(2:N-1)=(U(3:N)-U(2:N-1)).*phi(2:N-1)/dx;
% %     pGradF(3:N-1)=a*(U(3:N-1)-U(2:N-2)+0.5*(1-a*dt/dx)*dx*(sigma(3:N-1)-sigma(2:N-2)))/dx;
% end
% 
% if a<0
%     mu=a*dt/dx;
%     theta(1:N-2)=max(-1000,min(1000,(U(3:N)-U(2:N-1))./(U(2:N-1)-U(1:N-2))));
%     phi(1:N-2)=(abs(theta(1:N-2))+theta(1:N-2))./(1+abs(theta(1:N-2)));
% %     sigma(1:N-2)=(U(2:N-1)-U(1:N-2)).*phi(1:N-2)/dx;
% %     pGradF(1:N-3)=a*(U(2:N-2)-U(1:N-3)-0.5*(1+a*dt/dx)*dx*(sigma(2:N-2)-sigma(1:N-3)))/dx;
%     Flux(1:N-2)=a*U(2:N-1)+0.5*(-1-mu)*a*phi(1:N-2).*(U(2:N-1)-U(1:N-2));
%     pGradF(2:N-2)=(Flux(2:N-2)-Flux(1:N-3))/dx;
% end
% 
% return

%% dirichlet boundary condition
function Flux=VanLeerFlux(U,dt,dx,a)
U = U(:);
N=length(U);
% theta=zeros(N-1,1);
% phi=zeros(N-1,1);

Flux=zeros(N-1,1);
if a>0
    mu=a*dt/dx;
    UR=[U(2:N)];
    UL=[U(1:N-1)];
    theta=max(-10000,min(10000,(UR(1:end-1)-UL(1:end-1))./(UR(2:end)-UL(2:end))));
    theta = [0;theta];
    phi=(abs(theta)+theta)./(1+abs(theta));
    Flux = mu*UL + 0.5*mu*(1-mu)*phi.*(UR-UL);
end

if a<0
    mu=a*dt/dx;
    UR=[U(2:N)];
    UL=[U(1:N-1)];
    theta=max(-10000,min(10000,(UR(2:end)-UL(2:end))./(UR(1:end-1)-UL(1:end-1))));
    theta = [theta;0];
    phi=(abs(theta)+theta)./(1+abs(theta));
    Flux = mu*UR - 0.5*mu*(1+mu)*phi.*(UR-UL);
end

Flux = Flux';
return