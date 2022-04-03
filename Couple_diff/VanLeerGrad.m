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
function pGradF=VanLeerGrad(U,dt,dx,a)
N=length(U);
theta=zeros(N,1);
phi=zeros(N,1);

Flux=zeros(N,1);
pGradF=zeros(N,1);
if a>0
    mu=a*dt/dx;
    UR=[U(2:N)];
    UL=[U(1:N-1)];
    theta=max(-1000,min(1000,(U-UL)./(UR-U)));
    phi=(abs(theta)+theta)./(1+abs(theta));
    Flux=a*U+0.5*(1-mu)*a*phi.*(UR-U);
    pGradF=(Flux-[Flux(N);Flux(1:N-1)])/dx;
end

if a<0
    mu=a*dt/dx;
    UR1=[U(2:N);U(1)];
    UR2=[U(3:N);U(1);U(2)];
    theta(1:N)=max(-1000,min(1000,(UR2-UR1)./(UR1-U)));
    phi(1:N)=(abs(theta(1:N))+theta(1:N))./(1+abs(theta(1:N)));
    Flux(1:N)=a*UR1+0.5*(-1-mu)*a*phi(1:N).*(UR1-U);
    pGradF(1:N)=(Flux(1:N)-[Flux(N);Flux(1:N-1)])/dx;
end

return