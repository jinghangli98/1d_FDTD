%Computer project 2
clear; close all; clc;
dt = 2.5e-11;
dz = 1.5e-2;
T = 1e-9;
t = 0:dt:6*T;
Ex = zeros(5000,1000);
Ex(1:length(t),1) = exp(-(t-3*T).^2/T^2);
c = 299792458; %m/s
k = 1e3; %spatial grid points
%% 
mu = 1.257e-6 * ones(1,k); %(H/m)
e0 = 8.85e-12;
er = e0;
er = 7.21 * e0; %(F/m)
% er = 5.92 * e0;
% er = 20.093 * e0;
sigma_star = 0;
% sigma_star = 0.03687;
% sigma_star = 0.31684;

ep = [e0*ones(1,500), er*ones(1,100), e0*ones(1,400)];
sigma = [0*ones(1,500), sigma_star*ones(1,100), 0*ones(1,400)];

%% update coefficient
B = -(dt/dz)./(mu); 
B = repmat(B,5000,1);
C = -2*dt./(dz*(2*ep + dt*sigma));
C = repmat(C,5000,1);
D = ((2*ep-dt*sigma)./(2*ep+dt*sigma));  
D = repmat(D,5000,1);


%% 
Hy = zeros(5000,k);
for nt = 1:4900
   for nz = 1:998
       Ex(nt+1, nz+1) = D(nt+1,nz)*Ex(nt,nz+1) + C(nt,nz)*(Hy(nt,nz+1)-Hy(nt,nz));    
       Hy(nt+1, nz) = Hy(nt,nz) + B(nt,nz)*(Ex(nt+1,nz+1) - Ex(nt+1,nz));     
       Ex(nt,1000) = 0;
       Hy(nt,1000) = 0;
   end
   nz = 999;
   Ex(nt+1, nz+1) = 0;
   Hy(nt+1, nz) = Hy(nt,nz) + B(nt,nz)*(0 - Ex(nt+1,nz));     
end

counter = 1;
for i = [480,501,560,620]
    figure(1)
    set(gcf,'position',[10, 10,2000,500])
    subplot(1,4,counter)
    counter = counter + 1;
    plot(linspace(22.5,62.5,1600),Ex(900:2500-1,i))
    grid on
    xlabel("time (nano seconds)")
    ylabel("Amplitude")
    title(sprintf("At index: %d", i))
    axis([20,65,-1,1])
end

playVideo = "true";
if (playVideo == "true")
    figure(2)
    for i = 1:5000
        hold off
        plot(Ex(i,:),"b")
        hold on
        plot(377*Hy(i,:),"r")
        grid on
        axis([0,1000,-2,2])
        getframe();
    end
end




