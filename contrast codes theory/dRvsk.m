%close all
R0 = reflectivitystack(2.7,1.8,0);
N=20;permax = 3;n= linspace(2,3,N);
k=1.7;%linspace(0.5,1.7,N);%per=linspace(-permax,permax,N);
R=zeros(1,N); dRbyR0=zeros(1,N);
for i=1:N
    R(i) = reflectivitystack(n(i),k,0.001);
    dRbyR0(i) = (R(i)-R0)./ R(i);
end
figure
plot(n,dRbyR0)
title('dR/R0');
xlabel('n(real part of graphene ref index)');
ylabel('dR/R0');
grid on
% figure
% plot(100*per,dRbyR0)
% title('dR/R0');
% xlabel('dk(change in imaginary part of Si as percent)');
% ylabel('dR/R0');
% grid on