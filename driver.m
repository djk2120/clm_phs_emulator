clear
close all

kmax  = 7e-9;
krmax = 1e-9;
lai   = 4;
ck    = 3.95;
p50   = -2.5e5;
z     = 25;
soillayers = [     0  , 0.0200  , 0.0600  , 0.1200  , 0.2000  , 0.3200  , 0.4800,...
    0.6800  , 0.9200  , 1.2000  , 1.5200  , 1.8800  , 2.2800  , 2.7200,...
    3.2600  , 3.9000  , 4.6400  , 5.4800  , 6.4200  , 7.4600  , 8.6000]';
rai=5*[0,2.73e-2,3.96e-2,5.02e-2,7.02e-2,...
    8.49e-2,9.36e-2,9.62e-2,9.36e-2,8.67e-2,...
    7.68e-2,6.54e-2,5.36e-2,4.67e-2,3.67e-2,...
    2.62e-2,1.71e-2,1.03e-2,5.70e-3,2.92e-3]';
dz = soillayers(2:end)-soillayers(1:end-1);


ns    = length(soillayers)-1;
smp   = zeros(ns,1)-10000;
q     = 1e-4;

[smp,~,~,hk] = bucket( smp,0*smp,dz,0,0 );


x = {kmax,krmax,lai,ck,p50,z,...
    soillayers,smp,q,rai,hk};


qout= zeros(100,1);
vout= zeros(100,1);
for i=1:1000
[vwp,qout(i),qr] = phs(x);
vout(i) = vwp(2);
[smp,~,~,hk] = bucket( smp,qr,dz,0,0 );
x(8)={smp};
x(11)={hk};

if mod(i,10)==0
subplot(1,3,1)
barh(-1:-1:-20,smp)
xlabel('Soil Potential')
ylabel('Soil Layer')
ylim([-21,0])
xlim([-1e5,0])

subplot(1,3,2)
plot(qout(qout>0))
xlabel('Timestep')
ylabel('Transpiration')
xlim([0,1000])
ylim([0,1e-4])

subplot(1,3,3)
plot(vout(vout<0))
xlabel('Timestep')
ylabel('Leaf Potential')
xlim([0,1000])
ylim([-2e5,0])
pause(0.01)
end
end






