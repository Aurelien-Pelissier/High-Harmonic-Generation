function [solf,metaf] = solvepde(ethai,Pm,Tev,tpp,gas,rmax,nres)
%look matlab pdesolver on internet for more information on pde


%return the ionization fraction(r) and the metastable atoms fraction after tpp considering velocity, diffucion and recombination
%solf = ionization fraction after tpp, array (r)
%metaf = metastable atoms fraction after tpp, array (r)

%ethai = imitial ionization fraction, array(r)
%Pm = maximal pressure
%TeV = electron temperature in eV
%tpp = time between 2 pulses
%rmax and nres = mappimg imformations

%initial condition
P0 = 1013; %mbar
Natm = Ndens(gas);
etha0 = max(ethai); %ionization at r=0
rho0 = etha0*Natm*Pm/P0; %electron density after one pulse (e-/m3)

%recombination
kb = 1.38064852e-23; %boltzman cst
Te = Tev/kb/6.242e18; %free electron temperature in K
gamma0 = 1.1e-20*(Te)^(-9/2)*rho0^2*1e-9; %we add 1e-9 because we work in ns, if we dont do that it doesnt work

%diffusion
e  = 1.602e-19;
Ti = 300; %temperature of ions
TXev = Ti*kb*6.242e18;%temperature of ions in eV
Diatm = mobility(gas)*kb*Ti/e*1e4*1e-9;%cm2/s  %we add 1e-9 because we work in ns, if we dont do that it doesnt work

%defining mesh
rmin	= -rmax;
tmin	= 0;
tmax	= tpp*1e9; %time between 2 pulses %we add 1e9 because we work in ns, if we dont do that the integration fail

xmesh = [rmin : (rmax-rmin)/nres : rmax];
tmesh = [tmin : (tmax-tmin)/200 : tmax];


% solving pde
m = 0;
sol = pdepe(m,@(x,t,u,DuDx)pdex1pde(x,t,u,DuDx,gamma0,Diatm,Tev,TXev,Pm),@(x)pdex1ic(x,xmesh,ethai),@pdex1bc,xmesh,tmesh);
u = sol(:,:,1);

solf = u(end,:);


% ions which recombined into metastable xenon

% -------- if we are interrested by the evolution of metastable xenon
% rhometa = zeros(length(xmesh),length(tmesh));
% 
% for i = 1:length(tmesh)
%     ti = tmesh(1:i);
%     for j = 1:length(xmesh)
%         if i==1
%             rhometa(j,i) = 0;
%         else
%             uj = u(:,j);
%             rhometa(j,i) = gamma0*trapz(ti,(uj(1:i)).^3);
%         end
%     end
% end
% 
% metaf = rhometa(:,end);

% -------- if we are interrested only by the metastable xenon at t=tpp
metaf = zeros(length(xmesh),1);
for j = 1:length(xmesh)
    uj = u(:,j);
    if max(uj) == 0
        metaf(j) = 0;
    else
        metaf(j) = gamma0*trapz(tmesh,uj.^3);
    end
end




% ------------------ surface plot (not mandatory) ----------------- %
% up = permute(u,[2 1]);
% figure;
% surf(tmesh*1e9,xmesh*1e6,up);
% title('Diffusion & Recombination');
% xlabel('Time t [ns]');
% ylabel('Distance r [um]');
% zlabel('ion density [m-3]');
% colorbar
% zlim([0,1]);

% figure;
% surf(tmesh*1e9,xmesh*1e6,rhometa);
% title('metastable xenon');
% xlabel('Time t [ns]');
% ylabel('Distance r [um]');
% zlabel('Xe* density [m-3]');
% colorbar
% zlim([0,1]);
% ------------------------------------------------------------------ %

function [c,f,s] = pdex1pde(x,t,u,DuDx,gamma0,Diatm,Tev,TXev,Pm)
P0 = 1013; %mbar

    %If we suppose the pressure stays constant during the process
    Di = Diatm*P0/Press(0,x,Pm,1,'gauss',0);%we need only yhte pressure on r axis, so ze dont care about the z-parameters
 
    % If we consider that the pressure change with the diffusion, and the ionization fraction is modified by the decay
    % rhom = P(x,Pm) * Natm/P0;
    % rhodecay = (rho0x(x,rho0) - ((2^(1/2)*(1/(gamma0*t + 1/(2*rho0x(x,rho0)^2)))^(1/2))/2))/3; %The decay is actually less than that since diffusion lead to a slower decay, expecially for, but couldnt do beter approximaition
    % rhoat = rhom - rho0x(x,rho0) +rhodecay;
    % rhotot = rhoat + u;
    % Di = Diatm*Natm/rhotot;
    
Da = Di*(1+Tev/TXev);



c = 1;
f = Da*1e-4*DuDx;
s = -gamma0*u^3;

% ------------------------- fct for pde solver ------------------------- %

function u0 = pdex1ic(x,xmesh,ethai)
[~, index] = min(abs(xmesh-x));
u0 = ethai(index);


function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur;
qr = 0;


% --------------------------------------------------------------------------

% this 2 fction are for the complicated Di
function eth = etha(x,etha0)
    R0 = 10e-6;
    eth = etha0*exp(-2*x^2/R0^2);

function rox = rho0x(x,rho0)
    R0 = 10e-6;
    rox = rho0*exp(-2*x^2/R0^2);   
% rhodecay = gamma0*integral(@(t)u^3,0,t);
