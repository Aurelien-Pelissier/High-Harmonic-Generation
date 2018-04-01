function [Sw2,PHI,PSI,T,X,K,qW2] = dipole(I0,lambda1,tp,gas,myfile,axes1)
%calculate the dipole response, adapted from TJ's code

%% basic constants
eta0 = 376.73031;
auElec = 5.142e11; %electric field in atomic unit
c		= 299792458;
cau     = 137; %celerity of light in au


tres = floor(3000*tp*1e15/20); %increasing the pulse length will increase the mandatory time resolution.
xres = 3^8;

%% experiment parameter in atomic units

w1		= cau*2*pi/lambda1*5.291e-11; %laser pulsation in au
tau     = tp/2.419e-17; %pulse length in au
E0 = sqrt(I0*1e4*2*eta0)/auElec; %driving field in au

X0 = potde(gas); %potential depth

%% defining everything to avoid loading SS_AN.mat 
% !! it works only for the ground state
% load SS_AN.mat



MAX_X = 2*E0/(w1^2) *2;
% X_CUT is used to set the boundary conditions to avoid interferences between short and long trajectory
% the potential is not the same for all gas, so we set different boundaries.
if gas == 'Xe' 
    XCUT = MAX_X/2.5;
elseif gas == 'Kr'
    XCUT = MAX_X/3.5;
elseif gas == 'Ar'
    XCUT = MAX_X/4.5;
end
Nx = xres;
x = 160*linspace(-1,1,Nx); dx = x(2)-x(1);
dk = 2*pi/(Nx*dx);
k = -Nx/2*dk+dk/2:dk:Nx/2*dk-dk/2;


% setting boundary conditions
Vi = 0*x;
for countx = 1:Nx
    if x(countx)>=XCUT
        Vi(countx) = (x(countx)-XCUT).^2;
    end
    if x(countx)<=-XCUT
        Vi(countx) = (x(countx)+XCUT).^2;
    end
end

Vi = Vi*5e-4;
Vr = -1./sqrt(x.^2+X0^2);
Vc = Vr - 1i*Vi;



%%
%TJ defined tau this way
% TAU = 2*pi/w1;
% tau = 5*TAU/(sqrt(2*log(2)));


% Nx = length(x);
Nt = tres   +1;
t = 2*linspace(-tau,tau,Nt); dt = t(2)-t(1);
dw = 2*pi/(Nt*dt);
w = -Nt/2*dw+dw/2:dw:Nt/2*dw-dw/2;

%w = 2*pi/(2*dt)*linspace(-1,1,Nt);

% env = exp(-(t-2.5*tau).^2/tau^2);
sig = tau/(2*sqrt(2*log(2)));
env = exp(-t.^2/sig.^2);

NC = 1;
CEPvec = 0*linspace(0,2*pi,NC+1) + pi/4;
CEPvec = CEPvec(1:end-1);

%%
for countC = 1:NC
%     CEP = CEPvec(countC);
    CEP = 1.15*pi;
    relphase = 0.4;
    wph = 2*pi/(1/(11.7*1.602e-19/(6.626e-34))/2.41888e-17);
    El = E0*env.*(exp(1i*(w1*t-CEP))+0.0*exp(2*1i*(w1*t-CEP+relphase*pi)));
%     El = 0.1*E0*env.*exp(1i*(wph*t));
%     figure(1)
%     plot(t,El)
    
    Elmat(countC,:) = El;
    
    %%
    YYN = 1/(2*X0)*sech(x./X0);
    Y0 = YYN(1,:);
    normY = sum(abs(Y0.^2))*dx;
    psi0A  = Y0/sqrt(normY);
    psi0 = psi0A;
    
    %% For the momentum
    raise1 = circshift(eye(Nx),[0,1]);
    lower1 = circshift(eye(Nx),[0,-1]);
    
    dUsqbdrsq = (raise1*psi0' + lower1*psi0' - 2*psi0')/(dx^2); % this is accurate to dr^2; you don't need to worry about dUsqbdrsq(1) because U0 = rvec*R so U0(1)=0 always
    dUsqbdrsq = dUsqbdrsq';
    dUsqbdrsq(1) = 0;
    dUsqbdrsq(Nx) = 0;
    
    HU0 = -1/2*dUsqbdrsq + Vr.*psi0;
    
    %% Finding the energy
    
    Egs_sum = sum(conj(psi0).*HU0)*dx;
    Egs_sum*27.211;
    
    %%
    skipper = 5;
    phimat = zeros(round(Nt/skipper),Nx);
    psimat = zeros(round(Nt/skipper),Nx);
    D = zeros(size(Y0,1),round(Nt/10));
    dip = zeros(1,Nt);
    tot = zeros(1,Nt);
    tic
    
    for countT = 1:Nt
        %     countT
        Vl = -real(El(countT)).*x;
        Va = Vc + Vl;
        phi0 = fftshift(fft(fftshift(psi0.*exp(-i*Va*dt)))); %*(-1i)
        psi0 = ifftshift(ifft(ifftshift(phi0.*exp(-i*k.^2*1/2*dt)))); %*(-1i)
        
        if mod(countT,10)==0
            for countNY = 1:size(YYN,1)
                D(countNY,countT/10) = sum(conj(YYN(countNY,:)).*psi0)*dx;
            end
        end
        
        if mod(countT,skipper)==0
            phimat(countT/skipper,:) = abs(phi0).^2;
            psimat(countT/skipper,:) = abs(psi0).^2;
            
            if mod(countT,floor(Nt/5))==0
                axes(axes1)
                imagesc(x.*5.29e-2,t.*2.418e-2,log10(abs((1/sqrt(5.29e-11)).*psimat)))
                colormap(axes1,jet)
                axis xy
                caxis([-5,4])
                drawnow
            end
            
        end
        
        dip(countT) = (psi0.*x)*psi0'*dx;
        
        tot(countT) = sum(conj(psi0).*psi0)*dx;
        
    end
    toc
%     close all
    
    %calculating harmonic amplitude dipole response
    dipa = (circshift(dip, [1 1]) + circshift(dip, [-1 -1]) - 2*dip)/dt^2;
    smoother = exp(-t.^20/(t(end)/1.25).^20);
    dipa = real(dipa).*smoother;
    Sw = abs(fftshift(fft(fftshift(dipa)))).^2;
    Sw2 = Sw(floor(length(Sw)/2):floor(7*length(Sw)/8));
    
    %back to SI
    K = k./5.29e-11;
    T = t.*2.418e-17;
    X = x.*5.29e-11;
    qW = w.*27.211;
    qW2 = qW(floor(length(qW)/2):floor(7*length(qW)/8));
    PSI = 1/sqrt(5.29e-11).*psimat;
    PHI = 1/sqrt(5.29e-11).*phimat;
    
    
    %write the resulting dipole amplitude to a file
    fid = fopen(myfile,'wt');
    for j=1:size(Sw2,1)
        fprintf(fid, '%d ', Sw2(j,:));
        fprintf(fid, '\n');
    end
    fclose(fid);
   
    
    
    
    %checking plot (not mandatory)
    
%     figure(2);
%     imagesc(x.*5.29e-2,t.*2.418e-2,log10(abs((1/sqrt(5.29e-11)).*psimat)))
%     title('Dipole response in space domain')
%     xlabel('x[nm]'); ylabel('t[fs]')
%     h = colorbar;
%     ylabel(h,'log(|Psi|^2)');
%     title('Dipole response in momentum domain')
%     colormap(jet)
%     caxis([1,9])
%     xlim([-1e11 1e11])
%     axis xy
    
%     figure(3);
%     imagesc(k./5.29e-11,t.*2.418e-2,log10(abs((1/sqrt(5.29e-11)).*phimat)))
%     xlabel('k[m-1]'); ylabel('t[fs]')
%     h = colorbar;
%     ylabel(h,'log(|Phi|^2)');
%     title('Dipole response in momentum domain')
%     colormap(jet)
%     caxis([1,9])
%     xlim([-1e11 1e11])
%     axis xy
%     %%
%     
%     
%     % dipe = dip.*exp(-(t-2.5*tau).^8/(1.9*tau).^8);
%     
% %     dipa = (circshift(dip, [1 1]) + circshift(dip, [-1 -1]) - 2*dip)/dt^2;
% %     % dipa(1) = dipa(2); dipa(end)=dipa(end-1);
% %     smoother = exp(-t.^20/(t(end)/1.25).^20);
% %     dipe = real(dip).*smoother;
% %     dipa = real(dipa).*smoother;
%     
%     figure(4)
% %     semilogy(w*27.211,w.^2.*abs(fftshift(fft(dip))),w*27.211,w.^2.*abs(fftshift(fft(dipe))),w*27.211,abs(fftshift(fft(dipa))))
%     %p1 = semilogy(w*27.211,w.^2.*abs(fftshift(fft(dipe))),'linewidth', 2)%,'color',rand(3,1));
%     %hold on
% %     plot(-Egs_mat(1)+Egs_mat,0*Egs_mat+1,'.r')
% %     plot(-Egs_mat(1)+Egs_mat,0*Egs_mat+10,'.r')
% %     dipamat(countC,:) = dipa;
% %     dipemat(countC,:) = dipe;
%     
% %     Cend(countC,:) = D(:,end);
%     
% end
% 
% %%
% % figure(5)
% % % surf(t(1:10:end-1)*2.41888e-17*1e15,2:20,log10(abs(D(2:end,:)).^2)), view([0 0 1]), axis tight, shading interp
% % imagesc(t(1:10:end-1)*2.41888e-17*1e15,2:20,log10(abs(D(2:end,:)).^2)), axis xy
% % ca = caxis;
% % % -Egs_mat(1)+Egs_mat(2:end)
% % caxis([ca(2)-4. ca(2)])
% % xlim([-20 20])
% % set(gca,'fontsize',16)
% % xlabel('Time (fs)')
% % ylabel('State Number')
% % % zlim([ca(2)-4 ca(2)])
% %%
% figure(4)
% wa = 2*pi/(2*dt)*linspace(-1,1,length(t));
% p2 = semilogy(w*27.211,abs(fftshift(fft(fftshift(dipa)))),'linewidth', 2, 'color', 'r');
% legend([p2],'da(w)')
% xlabel('harmonic order q'); ylabel('S(w) [arb. unit]')
% ylim([1e-7 1e1])
% xlim([3 30])
% title('Harmonic amplitude S(w)')
% 
% % figure(7)
% % semilogy(-Egs_mat(1)+Egs_mat(2:end),abs(D(2:end,end)).^2,'.--r')
% %%
% 
% % figure(8)
% % dipW = fftshift(fft(fftshift(dipamat')))';
% % redip = [dipW; dipW; dipW;];
% % imagesc(w*27.211,CEPvec,log10(tj_smooth_2d(abs(redip).^2,3,3))), axis xy, xlim([6 15])
% % ca = caxis;
% % % caxis([ca(2)*0 ca(2)*0.1])
% % caxis([ca(2)-4.0 ca(2)])

end