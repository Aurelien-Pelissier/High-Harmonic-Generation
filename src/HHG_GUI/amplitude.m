function [z,Iq] = amplitude(dipq,ethaf,Dp,Pm,lp,profile,q,lambda1,znozzle,zmax,nres,gas)
%solving the maxwell equation to calculate the final harmonic amplitude

%basic constants
c		= 299792458;
mu0     = 1.2566e-6;
eta0 = 376.73031;
w1		= c*2*pi/lambda1;
wq      = q*w1;
T0 = 273;
Kb = 1.38064852e-23;

%defining mesh
zmin	= -zmax*(1-0.01); %to avoid z=0 which cause some divergence problem
zmesh	= [zmin : (zmax-zmin)/nres : zmax];

if lp <= 2*zmax
    nres2 = floor(nres*lp/2/zmax);
    zmeshs = [-lp/2 + lp/20 + znozzle : lp/nres2 : lp/2 - lp/20 + znozzle];
    %the lp/20 is added to avoid calculations at the square profile border
else
    zmeshs = zmesh;
end


%solving diff.eq.
sigma = absorb(lambda1,q,gas);    %absorbtion cross section in m2
% sigma = 0;
E0 = 0;

if profile == 'gauss' | profile == 'squ_g'
    a = @(z)-Press(z,0,Pm,lp,profile,znozzle)/(T0*Kb)*sigma/2*100;  %100 because Press is in mbar and we want in Pa
    b = @(z)i*mu0*q*w1*c/2*Press(z,0,Pm,lp,profile,znozzle)/(T0*Kb)*100*(1-ionizf(ethaf,z,zmesh))*sqrt(dipoq(dipq,znozzle,zmesh))*exp(i*Dphi(Dp,z,zmesh));

    [z,E] = ode23s(@(z,E) a(z)*E + b(z), zmesh, E0);
    
    %for the square profile, we have to go through a complicated process
    %because the discontinuity cause problems in the pde solving.

    
elseif profile == 'squar' 
    if nres2 > 1
        a = @(z)-Press(z,0,Pm,lp,profile,znozzle)/(T0*Kb)*sigma/2*100;  %100 because Press is in mbar and we want in Pa
        b = @(z)i*mu0*q*w1*c/2*Press(z,0,Pm,lp,profile,znozzle)/(T0*Kb)*100*(1-ionizf(ethaf,z,zmesh))*sqrt(dipoq(dipq,znozzle,zmesh))*exp(i*Dphi(Dp,z,zmesh));

        [zs,Es] = ode23s(@(z,E) a(z)*E + b(z), zmeshs, E0); %first solving in the interaction region
    
        frac = (zmax + znozzle-lp/2)/(zmax - znozzle -lp/2);%then reconstruct the full z and E
        res1 = floor((nres-nres2)/(1/frac+1));    
        res2 = nres - res1 - nres2;
    
        z1 = [zmin : abs(zmin+lp/2-lp/20-znozzle)/res1 : -lp/2 + lp/20 + znozzle]; %lp/20 to remove the border
        z2 = [lp/2 - lp/20 + znozzle : abs(lp/2-lp/20 + znozzle - zmax)/res2 : zmax];      
        z = [z1';zs;z2'];
    
        E1 = linspace(0,0,res1+1);
        E2 = linspace(Es(end),Es(end),res2+1);
        E =  [E1';Es;E2'];
    else
        z = zmesh;
        E = linspace(0,0,length(z));        
    end
end
    
Iq = E.^2./(2*eta0)*1e-60; %1e-60 cause its arb.unit anyway

end

function dipo = dipoq(dipq,znozzle,zmesh)
    [~, index] = min(abs(zmesh-znozzle));
    dipo = dipq(index);
end

function Dpz = Dphi(Dp,z,zmesh)
    [~, index] = min(abs(zmesh-z));
    Dpz = Dp(floor(end/2)+1,index);
end

function eta = ionizf(ethaf,z,zmesh)
    [~, index] = min(abs(zmesh-z));
    eta = ethaf(floor(end/2)+1,index);
end



