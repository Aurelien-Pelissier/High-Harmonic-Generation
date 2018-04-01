function [KZ,PhiZ,ethai,ethaf] = phase_matching_1D(I0,V,Pmax,lp,profile,gas,q,f,R0,lambda1,tp,Te,alpha,zmax,nres,znozzle,approx)

%calculate phase matching and ionization fraction

zmin	= -zmax*(1-0.01); %to avoid z=0 which cause some divergence problem
Pmin  = 0;
t = 0;
r=0;
rmax = 50e-6;

z	  = [zmin : (zmax-zmin)/nres : zmax];
Pm    = [Pmin : (Pmax-Pmin)/nres : Pmax];

    
%----------------- taking the ionization fraction file -------------------%
    ethai = detectetha(I0,tp,R0,lambda1,gas,zmax,rmax,nres);

%----------- calculating phase matching for each point (Pm,z) ------------%

PhiZ = zeros(length(Pm),length(z));
KZ = zeros(length(Pm),length(z));


    for j = 1:length(Pm)
        if approx == 0
            ethaf = detectethan(I0,tp,V,Pm(j),lp,profile,znozzle,Te,f,R0,lambda1,gas,zmax,rmax,nres,ethai);
        else
            K = sqrt(Pm(j)/1050)*7*exp(-(V-100)/150)+1;
            ethaf = K.*ethai;
        end
        for i = 1:length(z)
            [PhiZ(j,i),KZ(j,i)] = phase(z(i),r,t,I0,tp,Pm(j),lp,profile,znozzle,q,alpha,lambda1,R0,gas,zmax,rmax,nres,ethaf);
        end
    end

end
