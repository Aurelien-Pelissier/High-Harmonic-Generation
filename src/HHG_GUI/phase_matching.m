function [KZ,PhiZ,ethai,ethaf,ethafb,ethaft] = phase_matching(I0,V,Pm,lp,profile,znozzle,gas,q,f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t)

%calculate and phase matching and ionization fraction

zmin	= -zmax*(1-0.01); %to avoid z=0 which cause some divergence problem
rmin	= -rmax;



z		= [zmin : (zmax-zmin)/nres : zmax];
r		= [rmin : (rmax-rmin)/nres : rmax];

    ethai = detectetha(I0,tp,R0,lambda1,gas,zmax,rmax,nres,tp); %detecting and calculating ionization fraction after one pulse, t = 2*tp is like t=inf
    if V/f >= 2.5*R0 %if the ions leave the focus area between 2 pulses
        ethaf = ethai;
        ethafb = zeros(size(ethai,1),size(ethai,2));
        ethaft = detectetha(I0,tp,R0,lambda1,gas,zmax,rmax,nres,t);
    else
        [ethaf,ethafb,ethaft] = detectethan(I0,tp,V,Pm,lp,profile,znozzle,Te,f,R0,lambda1,gas,zmax,rmax,nres,ethai,t); %detecting and calculating ionization fraction after last pulse
    end
    
    %------ calculating phase matching for each point (z,r)-------%

    PhiZ = zeros(length(r),length(z));
    KZ = zeros(length(r),length(z));

    for l = 1:length(z)*length(r) %to avoid a double loop "for"
        i = floor(l/(length(z)+0.0000000000001))+1;
        j = mod(l,length(r))+1;
        [PhiZ(j,i),KZ(j,i)] = phase(z(i),r(j),t,I0,tp,Pm,lp,profile,znozzle,q,alpha,lambda1,R0,gas,zmax,rmax,nres,ethaft);
    end

    
end