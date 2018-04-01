function [ethaf,metaf] = ethabp(ethai,Pm,lp,profile,znozzle,Te,tpp,gas,zmax,rmax,nres)
%return the ionization fraction (r,z) after 16.6ns from an initial ionization ethai(r,z)

zmin = -zmax;
rmin = -rmax;

r		= [rmin : (rmax-rmin)/nres : rmax];
z       = [zmin : (zmax-zmin)/nres : zmax]; 
ethaf = zeros(size(ethai,1),size(ethai,2));

for i = 1:size(ethai,2)%callculate diffusion for each z position
    %check initial condition to send it in pde solver, if its 0, we return 0
    Ai = max(ethai(:,i));
    if Ai == 0
        soli = zeros(length(r),1);
        metai = zeros(length(r),1);
    else
        Pz = Press(z(i),0,Pm,lp,profile,znozzle);
        if Pz <= 1 %a too small pressure cause problems to the pde solver, a very small pressure leads to very fast diffusion
            soli = zeros(length(r),1);
            metai = zeros(length(r),1);
        else
            %solve the differential equation until tpp=16.6ns
            [soli,metai] = solvepde(ethai(:,i),Pz,Te,tpp,gas,rmax,nres);
        end
    end
    
    %save the result into the matrix ethaf
    ethaf(:,i) = soli;
    metaf(:,i) = metai;
end
end
