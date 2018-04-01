function [etha,dethaz,dethar] = ioniz(z,r,zmax,rmax,nres,ethaft);
%return Ionization fraction and its derivative at position (z,r) and time t%

zmin = -zmax;
rmin = -rmax;


%-------------- calculating the ionization rate on (r,z) --------------%
etha = valetha(z,r,zmax,rmax,ethaft);

%-------------- calculating derivative dethaz --------------%
zstep = (zmax-zmin)/nres;
ethaz1 = valetha(z+zstep,r,zmax,rmax,ethaft);
dethaz = (ethaz1 - etha)/zstep;

%-------------- calculating derivative dethar --------------%
rstep = (rmax-rmin)/nres;
ethar1 = valetha(z,r+rstep,zmax,rmax,ethaft);
dethar = (ethar1 - etha)/rstep;
dethar=0;
end

function ethafrz = valetha(z,r,zmax,rmax,ethaf)
%assign value to etha(r,z) using the matrix ethaf
%it works only if nres is a multiple of size(etha).


zmin = -zmax;
rmin = -rmax;


fr = floor((r-rmin)/(rmax-rmin)*size(ethaf,1))+1;
fz = floor((z-zmin)/(zmax-zmin)*size(ethaf,2))+1;
if fr > size(ethaf,1) %to avoid dimension problem on the last index
   fr = size(ethaf,1);
end
if fz > size(ethaf,2)
   fz = size(ethaf,2);
end

ethafrz = ethaf(fr,fz);

end