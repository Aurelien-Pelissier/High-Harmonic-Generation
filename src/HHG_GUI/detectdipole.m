function [dipqz,dipqt,qW,dipo0] = detectdipole(I0,lambda1,znozzle,tp,R0,gas,q,t,zmax,nres,axes1)
%detecting dipole file to avoid doing the calculation everytime%
%if the file doesnt exist, a new one is created%
%doing this for every intensity to calculate dipq(z)

%creating the folder if doesnt exist
myfolder = fullfile(sprintf('%.0fnm', lambda1*1e9),'dipole',gas,sprintf('%.0ffs',tp*1e15),'Intensity_files');
if exist (myfolder,'dir')  == 0 %checking if the folder exist
    mkdir(myfolder) %create the folder if doesnt exist
end


%defining zmesh
zmin	= -zmax;
z	= [zmin : (zmax-zmin)/nres : zmax];
[~, indexz] = min(abs(z-znozzle));

%defining tmesh
tmesh	= [-tp : (2*tp)/nres : tp];
[~, indext] = min(abs(t-tmesh));

%finding the index corresponding on the harmonic
qW = harm_axis(tp,lambda1);
[~, indexq] = min(abs(qW-q));

%calculating dipole response for every positions z and every time t
for j = 1:length(tmesh)
for i = 1:length(z)
    
    if j == indext
        I = Igauss(I0,tp,R0,lambda1,z(i),0,tmesh(indext));
        myfile = fullfile(sprintf('%.0fnm', lambda1*1e9),'dipole',gas,sprintf('%.0ffs',tp*1e15),'Intensity_files',sprintf('dipole_tp%.0ffs_I0%.2e.txt',tp*1e15, I));
        if exist(myfile,'file') ~= 0 %checking if the file exist
            dipo = dlmread(myfile); %reading the ionization matrix
            %disp(sprintf('using dipole_tp%.0ffs_I0%.2e.txt',tp*1e15, I)) % too much disp make the program slower
        else
            disp(sprintf('calculating new dipole file_%.0ffs_I0%.2e.txt',tp*1e15, I))
            dipo = dipole(I,lambda1,tp,gas,myfile,axes1); %creating a new one
        end
        a = floor(length(dipo)/120);
        dipqz(i) = max(dipo(indexq-a:indexq+a));
%     dipq(i) = dipo(index);
%     dipq(i) = mean(dipo(indexq-a:indexq+a));
    end
    
      if i == indexz
        I = Igauss(I0,tp,R0,lambda1,z(indexz),0,tmesh(j));
        myfile = fullfile(sprintf('%.0fnm', lambda1*1e9),'dipole',gas,sprintf('%.0ffs',tp*1e15),'Intensity_files',sprintf('dipole_tp%.0ffs_I0%.2e.txt',tp*1e15, I));
        if exist(myfile,'file') ~= 0 %checking if the file exist
            dipo = dlmread(myfile); %reading the ionization matrix
            %disp(sprintf('using dipole_tp%.0ffs_I0%.2e.txt',tp*1e15, I))  % too much disp make the program slower
        else
            disp(sprintf('calculating new dipole file_%.0ffs_I0%.2e.txt',tp*1e15, I))
            dipo = dipole(I,lambda1,tp,gas,myfile,axes1); %creating a new one
        end
        a = floor(length(dipo)/120);
        dipqt(j) = max(dipo(indexq-a:indexq+a));
%     dipq(i) = dipo(index);
%     dipq(i) = mean(dipo(indexq-a:indexq+a));
    end

    if i == indexz
        if j == indext
            dipo0=dipo;
        end
    end
    
end
end

end


function qW = harm_axis(tp,lambda1)
    tau     = tp/2.419e-17;
    Nt = 3000*tp*1e15/20   +1;
    t = 2*linspace(-tau,tau,Nt); dt = t(2)-t(1);
    dw = 2*pi/(Nt*dt);
    w = -Nt/2*dw+dw/2:dw:Nt/2*dw-dw/2;    
    qWl = w.*27.211*61/72*lambda1/1050e-9; %*61/72 dunno why but need it, its to make sure that everymaxium correspond to an odd harmonic
    qW = qWl(floor(length(qWl)/2):floor(7*length(qWl)/8));
end