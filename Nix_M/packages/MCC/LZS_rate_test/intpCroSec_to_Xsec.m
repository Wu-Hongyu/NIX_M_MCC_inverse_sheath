clear
name=['intpCroSec.mat']; % Ωÿ√Ê ¶“-E
s=load(name);  
sc=struct2cell(s);  
CroSec=cell2mat(sc); 

Xsec = ImportLXCat;
Xsec.name{1}='for test';
Xsec.dir{1}='for test';
Xsec.ela{1}{1}=[CroSec(:,1),CroSec(:,2)];
% Xsec.add_key_points2();
Xsec.info.ela{1}{1}='From LZS HeavyRate';
Xsec.ion{1}{1}=[0,0;1e4,0];
Xsec.ionThresh{1}{1}=0;
Xsec.info.ion{1}{1}='no';
Xsec.exc{1}{1}=[0,0;1e4,0];
Xsec.excThresh{1}{1}=0;
Xsec.info.exc{1}{1}='no';
Xsec.att{1}{1}=[0,0;1e4,0];
% Xsec.attThresh{1}{1}=0;
Xsec.info.att{1}{1}='no';
Xsec.eff{1}{1}=[0,0;1e4,0];
Xsec=Xsec.getEnergy();
Xsec.info.eff{1}{1}='no';
Xsec=Xsec.totalXsection();
figure

Xsec.plotXsections(14,'log','log')