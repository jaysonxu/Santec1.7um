clear;clc
close all;
tic
tmp=dir('.\*StruC0.dcm');
for Nfile=1:size(tmp,1), fprintf(['[%d]\t',tmp(Nfile).name(1:end-5),'\n'],Nfile);end
%%
for Nfile=1:size(tmp,1)
    filename=tmp(Nfile).name;
    STR=double(squeeze(dicomread(filename)));
    nX=size(STR,2);nY=size(STR,3); nZ=size(STR,1);
    A_mean=mean(STR,3);
    All_mean(:,Nfile)=mean(A_mean,2);
    f1=figure(1);plot(All_mean(:,Nfile)); 
    print(f1,[filename,'.png'],'-dpng');
end
    %%
xfactor=1.5e3/166;   
Depth=(1:nZ)*xfactor;
FS=20;
LW=2;
f1=figure(1);plot(Depth,All_mean(:,1),'r',Depth,All_mean(:,2),'b','Linewidth',LW);  %plot(1:500); %
legend('VSCSEL 1.3um','Santec 1.7um'); legend boxoff;
set(gca,'fontsize',FS);
xlabel('Depth (um)');
ylabel('Intensity');
axis tight;
print(f1,['All.png'],'-dpng');

%%
All_Norm=zeros(size(All_mean));
All_Norm(:,1)=(All_mean(:,1)-min(All_mean(:,1)))/(max(All_mean(:,1))-min(All_mean(:,1)));
All_Norm(:,2)=(All_mean(:,2)-min(All_mean(:,2)))/(max(All_mean(:,2))-min(All_mean(:,2)));

f2=figure(2);plot(Depth,All_Norm(:,1),'r',Depth,All_Norm(:,2),'b','Linewidth',LW);  %plot(1:500); %
legend('VSCSEL 1.3um','Santec 1.7um'); legend boxoff;
set(gca,'fontsize',FS);
xlabel('Depth (um)');
ylabel('Intensity (Norm.)');
axis tight;
print(f2,['All_Norm.png'],'-dpng');
toc
%%
fprintf('Test1');
%%
fprintf('Test2'); 
%%
fprintf('Test3'); 
%%
fprintf('Test4'); 
%%
<<<<<<< HEAD
fprintf('Test8'); 
%%
fprintf('Test9'); 
=======
fprintf('Test5'); 
%%
fprintf('Test6'); 
%%
fprintf('Test7'); 
>>>>>>> Jayson
%%
fprintf('Test10'); 