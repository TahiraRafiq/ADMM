iteration_alp2=5;
x=zeros(Nx,Ny,iteration_alp2+1);
x(:,:,1)=SRSoS;
U2=zeros(Nx,Ny,iteration_alp2+1);
U2(:,:,1)=x(:,:,1);
U0=zeros(Nx*Nc,Ny,iteration_alp2+1);
meu=10e-10;
v1=1;
v2=1;
cg_ite=4;
 Threshold_WT=2.5e-7;
Threshold_TV=1e-7;
%% neta
neta_twenty=zeros(Nx*Nc,Ny,iteration_alp2+1);
neta_twentyone=zeros(Nx*3,Ny,iteration_alp2+1);
neta_twentytwo=zeros(Nx,Ny);


%%

for j=1:iteration_alp2
%% U0
for n=1:Nc
SX(:,:,n)=c_sens(:,:,n).*x(:,:,j);
end

permute_SX=permute(SX,[1 3 2]);
reshape_SX=reshape(permute_SX,Nx*Nc,[]);

permute_d=permute(d,[1 3 2]);
reshape_d=reshape(permute_d,Nx*Nc,[]);


Uzeros_cg=zeros(Nx*Nc,Ny,cg_ite+1);
Uzeros_cg(:,:,1)=reshape_SX;
p=reshape_SX;
r=zeros(Nx*Nc,Ny,cg_ite+1);
r(:,:,1)=reshape_SX;



for i=1:cg_ite
   
    data_cons=(ifftshift(ifft2(ifftshift(reshape_d-(fftshift(fft2(fftshift(p))))))));
    q=-data_cons+meu*(p-reshape_SX-neta_twenty(:,:,j));
    
permute_r=permute(r,[1 2 3]);
reshape_r=reshape(permute_r,Nx*Ny*Nc,[]);
    alpha=reshape_r(:,i)'*reshape_r(:,i)./(p(:)'*q(:));
    Uzeros_cg(:,:,i+1)=Uzeros_cg(:,:,i)+alpha*p;
    r(:,:,i+1)=r(:,:,i)-alpha*q;
    beta=reshape_r(:,i+1)'*reshape_r(:,i+1)./(reshape_r(:,i)'*reshape_r(:,i));
    p=r(:,:,i+1)+beta*p;
end

 U0(:,:,j+1)=Uzeros_cg(:,:,cg_ite+1);
%% U1

WT_full=XFM*x(:,:,j);
    summation_of_WT_and_neta=WT_full+neta_twentyone(1:256,:,j);
    shrink_WT=zeros(Nx,Ny);
for k=1:Nx
for i=1:Ny
if (summation_of_WT_and_neta(k,i))>Threshold_WT
   shrink_WT(k,i)=((abs(summation_of_WT_and_neta(k,i)))-Threshold_WT)*(summation_of_WT_and_neta(k,i))/(abs(summation_of_WT_and_neta(k,i)));
else
    shrink_WT(k,i)=0;
end
end
end


% % TVX
TV_full=TVOP*x(:,:,j);
summation_of_TVX_and_neta=TV_full(:,:,1)+neta_twentyone(257:512,:,j);
shrink_TVx=zeros(Nx,Ny);

for k=1:Nx
for i=1:Ny
 if (summation_of_TVX_and_neta(k,i))>Threshold_TV
  shrink_TVx(k,i)=(abs(summation_of_TVX_and_neta(k,i))-Threshold_TV)*(summation_of_TVX_and_neta(k,i))./abs(summation_of_TVX_and_neta(k,i));
else
    shrink_TVx(k,i)=0;
end
end
end

% % TVY
shrink_TVy=zeros(Nx,Ny);
summation_of_TVY_and_neta=TV_full(:,:,2)+neta_twentyone(513:768,:,j);
for k=1:Nx
for i=1:Ny

if (summation_of_TVY_and_neta(k,i))>Threshold_TV
  
     shrink_TVy(k,i)=(abs(summation_of_TVY_and_neta(k,i))-Threshold_TV)*(summation_of_TVY_and_neta(k,i))./abs(summation_of_TVY_and_neta(k,i));
else
    shrink_TVy(k,i)=0;
end
end
end

U1=[shrink_WT;shrink_TVx;shrink_TVy];


%% U2
Utw0_cg=zeros(Nx,Ny,cg_ite+1);
Utw0_cg(:,:,1)=U2(:,:,1);
p2=U2(:,:,1);
r2=zeros(Nx,Ny,cg_ite+1);
r2(:,:,1)=p2;


for i=1:cg_ite
   
    ptv=TVOP*p2;
    ru2=(U1-[XFM*p2;ptv(:,:,1);ptv(:,:,2)]-neta_twentyone(:,:,j));
    RH_wt=XFM'*ru2(1:256,:);
    RH_TV=zeros(Nx,Ny,2);
    RH_TV(:,:,1)=ru2(257:512,:);
    RH_TV(:,:,2)=ru2(513:768,:);
    R=TVOP'*RH_TV;
    RH=-meu*v1*(RH_wt+R);
    q=RH+meu*v2*(p2-x(:,:,j)-neta_twentytwo(:,:,j));
  
permute_r=permute(r2,[1 2 3]);
reshape_r=reshape(permute_r,Nx*Ny,[]);
    alpha=reshape_r(:,i)'*reshape_r(:,i)./(p2(:)'*q(:));
    Utw0_cg(:,:,i+1)=Utw0_cg(:,:,i)+alpha*p2;
    r2(:,:,i+1)=r2(:,:,i)-alpha*q;
    beta=reshape_r(:,i+1)'*reshape_r(:,i+1)./(reshape_r(:,i)'*reshape_r(:,i));
    p2=r2(:,:,i+1)+beta*p2;
end

 U2(:,:,j+1)=Utw0_cg(:,:,cg_ite+1);


%% x

x_cg=zeros(Nx,Ny,cg_ite+1);
x_cg(:,:,1)=x(:,:,1);
p_x=x(:,:,1);
r_x=zeros(Nx,Ny,cg_ite+1);
r_x(:,:,1)=x(:,:,1);


for i=1:cg_ite
   
    for n=1:Nc
Sx(:,:,n)=c_sens(:,:,n).*p_x;
end

permute_Sx=permute(Sx,[1 3 2]);
reshape_Sx=reshape(permute_Sx,Nx*Nc,[]);


first_bracket=(U0(:,:,j+1)-reshape_Sx-neta_twenty(:,:,j));

permute_first_bracket=permute(first_bracket,[1 2 3]);
reshape_first_bracket=reshape(permute_first_bracket,Nx,Ny,[]);

first=-meu*(conjugate_sens.*reshape_first_bracket);

add=zeros(Nx,Ny);
for n=1:Nc
add=add+first(:,:,n);
end
    
    q=add-meu*v2*(U2(:,:,j+1)-p_x-neta_twentytwo(:,:,j));
    
    
permute_r=permute(r_x,[1 2 3]);
reshape_r=reshape(permute_r,Nx*Ny,[]);
    alpha=reshape_r(:,i)'*reshape_r(:,i)./(p_x(:)'*q(:));
    x_cg(:,:,i+1)=x_cg(:,:,i)+alpha*p_x;
    r_x(:,:,i+1)=r_x(:,:,i)-alpha*q;
    beta=reshape_r(:,i+1)'*reshape_r(:,i+1)./(reshape_r(:,i)'*reshape_r(:,i));
    p_x=r_x(:,:,i+1)+beta*p_x;
end

x(:,:,j+1)= x_cg(:,:,cg_ite+1);


%%
for n=1:Nc
for_neta(:,:,n)=c_sens(:,:,n).*x(:,:,j+1);
end

permute_for_neta=permute(for_neta,[1 3 2]);
reshape_for_neta=reshape(permute_for_neta,Nx*Nc,[]);

neta_twenty(:,:,j+1)=neta_twenty(:,:,j)-(U0(:,:,j)-reshape_for_neta);

for_neta21_wt=XFM*U2(:,:,j+1);
for_neta21_TV=TVOP*U2(:,:,j+1);
neta_twentyone(:,:,j+1)=neta_twentyone(:,:,j)-(U1-[for_neta21_wt;for_neta21_TV(:,:,1);for_neta21_TV(:,:,2)]);
neta_twentytwo(:,:,j+1)=neta_twentytwo(:,:,j)-(U2(:,:,j+1)-x(:,:,j+1));

end

figure,
imshow(abs(x(:,:,iteration_alp2+1)),[])



scalar=zeros(1,iteration_alp2);

j=[1:iteration_alp2];
for n=1:iteration_alp2

gfdg=x(:,:,n);
min_a = min(gfdg(:));
max_a =(max(gfdg(:)));

    normhjgh = (gfdg-min_a)./abs(max_a-min_a); 
    
gfdg=ref;
min_a = min(gfdg(:));
max_a =(max(gfdg(:)));

    normhjghnbvnv = (gfdg-min_a)./abs(max_a-min_a);
    
error = (abs(normhjghnbvnv)-abs(normhjgh)).^2;
    RMSE = sqrt(sum(sum(error))/(Nx * Ny));
    NRMSE = RMSE/(Nx*Ny);
    scalar(1,n)=NRMSE;
end
 figure,
% plot(2:50,scal(:,2:50))
% hold on
%  plot(2:50,scala(:,2:50))
% hold on
 plot(j,scalar)
xlabel('Iteration')
ylabel('Error')