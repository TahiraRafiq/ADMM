XFM = Wavelet('Daubechies',4,4);

meu=10e-20;

iteration=2;


neta=zeros(Nx*3,Ny,iteration);
ini=zeros(Nx,Ny,iteration);
ini(:,:,1)=SRSoS;

for m=1:iteration
    full_coeff_WT=XFM*ini(:,:,m);
thre_wt= 0.9;
shrink_WT=zeros(Nx,Ny);
for j=1:Ny
for i=1:Nx
if (full_coeff_WT(i,j)+neta(i,j,m))>thre_wt
   shrink_WT(i,j)=(((abs(full_coeff_WT(i,j)+neta(i,j,m)))-thre_wt)*(full_coeff_WT(i,j)+neta(i,j,m))/abs(full_coeff_WT(i,j)+neta(i,j,m)));
else
    shrink_WT(i,j)=0;
end
end
end

% figure,
% imshow(abs(shrink_WT),[])
%% TV
full_coeff_TV=TVOP*ini(:,:,m);

% TVX

thre_TVx=0.09;
shrink_TVx=zeros(Nx,Ny);
for j=1:Ny
for i=1:Nx
 if (full_coeff_TV(i,j,1)+neta(i+256,j,m))>thre_TVx
%  if sqrt((full_coeff_TV(i,j,1).^2)+(neta(i+256,j,m).^2))>thre_TVx
%   shrink_TVx(i,j)=(abs(full_coeff_TV(i,j,1))-thre_TVx)*(full_coeff_TV(i,j,1)/abs(full_coeff_TV(i,j,1)));
  shrink_TVx(i,j)=((sqrt((full_coeff_TV(i,j,1).^2)+(neta(i+256,j,m).^2)))-thre_TVx)*(full_coeff_TV(i,j,1)+neta(i+256,j,m))./(sqrt((full_coeff_TV(i,j,1).^2)+(neta(i+256,j,m).^2)));
else
    shrink_TVx(i,j)=0;
end
end
end

% figure,
% imshow(abs(shrink_TVx),[])

% TVY
shrink_TVy=zeros(Nx,Ny);
for j=1:Ny
for i=1:Nx
  if (full_coeff_TV(i,j,2)+neta(i+512,j,m))>thre_TVx
%  if sqrt((full_coeff_TV(i,j,2).^2)+(neta(i+512,j,m).^2))>thre_TVx
  % shrink_TVy(i,j)=(abs(full_coeff_TV(i,j,2))-thre_TVx)*(full_coeff_TV(i,j,2)/abs(full_coeff_TV(i,j,2)));
     shrink_TVy(i,j)=((sqrt((full_coeff_TV(i,j,2).^2)+(neta(i+512,j,m).^2)))-thre_TVx)*(full_coeff_TV(i,j,2)+neta(i+512,j,m))./(sqrt((full_coeff_TV(i,j,2).^2)+(neta(i+512,j,m).^2)));
else
    shrink_TVy(i,j)=0;
end
end
end



%% step 3
U1=[shrink_WT;shrink_TVx;shrink_TVy];

%% Step 4
ini(:,:,1)=a;

b=zeros(Nx,Ny,11);
b(:,:,1)=0;
p=ini(:,:,m);
r=zeros(Nx,Ny,11);
r(:,:,1)=ini(:,:,m);


for i=1:10
    rr=TVOP*p;
R2=rr(:,:,1);
R3=rr(:,:,2);
R1=XFM*p;
R=[R1;R2;R3];
last_braket=U1-R-neta(:,:,m);
RH_U11=(XFM'*last_braket(1:256,:));
U21=zeros(Nx,Ny,2);
U21(:,:,1)=last_braket(257:512,:);
U21(:,:,2)=last_braket(513:end,:);
RH_U21=(TVOP'*U21);
RH_all_u=meu*(RH_U11)+meu*(RH_U21);
    
    for n=1:Nc
       Ewithp(:,:,n)=c_sens(:,:,n).*p;
    end
    k_space=fftshift(fft2(fftshift(Ewithp)));
    k_sample=k.*k_space;
   
    Eh_half=ifftshift(ifft2(ifftshift(k_sample)));
    Eh=conjugate_sens.*Eh_half;
      q=Eh;
      qsum=zeros(Nx,Ny);
    for n=1:Nc
        qsum=qsum+q(:,:,n);
    end
%      qsum=-qsum-(10000.*RH_all_u);
     qsum=-qsum-(RH_all_u);
    
permute_r=permute(r,[1 2 3]);
reshape_r=reshape(permute_r,256*256,[]);
    alpha=reshape_r(:,i)'*reshape_r(:,i)./(p(:)'*qsum(:));
    b(:,:,i+1)=b(:,:,i)+alpha*p;
    r(:,:,i+1)=r(:,:,i)-alpha*qsum;
    beta=reshape_r(:,i+1)'*reshape_r(:,i+1)./(reshape_r(:,i)'*reshape_r(:,i));
    p=r(:,:,i+1)+beta*p;
end


ini(:,:,m+1)=b(:,:,11);
%% step 5
rr=TVOP*ini(:,:,m+1);
R2=rr(:,:,1);
R3=rr(:,:,2);
R1=XFM*ini(:,:,m+1);
R=[R1;R2;R3];
neta(:,:,m+1)=neta(:,:,m)-(U1-R);
end
  
figure,
imshow(abs(ini(:,:,2)),[])

error = (ref-abs(ini(:,:,1)) ).^2;
    RMSE = sqrt(sum(sum(error))/(Nx * Ny));
    NRMSE = RMSE/(Nx*Ny)
% %     




iteration=2;
sca=zeros(1,iteration);

for j=1:iteration
error = (ref-abs(ini(:,:,j)) ).^2;
    RMSE = sqrt(sum(sum(error))/(Nx * Ny));
    NRMSE = RMSE/(Nx*Ny)
    sca(1,j)=NRMSE;
end
j=[1:iteration];
figure,
plot(j,sca)
xlabel('Iteration')
ylabel('Error')