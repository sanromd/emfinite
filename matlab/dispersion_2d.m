% Simulating the propagation of a pulse inside a waveguide 2D

% Computational Domain
Dist_a=500.0e-9; % Horizantal Distance
Dist_b=300.0e-9; % Vertical Distance
time=2*Dist_a/3e8; % Time

dx=5.0e-9;     % Discretizing Space
dy=5.0e-9;
dt=1/(3e8*sqrt((1/dx)^2+(1/dy)^2));    %Discretizing Time
dt2=dt^2;

Nt=time/dt;        %Number of points in X axis
Nx=round(Dist_a/dx);     %Number of points in X axis
Ny=round(Dist_b/dy);
Nd=2;

% Material Properties
pi=atan(1.0)*4.0 ;
epo=8.854e-12; %farad per meter
muo=4*pi*1e-7; %Hennery per meter


%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c--------------------- PML Walls Setup
%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c--------------------- PML Walls from all sides
delta=8; %Number of PML Cells
jb=Ny-delta; ib=Nx-delta; ja=delta; ia=delta;

sgdy=zeros(Nx,Ny); sgdx=zeros(Nx,Ny); sgby=zeros(Nx,Ny); sgbx=zeros(Nx,Ny); 
stdy=zeros(Nx,Ny); stdx=zeros(Nx,Ny); stby=zeros(Nx,Ny); stbx=zeros(Nx,Ny);
cdx1=zeros(Nx,Ny); cdx2=zeros(Nx,Ny); cdy1=zeros(Nx,Ny); cdy2=zeros(Nx,Ny);

cbzx1=zeros(Nx,Ny); cbzx2=zeros(Nx,Ny); cbzy1=zeros(Nx,Ny); cbzy2=zeros(Nx,Ny);

%c--------------------- PML Sigmas every where " Related to PML "
sgmax=-2*3e8*log10(1e-8);
sgdy_max=sgmax/(2*delta*dy);
sgby_max=sgdy_max;
sgdx_max=sgmax/(2*delta*dx);
sgbx_max=sgdx_max;

%---------------------	calculatig sigmas
for i=1:Nx
    for j=jb:Ny  % Top side
        sgdy(i,j)=sgdy_max*((j-jb)/delta)^3;
        sgby(i,j)=sgby_max*((j-jb+0.5)/delta)^3;      
    end
end

for i=1:Nx;
    for j=1:ja;	% Bottom side
        sgdy(i,j)=sgdy_max*((ja-j+1)/delta)^3;
        sgby(i,j)=sgby_max*((ja-j+0.5)/delta)^3;        
    end
end

for i=ib:Nx;
    for j=1:Ny;	%Right side
        sgdx(i,j)=sgdx_max*((i-ib)/delta)^3;
        sgbx(i,j)=sgbx_max*((i-ib+0.5)/delta)^3;
    end
end

for i=1:ia;
    for j=1:Ny;	 %Left side
        sgdx(i,j)=sgdx_max*((ia-i+1)/delta)^3;
        sgbx(i,j)=sgbx_max*((ia-i+0.5)/delta)^3;
    end
end

%--------------------- TE Loop Constant for PML Walls

for i=1:Nx;
    for j=1:Ny;        
        stdy(i,j)=dt*sgdy(i,j);
        stdx(i,j)=dt*sgdx(i,j);
        stby(i,j)=dt*sgby(i,j);
        stbx(i,j)=dt*sgbx(i,j);
                
        cdx1(i,j)=(1-stdy(i,j))/(1+stdy(i,j));
        cdx2(i,j)=dt/(dy*(1+stdy(i,j)));
        
        cdy1(i,j)=(1-stdx(i,j))/(1+stdx(i,j));
        cdy2(i,j)=dt/(dx*(1+stdx(i,j)))*(-1);
        
        cbzx1(i,j)=(1-stbx(i,j))/(1+stbx(i,j));
        cbzx2(i,j)=dt/(dx*(1+stbx(i,j)))*(-1);
        
        cbzy1(i,j)=(1-stby(i,j))/(1+stby(i,j));
        cbzy2(i,j)=dt/(dy*(1+stby(i,j)));
        
    end
end

% Specifiyng material properties along the x axis
mus=zeros(Nx,Ny);
eps=zeros(Nx,Ny);
eps(:,:)=epo;
mus(:,:)=muo;

%Source parameters
Freq=3e9;  % oscillation frequency
T0=3e-10;   % Gaussian pulse t0
TT=5e-11;   % Gaussian pulse waist

% preparing matricies for calculations
Hz=zeros(Nx,Ny);
Ex=zeros(Nx,Ny);
Ey=zeros(Nx,Ny);

Py=zeros(Nx,Ny,Nd);
Px=zeros(Nx,Ny,Nd);

Py1=zeros(Nx,Ny,Nd); Px1=zeros(Nx,Ny,Nd); Py2=zeros(Nx,Ny,Nd); Px2=zeros(Nx,Ny,Nd);

C1=zeros(Nx,Ny,Nd);
C2=zeros(Nx,Ny,Nd);
C3=zeros(Nx,Ny,Nd);

Dex=zeros(Nx,Ny);
Dey=zeros(Nx,Ny);

Bhzx=zeros(Nx,Ny);
Bhzy=zeros(Nx,Ny);
Bhz=zeros(Nx,Ny);

strctr=zeros(Nx,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ZnOC1 = 2.81418;
ZnOC2 = 0.87968;
ZnOC3 = 0.3042^2;
ZnOC4 = 0.00711 ;


epf=2.81418;
ZnOC3=ZnOC3/(36*pi^2*1e28);
ZnOC4=ZnOC4*36*pi^2*1e28;

ac(1)=ZnOC2*epo;
ac(2)=ZnOC4*epo;

bc(1)=1;
bc(2)=0;


dc(1)=ZnOC3;
dc(2)=1;


cc(1)=0;
cc(2)=0;



for n=1,2;
    C1_ZnO(n)=(4*dc(n)-2*bc(n)*dt2)/(2*dc(n)+cc(n)*dt);
    C2_ZnO(n)=(cc(n)*dt-2*dc(n))/(2*dc(n)+cc(n)*dt);
    C3_ZnO(n)=(2*dt2*ac(n))/(2*dc(n)+cc(n)*dt);
    ep_ZnO=epo*epf;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NoW=1;
Wlength=200e-9	;   
WThickness=100e-9;

WTdy=round(0.5*WThickness/dy);
WLdx=round(0.5*Wlength/dx);

i_p=zeros(1,NoW);
j_p=zeros(1,NoW);
jtop=zeros(1,NoW);
jbot=zeros(1,NoW);
irgs=zeros(1,NoW);
ilfs=zeros(1,NoW);

for i=1:NoW;
i_p(i)=round(0.5*Nx)
j_p(i)=round(i*Ny/(NoW+1))
jtop(i)=j_p(i)+WTdy;
jbot(i)=j_p(i)-WTdy;

irgs(i)=i_p(i)+WLdx;
ilfs(i)=i_p(i)-WLdx;


strctr(i_p(i),j_p(i))=150;
end 

 


    for k=1:NoW;
        for  i=ilfs(k):irgs(k)
            for  j=jbot(k):jtop(k)
               
                
                strctr(i,j)=100;
                for n=1,2;
                    C1(i,j,n)=C1_ZnO(n);
                    C2(i,j,n)=C2_ZnO(n);
                    C3(i,j,n)=C3_ZnO(n);
                    eps(i,j)=ep_ZnO	;
                end
            end
        end
    end


%Hz(Nx/2,Ny/2)=20/377;
Rec=zeros(1,Nt);
b=0; % a variable for showing the movie it need to be zero .


figure(1)
mesh(strctr')
figure(2)
% 
% 
% 
for n=1:Nt;       % Start of the FDTD loop in Time
    
    for i=1:Nx-1;
        for j=2:Ny-1;
            Dex(i,j)=cdx1(i,j)*Dex(i,j)+cdx2(i,j)*(Hz(i,j)-Hz(i,j-1));
            
            Sumx=sum(Px(i,j,:));
            Ex(i,j)=(Dex(i,j)-Sumx)/eps(i,j);
            
            for k=1:2;
            Px(i,j,k)=C1(i,j,k)*Px1(i,j,k)+C2(i,j,k)*Px2(i,j,k)+C3(i,j,k)*Ex(i,j);
            
            Px2(i,j,k)=Px1(i,j,k);
            Px1(i,j,k)=Px(i,j,k);
            end
        end
    end
     
   % Ex(Nx/5,Ny/2)=Ex(Nx/5,Ny/2)+(20/377)*exp(-((n*dt-T0)/TT)^2);
    
    
    %---------------finding Ey
    for i=2:Nx-1;
        for j=1:Ny;
            Dey(i,j)=cdy1(i,j)*Dey(i,j)+cdy2(i,j)*(Hz(i,j)-Hz(i-1,j));
            
            Sumy=sum(Py(i,j,:));
            Ey(i,j)=(Dey(i,j)-Sumy)/eps(i,j);
            
            for k=1:2;
            Py(i,j,k)=C1(i,j,k)*Py1(i,j,k)+C2(i,j,k)*Py2(i,j,k)+C3(i,j,k)*Ex(i,j);
            
            Py2(i,j,k)=Py1(i,j,k);
            Py1(i,j,k)=Py(i,j,k);
            end
            
        end
    end
    
    
    
    Ey(Nx/5,:)=Ey(Nx/5,:)+sin(2*pi*Freq*n*dt);%*exp(-((n*dt-T0)/TT)^2);
    
    %---------------finding Hz
    for  i=1:Nx-1;
        for  j=1:Ny-1;
            Bhzx(i,j)=cbzx1(i,j)*Bhzx(i,j)+cbzx2(i,j)*(Ey(i+1,j)-Ey(i,j));
            Bhzy(i,j)=cbzy1(i,j)*Bhzy(i,j)+cbzy2(i,j)*(Ex(i,j+1)-Ex(i,j));
            Bhz(i,j)=Bhzx(i,j)+Bhzy(i,j);
            
            
        end
    end
    Hz=Bhz./mus;
    
   % Hz(Nx/5,Ny/2)=Hz(Nx/5,Ny/2)+(20/377)*exp(-((n*dt-T0)/TT)^2);
    
    Rec(n)=Ey(Nx/2,jb-2);
    pcolor(log10(abs(Hz')+0.01))
    shading interp
    %mesh(abs(0.5*Ex.*Hz)');
    %axis([1 Nx 1 Ny 0.0 0.001]);
    b=b+1;
    M(b) = getframe;
    
end

%movie(M)
