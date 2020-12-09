clr
%% Input parameter values
Rabx=0:0.1:1; %dimensionless paramater, Rab = qin/(beta*nu)

for ii=1:length(Rabx),
    Rab=Rabx(ii);
m=0; %set m = 1 to display movie
%{
Qriver = 100; %m3/yr;
nu = 1e5; %m2/yr
w = 10; %m, delta width

beta = 1e-4; %basement slope

Rab = Qriver./(w.*beta.*nu);
%}
%% setting up the domain and time
nx=1000;
p=500; %initial shoreline position
pos=p;
dx=0.01; %m

dt=0.00005; %yr
Tmax=1;
t=0:dt:Tmax;
nt=length(t);

%courant check:

%% Initialize time vectors
s=zeros(1,nt); % Shoreline position
r=zeros(1,nt); % Alluvial-basement position
Z = zeros(1,nt);   %sea level through time
%Z = sin(pi*t./4); %(m)
%% Initialize space vectors

%cycle values
x = zeros(1,nx);%basement elevation
x=dx.*((1:nx)-p+0.5);
E = -x;
h = max(E,0);%height above current sea-level
%E = max(E,-0.000001);

F=zeros(1,nx); % Sediment flux (m^2/y)
%Fx = zeros(1,nx);
H=zeros(1,nx); %Enthalpy
Hnew=zeros(1,nx);
di = zeros(nt,nx); % Elevation stored for frames of movie

% To define the dimension of F(flux)
F(1) = Rab; %Upstream boundary condition

for j=1:nt
    
    %s(j)=(pos-p)*dx - H(pos)*dx/(E(pos)-Z(j)); %shoreline position
    rcount = 0; % used to track ABT later
    
 %   if t(j)>1,
 %       disp('blub')
 %   end
    for i=2:nx-1
        F(i) = min((h(i) - h(i+1))/dx, H(i)*dx/dt+F(i-1)); %Flux Condition
        
        % Locate first cell where sediment deposition occurs
   %     if rcount == 0 && F(i) ~= F(i-1)
    %        rcount = 1;
    %        r(j)=(h(i) + Z(j) + Rab*x(i))/(1-Rab); %Interpolated ABT position
    %    end
        
        % Exner equation
        Hnew(i) = H(i) + dt/dx*(F(i-1)-F(i));
        
        % tracks SH regression
     %   if Hnew(pos) + E(pos)> Z(j)
     %       pos = pos + 1;
      %  end
        
        % tracks SH transgression
    %    while H(pos-1) + E(pos-1) < Z(j)
     %       pos = pos-1;
    %    end
 
        H(i)=Hnew(i); %seed new enthalpy values
        h(i) = max(H(i)+E(i)-Z(j),0); %determine height above current SL
        
        
    end
        
    %% Movie
    if m==1 && mod(j,.05/dt) == 0
        di(j,:) = H+E;
    end
    
    %% Counter
    if mod(j,.5/dt) == 0
        disp(['t=', num2str(j*dt)]);
    end
    
    %% Checks if boundaries of domain are hit
    if j>1000 %skips roughness at first few steps
        if pos == nx-1 || r(j)/dx >= p || abs(r(j) - r(j-1)) >= 100*dx
            disp('end of domain reached')
            break
        end
    end    
end

Fx(ii) = F(pos-1); %.*beta.*nu.*w; %m3/yr
ret(ii) = (F(1)-F(pos-1))./F(1);
end
plot(ret), hold on

%movie
if m ==1
    run('Enthalpy_Code_Movie')
end
