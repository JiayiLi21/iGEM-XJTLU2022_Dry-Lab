%%%%% Define the x-domain and x-grid %%%%%%%%%%%%%
Lx = 1;%domain:-Lx < x < Lx,unit:m 
Nx = 500;%# of intervals
nx = Nx + 1;%# of gridpoints in x
dx = 2*Lx/Nx;% grid length in x
x = (0:Nx)*dx;%x values on the grid
v= 100;% the speed of the water,unit:m/h
Kd = 6.655;% binding affinity,unit:μM  
n =0.5086;% hill coefficient  
%%% time step parameters %%%%%
nsteps = 2250000;%number of time steps
nouts = 8000;%plot every nout time steps
dt = (dx)^2;%borderline stability of FYCS unit：h
alpha = v*dt/dx;%equation parameter
%%%%% Construct the matrix %%%%%%%%
% Mainly use finite difference method 
% Assume polluted water firstly pass through the device
u  = zeros(Nx+1,1);
u(1,1) = r2;% 0~12 μM random input metal ions  
MBP = 1000000*ones(Nx+1,1);% initial concentration of MBP evenly in the device,unit:μM  
z = zeros(Nx+1,1);% initial concentration of MBP·AgNP complex in each position,unit:μM  
u1= zeros(Nx+1,1);
for j = 1:nsteps
    u(1,1) = r2;
    u1(1,1) = u(1,1);
    % calculate the distribution of the metal ions in each dt
    for i = 2:Nx+1
       
        temp1 = ((1 - alpha)*u(i) + alpha*u(i - 1)); 
        temp =dt*(n*((temp1)^n)*MBP(i))/((Kd + (temp1)^n)*36) ;
        MBP(i) = MBP(i) - temp;
        z(i) = z(i) + temp;
        u1(i) = temp1-temp;
       

        if MBP(i) <=0
            MBP(i) = 0;
        end

        if u1(i) <= 0
            u1(i) = 0;
        end
    end
   u = u1;
   

    %plot
    if mod(j,nouts) == 0 && u(251) <= 0.9
        plot(x(1:251),u(1:251),'linewidth', 1.1); 
        ylim = [0,12];
        set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times');
        title("The distribution of concentration of metal ions in the device with random input concentration of metal ions","FontSize",16);
        xlabel("Distance(m)");
        ylabel("The concentration of the metal ions(μM)");
       
        pause(0.01);
    end

    if u(251) >0.9
       break;
    end
end
   



