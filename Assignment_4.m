%% Defining the problem
length = 0.2;
h = 0.005;
imax = length/h + 1; % no_of_points
%As number of grid points are same in both direction we denote by imax for
%simplicity

x = 0:h:length; %Length of domain
y = 0:h:length; %Height of domain
U = input('Input the velocity of upper lid = ') ; %velocity of lid
nu = 0.004; %Kinematic viscosity

Re = U*length/nu;  %Reynold's number

fprintf('The Reynolds number for given flow condition is = %4.2f \n',Re) 

% Under-relaxation parameter
alpha = 0.8;
%% Initializing the variables

%Final collocated variables
u_final(imax,imax) = 0;
v_final(imax,imax) = 0;
p_final(imax,imax) = 1;
u_final(1,:) = U;

%Staggered variables
u(imax+1,imax) = 0;
u_star(imax+1,imax) = 0;
d_e(imax+1,imax) = 0;
v(imax,imax+1) = 0;
v_star(imax,imax+1) = 0;
d_n(imax,imax+1) = 0;
p(imax+1,imax+1) = 1;
p_star(imax+1,imax+1) = 1;
pc(imax+1,imax+1) = 0;
b(imax+1,imax+1) = 0;
u(1,:) = 2*U;

u_new(imax+1,imax) = 0;
v_new(imax,imax+1) = 0;
p_new(imax+1,imax+1) = 1;
u_new(1,:) = 2*U;

%% Solving using SIMPLE Algorithm

error = 1;
iterations = 0;
convergence = 1e-5; %final required error residual


while error > convergence
    % x-momentum eq. - Interior
    for i = 2:imax
        for j = 2:imax - 1
            u_E = 0.5*(u(i,j) + u(i,j+1));
            u_W = 0.5*(u(i,j) + u(i,j-1));
            v_N = 0.5*(v(i-1,j) + v(i-1,j+1));
            v_S = 0.5*(v(i,j) + v(i,j+1));
            
            a_E = -0.5*u_E*h + nu;
            a_W = 0.5*u_W*h + nu;
            a_N = -0.5*v_N*h + nu;
            a_S = 0.5*v_S*h + nu;
            
            a_e = 0.5*u_E*h - 0.5*u_W*h + 0.5*v_N*h - 0.5*v_S*h + 4*nu;
            
            A_e = -h;
            d_e(i,j) = A_e/a_e;
            
            u_star(i,j) = (a_E*u(i,j+1) + a_W*u(i,j-1) + a_N*u(i-1,j) + a_S*u(i+1,j))/a_e + d_e(i,j)*(p(i,j+1) - p(i,j));
        end
    end
    
    % x-momentum eq. - Boundary
    u_star(1,:) = 2*U - u_star(2,:);
    u_star(imax + 1,:) = -u_star(imax,:);
    u_star(2:imax,1) = 0;
    u_star(2:imax,imax) = 0;
    
    % y-momentum eq. - Interior
    for i = 2:imax - 1
        for j = 2:imax
            u_E = 0.5*(u(i,j) + u(i+1,j));
            u_W = 0.5*(u(i,j-1) + u(i+1,j-1));
            v_N = 0.5*(v(i-1,j) + v(i,j));
            v_S = 0.5*(v(i,j) + v(i+1,j));
            
            a_E = -0.5*u_E*h + nu;
            a_W = 0.5*u_W*h + nu;
            a_N = -0.5*v_N*h + nu;
            a_S = 0.5*v_S*h + nu;
            
            a_n = 0.5*u_E*h - 0.5*u_W*h + 0.5*v_N*h - 0.5*v_S*h + 4*nu;
            
            A_n = -h;
            d_n(i,j) = A_n/a_n;
            
            v_star(i,j) = (a_E*v(i,j+1) + a_W*v(i,j-1) + a_N*v(i-1,j) + a_S*v(i+1,j))/a_n + d_n(i,j)*(p(i,j) - p(i+1,j));
        end
    end
    
    % y-momentum eq. - Boundary
    v_star(:,1) = -v_star(:,2);
    v_star(:,imax + 1) = -v_star(:,imax);
    v_star(1,2:imax) = 0;
    v_star(imax,2:imax) = 0;
    
    % Zeroing the corrections
    pc(1:imax+1,1:imax+1)=0;
    
    % Continuity equation - pressure correction - Interior
    for i = 2:imax
        for j = 2:imax
            a_E = -d_e(i,j)*h;
            a_W = -d_e(i,j-1)*h;
            a_N = -d_n(i-1,j)*h;
            a_S = -d_n(i,j)*h;
            a_P = a_E + a_W + a_N + a_S;
            b(i,j) = -(u_star(i,j) - u_star(i,j-1))*h + (v_star(i,j) - v_star(i-1,j))*h;
            
            pc(i,j) = (a_E*pc(i,j+1) + a_W*pc(i,j-1) + a_N*pc(i-1,j) + a_S*pc(i+1,j) + b(i,j))/a_P;
        end
    end
    
    % Pressure correction
    for i = 2:imax
        for j = 2:imax
            p_new(i,j) = p(i,j) + alpha*pc(i,j);
        end
    end
    
    % Continuity eq. - Boundary
    p_new(1,:) = p_new(2,:);
    p_new(imax + 1,:) = p_new(imax,:);
    p_new(:,1) = p_new(:,2);
    p_new(:,imax + 1) = p_new(:,imax);
    
    % Velocity correction
    for i = 2:imax
        for j = 2:imax - 1
            u_new(i,j) = u_star(i,j) + alpha*d_e(i,j)*(pc(i,j+1) - pc(i,j));
        end
    end
    
    % x-momentum eq. - Boundary
    u_new(1,:) = 2 - u_new(2,:);
    u_new(imax + 1,:) = -u_new(imax,:);
    u_new(2:imax,1) = 0;
    u_new(2:imax,imax) = 0;
    
    for i = 2:imax - 1
        for j = 2:imax
            v_new(i,j) = v_star(i,j) + alpha*d_n(i,j)*(pc(i,j) - pc(i+1,j));
        end
    end
    
    % y-momentum eq. - Boundary
    v_new(:,1) = -v_new(:,2);
    v_new(:,imax + 1) = -v_new(:,imax);
    v_new(1,2:imax) = 0;
    v_new(imax,2:imax) = 0;
            
    
    % Continuity residual as error measure
    error = 0;
    for i = 2:imax
        for j = 2:imax
            error = error + abs(b(i,j));
        end
    end
    
   
    u = u_new;
    v = v_new;
    p = p_new;
    iterations = iterations + 1;
    
end

% After the converged solution, we map the staggered variables to
% collocated variables

for i = 1:imax
    for j = 1:imax
        u_final(i,j) = 0.5*(u(i,j) + u(i+1,j));
        v_final(i,j) = 0.5*(v(i,j) + v(i,j+1));
        p_final(i,j) = 0.25*(p(i,j) + p(i,j+1) + p(i+1,j) + p(i+1,j+1));
    end
end

%%
%Ghia's solution for centreline velocity

dist_y = [0 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 0.5 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1.0];
dist_x = [0 0.0625 0.0703 0.0781 0.0938 0.1563 0.2266 0.2344 0.5 0.8047 0.8594 0.9063 0.9453 0.9531 0.9609 0.9688 1.0];

u_vel_1 = [0 -0.03717 -0.04192 -0.04775 -0.06434 -0.1015 -0.15662 -0.2109 -0.20581 -0.13641 -0.00332 0.23151 0.68717 0.73722 0.78871 0.84123 1];
u_vel_2 = [0 -0.08186 -0.09266 -0.10338 -0.14612 -0.24299 -0.32726 -0.17119 -0.11477 0.02135 0.16256 0.29093 0.55892 0.61756 0.68439 0.75837 1];
v_vel_1 = [0 0.09233 0.10091 0.1089 0.12317 0.16077 0.17507 0.17527 0.05454 -0.24533 -0.22445 -0.16914 -0.10313 -0.08864 -0.07391 -0.05906 0];

%% Centreline u variation - Comparison with benchmark solution
figure(11);
plot(u_final(:,(imax+1)/2),1-y/length, 'LineWidth', 1)

figure(11);
hold on
plot(u_vel_1, dist_y, 'o', 'LineWidth', 1)
xlabel('u')
ylabel('y')
title('u-velocity Centerline variation (Re = 400)')
legend('SIMPLE solution', 'Ghia solution', 'location', 'southeast')

%% Centreline v variation - Comparison with benchmark solution
figure(12);
plot(v_final((imax+1)/2,:),x/length, 'LineWidth', 1)

figure(12);
hold on
plot(v_vel_1, dist_x, 'o', 'LineWidth', 1)
xlabel('v')
ylabel('y')
title('v-velocity Centerline variation (Re = 400)')
legend('SIMPLE solution', 'Ghia solution', 'location', 'southeast')


%% Contour and vector visuals.
x_dom = ((1:imax)-1).*h;
y_dom = 1-((1:imax)-1).*h;
[X,Y] = meshgrid(x_dom,y_dom);
figure(21);
contourf(X,Y,u_final, 21, 'LineStyle', 'none')
colorbar
colormap('jet')
title('Velocity Contour Plot (U=8)')
xlabel('x')
ylabel('y')


figure(22);
quiver(X, Y, u_final, v_final, 5, 'k')
title('Velocity Vector plot for U = 8 m/s')
xlim([0 0.2])

%%
u_avg = mean(u_final(2:(imax-1),(imax+1)/2));
v_avg = mean(v_final((imax+1)/2, 2:(imax-1)));

%% %Richardson extrapolation scheme

p = log((-0.0355+0.0218)/(-0.0218+0.0118))/log(2);

error_1 = (-0.0355+0.0218)/(1 - 2^p);
error_2 = (-0.0218+0.0118)/(1 - 2^p);
