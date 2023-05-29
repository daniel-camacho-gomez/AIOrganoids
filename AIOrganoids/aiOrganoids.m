function [n_c_evo] = aiOrganoids(t_t, n_ct, h_1, h_2, o_1, o_2)

tic

addpath(genpath('nnFunctions'))

%-------------------------------------------------------------------------
%~~~~~~~~~~~~~~~~~~~~~~~NUMERICAL PARAMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%-------------------------------------------------------------------------

%...Initial number of cells
n_c         = 1;
%...Counter for the target values
evo_t       = 1; %counter


%...Time steps
    t_sim     = t_t(end)*24*60;   %min  total time of simulation
    t_total   = 0;                %time to track elapsed time 
    t_plot    = 0;                %min
    t_shape   = 0;
    h_counter = 0; 
    
    t_mech    = 0.1;        %min   to track elapsed cell mechanics
    t_cell    = 6;          %min   to track elapsed cell cycle 

    Dt_mech   = 0.1;        %min   time step for mechanics
    Dt_cell   = 6;          %min   time step for cell cycle
    Dt_plot   = 2*15;       %min   
    Dt_shape  = 1*60;       %min

%...Domain limits 
    x_dom = [-100 100];    %um 
    y_dom = [-100 100];    %um 
    z_dom = [-100 100];    %um

   
%-------------------------------------------------------------------------
%~~~~~~~~~~~~~~~~~~~~~~~~~~~VARIABLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%-------------------------------------------------------------------------
   
   
%~~~~~~~~~~~~~~~~~CELLS~~~~~~~~~~~~~~~~~~~~~
%...At time n
    x_c = zeros(n_c,1);      %cell position
    y_c = zeros(n_c,1);      %cell position
    z_c = zeros(n_c,1);      %cell position
    
    v_x = zeros(n_c,1);      %cell velocity
    v_y = zeros(n_c,1);      %cell velocity
    v_z = zeros(n_c,1);      %cell velocity
    
    R_c = zeros(n_c,1);      %cell Radius     
    
%...At time n+1 
    x_c_n = zeros(n_c,1);    %cell position 
    y_c_n = zeros(n_c,1);    %cell position
    z_c_n = zeros(n_c,1);    %cell position
    
    v_x_n = zeros(n_c,1);    %cell velocity
    v_y_n = zeros(n_c,1);    %cell velocity
    v_z_n = zeros(n_c,1);    %cell velocity

  
%...At time n-1 
    v_x_nm = zeros(n_c,1);   %cell velocity
    v_y_nm = zeros(n_c,1);   %cell velocity
    v_z_nm = zeros(n_c,1);   %cell velocity
    
%...Forces 
    F_c_x    = zeros(n_c,1);
    F_c_y    = zeros(n_c,1);
    F_c_z    = zeros(n_c,1);    
    
%...cell cycle        
    state       = zeros(n_c,1);  %phase state of each cell

    
%-------------------------------------------------------------------------
%~~~~~~~~~~~~~~~~~~~~~~~~Initialization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%-------------------------------------------------------------------------

%...Initial position of the cells 
    %...cell 1
    x_c(1)   = 0;  
    y_c(1)   = 0;
    z_c(1)   = 0;

%-------------------------------------------------------------------------
%~~~~~~~~~~~~~~~~~~~~~~~~~~~PARAMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%------------------------------------------------------------------------- 

%-------------------------------
%----------CELL CYCLE----------- 
%-------------------------------
T_cycle     = 24*60;         % min
alpha       = 1/T_cycle;     % 1/min rate of nuclear biomass creation

%...initial state 
state(:)    = 2;             % 1=inactive  2=active

%-------------------------------
%----------CELL VOLUME---------- 
%-------------------------------
R_c(:) = 10;                   % um initial cells radii
V      = 4/3*pi*R_c^3;         % R_i = 10 um m
V_t    = 2*V;                  % target volume

%-------------------------------
%----------CELL MECHANICS------- 
%-------------------------------
F_comp_cc     = -300*0.016*1e3*3600;      %pN --> ug/um min^2 (1e3*3600)
F_adh_cc      =  1500*0.016*1e3*3600;     %PN --> ug/um min^2 (1e3*3600)

%potentials parameters 
lambda = 7;
x_o    = sqrt(1/(2*lambda)); 
v_o    = x_o*exp(-lambda*x_o^2); 

%-------------------------------
%-----------Matrix-------------- 
%-------------------------------
eta      = 20*1e3*60;                     % Pa s --> ug/um min (1e3*60)


%alpha-shapes
theta = 0:0.5:2*pi; 
gamma = 0:0.5:2*pi;

%-------------------------------------------------------------------------
%~~~~~~~~~~~~~~~~~~~~~~~~~SIMULATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%------------------------------------------------------------------------- 
counter        = 1;  
counter_steps  = 1;

while t_total < t_sim
    
    
%-------------------------------------------------------------------------    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~CELL CYCLE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%-------------------------------------------------------------------------
    
if t_total >= t_cell

  n_c_nm = n_c;   

   for s = 1:n_c_nm 

     switch state(s)
           
         case 1 %Inactivity/Quiescent
         
            %INPUTS
            state(s) = 0;
                 
            %state of cells    
            nci = sum(state(:)==1)/n_c;  %number of cells inactive
            nca = sum(state(:)==2)/n_c;  %number of cells activite

            I   = [(n_ct(evo_t) - n_c)/n_ct(evo_t), nci, nca]; 
            
            decision = nnForwardPropagation(I, h_1, h_2, o_1, o_2);
            
            [max_d, id_d] = max(decision);  
            
            state(s)  = id_d; 
            
            
            if id_d == 2         %if decision = active
               V_t(s) = 2*V(s);  %duplicate target volume 
            end

            
         case 2 %Active
             
%-------------------------------------------------------------------------   
%                             Growth
%-------------------------------------------------------------------------             
         if V(s) >=  V_t(s)       %Check if it's time to divide
             
%-------------------------------------------------------------------------   
%                             Division
%-------------------------------------------------------------------------               
                       
i = s;               
         
        %...create a new cell
        n_c = n_c + 1; 

        %...cell cycle 
             %...parent cell  
             
             %...daughter cell
                
        %...cell volume 
             %...daughter cell 2
                V(n_c)          = 0.5*V(i);    % total cell volume 
             %...daughter cell 1 
                V(i)            = 0.5*V(i);    % total cell volume    

                V_t(n_c)        = V(n_c);
                V_t(i)          = V(i);
                
            
          %Forces daughter cell
                F_c_x(n_c)   = F_c_x(i);
                F_c_y(n_c)   = F_c_y(i);
                F_c_z(n_c)   = F_c_z(i);                

                
           %random division
%                 u = [rand rand rand]; 
                  u = [-1+2*rand(1,1) -1+2*rand(1,1) -1+2*rand(1,1)]; 

                  orth = null(u);  
    
       x_cent = x_c(i); 
       y_cent = y_c(i);      
       z_cent = z_c(i); 

     
    %...daughter position and velocity
    x_c(n_c) = x_c(i) + ((R_c(i) - (1/2^(1/3))*R_c(i)))*orth(1,1);    %cell position
    y_c(n_c) = y_c(i) + ((R_c(i) - (1/2^(1/3))*R_c(i)))*orth(2,1);    %cell position 
    z_c(n_c) = z_c(i) + ((R_c(i) - (1/2^(1/3))*R_c(i)))*orth(3,1);    %cell position
    
    x_c(i)   = x_c(i) - ((R_c(i) - (1/2^(1/3))*R_c(i)))*orth(1,1);    %cell position
    y_c(i)   = y_c(i) - ((R_c(i) - (1/2^(1/3))*R_c(i)))*orth(2,1);    %cell position 
    z_c(i)   = z_c(i) - ((R_c(i) - (1/2^(1/3))*R_c(i)))*orth(3,1);    %cell position      %cell position
      
    
    rot = rand(1,1)*pi;       


x_c1 = x_c(i);
y_c1 = y_c(i);
z_c1 = z_c(i);

x_c2 = x_c(n_c);
y_c2 = y_c(n_c);
z_c2 = z_c(n_c);



 x_c(i) = x_cent  + (x_c1-x_cent)*(cos(rot) + u(1)^2*(1-cos(rot)))          ...
                  + (y_c1-y_cent)*(u(1)*u(2)*(1-cos(rot)) - u(3)*sin(rot))   ...
                  + (z_c1-z_cent)*(u(1)*u(3)*(1-cos(rot)) + u(2)*sin(rot)); 
        
 y_c(i) = y_cent  + (x_c1-x_cent)*(u(2)*u(1)*(1-cos(rot)) + u(3)*sin(rot))  ... 
                  + (y_c1-y_cent)*(cos(rot) + u(2)^2*(1-cos(rot)))        ...
                  + (z_c1-z_cent)*(u(2)*u(3)*(1-cos(rot)) - u(1)*sin(rot)); 
        
 z_c(i) = z_cent  + (x_c1-x_cent)*(u(3)*u(1)*(1-cos(rot)) - u(2)*sin(rot)) ...
                  + (y_c1-y_cent)*(u(3)*u(2)*(1-cos(rot)) + u(1)*sin(rot)) ...
                  + (z_c1-z_cent)*(cos(rot) + u(3)^2*(1-cos(rot)));
               
               
 x_c(n_c) = x_cent + (x_c2-x_cent)*(cos(rot) + u(1)^2*(1-cos(rot)))          ...
                   + (y_c2-y_cent)*(u(1)*u(2)*(1-cos(rot)) - u(3)*sin(rot))   ...
                   + (z_c2-z_cent)*(u(1)*u(3)*(1-cos(rot)) + u(2)*sin(rot)); 
        
 y_c(n_c) = y_cent  + (x_c2-x_cent)*(u(2)*u(1)*(1-cos(rot)) + u(3)*sin(rot))  ... 
                    + (y_c2-y_cent)*(cos(rot) + u(2)^2*(1-cos(rot)))        ...
                    + (z_c2-z_cent)*(u(2)*u(3)*(1-cos(rot)) - u(1)*sin(rot)); 
        
 z_c(n_c) = z_cent  + (x_c2-x_cent)*(u(3)*u(1)*(1-cos(rot)) - u(2)*sin(rot)) ...
                    + (y_c2-y_cent)*(u(3)*u(2)*(1-cos(rot)) + u(1)*sin(rot)) ...
                    + (z_c2-z_cent)*(cos(rot) + u(3)^2*(1-cos(rot)));
                            

                %...At time n     
                    v_x(n_c)    = v_x(i);      %cell velocity
                    v_y(n_c)    = v_y(i);      %cell velocity
                    v_z(n_c)    = v_z(i);      %cell velocity

                %...At time n-1 
                    v_x_nm(n_c) = v_x_nm(i);   %cell velocity 
                    v_y_nm(n_c) = v_y_nm(i);   %cell velocity 
                    v_z_nm(n_c) = v_z_nm(i);   %cell velocity               
                
                    
                    R_c(i)      = (3*V(i)/(4*pi)).^(1/3);      
                    R_c(n_c)    = (3*V(n_c)/(4*pi)).^(1/3);    


                    
         % After division 
         % NN call to decide for each cell 
         
         %i-cell
            %INPUTS
            state(i) = 0;
            
            
            %state of cells    
            nci = sum(state(:)==1)/n_c;  %number of cells inactive
            nca = sum(state(:)==2)/n_c;  %number of cells activite

            I   = [(n_ct(evo_t) - n_c)/n_ct(evo_t), nci, nca]; 
            
            decision = nnForwardPropagation(I, h_1, h_2, o_1, o_2);
            
            [max_d, id_d] = max(decision);  
            
            state(i)  = id_d; 
            
            
            if id_d == 2         %if decision = active
               V_t(i) = 2*V(i);  %duplicate target volume 
            end

            
         %n_c-cell
         
            %INPUTS
            state(n_c) = 0;
            

            %state of cells    
            nci = sum(state(:)==1)/n_c;  %number of cells inactive
            nca = sum(state(:)==2)/n_c;  %number of cells activite


            I   = [(n_ct(evo_t) - n_c)/n_ct(evo_t), nci, nca]; 
            
            decision = nnForwardPropagation(I, h_1, h_2, o_1, o_2);
            
            [max_d, id_d] = max(decision);  
            
            state(n_c)  = id_d; 
            
            
            if id_d == 2         %if decision = active
               V_t(n_c) = 2*V(n_c);  %duplicate target volume 
            end

            
            
         else %If it hasn't duplicate its volume, keep growing
            
%-------------------------------------------------------------------------   
%                               GROWTH
%-------------------------------------------------------------------------               
                 %Volume ODE
                 V(s)   =  V(s) + Dt_cell*alpha*V(s);


         end  %division/proliferation, activity
         
         
                 
     end %switch 
   end %for every cell

    
   
   t_cell = t_cell + Dt_cell;        

end %cell cycle
    
%...Cell radius 
R_c = (3/4*V/pi).^(1/3); 
      
    

%-------------------------------------------------------------------------    
%~~~~~~~~~~~~~~~~~~~~~~POSITION CALCULATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%-------------------------------------------------------------------------
   %...Position calculation
   for i=1:n_c
        x_c_n(i) = x_c(i) + 0.5*Dt_mech*(3*v_x(i) - v_x_nm(i)); 
        y_c_n(i) = y_c(i) + 0.5*Dt_mech*(3*v_y(i) - v_y_nm(i)); 
        z_c_n(i) = z_c(i) + 0.5*Dt_mech*(3*v_z(i) - v_z_nm(i));
   end

   %...Velocity reallocation for next time step
   v_x_nm = v_x; 
   v_y_nm = v_y;
   v_z_nm = v_z;


%-------------------------------------------------------------------------   
%~~~~~~~~~~~~~~~~~~~~~INTERACTING FORCES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%-------------------------------------------------------------------------

%----------------- 
%---CELL FORCES---
%-----------------

for i=1:n_c  
    
    
    %...Forces initialization for i-cell
        F_cc_x  = 0;
        F_cc_y  = 0;
        F_cc_z  = 0;   
  
    
    min_dist = -0.1*R_c(i);         
        
    %...CELL - CELL INTERACTION            
    
       for j=1:n_c
           if j ~= i  
                 
                 chi = (R_c(i)/2)*(1/R_c(i) + 1/R_c(j));               
               
                 %...distance between cell i and j 
                 r    = - [x_c(i)-x_c(j);y_c(i)-y_c(j);z_c(i)-z_c(j)];
                 r_m  = sqrt(r(1)^2+r(2)^2+r(3)^2); 

                 d    = r_m - R_c(i) - R_c(j); 

                 xx   = (d - min_dist)/R_c(i);
                 
                 
                    if xx < 0
                        F_cc =   F_comp_cc*chi*(-xx)^(3/2);
                    else 
                        F_cc = - F_adh_cc*chi*((xx+x_o)*exp(-lambda*(xx+x_o)^2)...
                               - v_o*exp(-lambda*xx^2));        
                    end  
                                     
                F_cc_x = F_cc_x + F_cc*r(1)/r_m; 
                F_cc_y = F_cc_y + F_cc*r(2)/r_m; 
                F_cc_z = F_cc_z + F_cc*r(3)/r_m;                    

           end
       end
    
               %CELL FORCES   
               F_c_x(i) =  F_cc_x; 
               F_c_y(i) =  F_cc_y;  
               F_c_z(i) =  F_cc_z;     
              
end %end interacting forces




%-------------------------------------------------------------------------    
%~~~~~~~~~~~~~~~~~~~~~~~~~VELOCITY UPDATE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%-------------------------------------------------------------------------

%...CELLS 
 for i=1:n_c
    v_x(i) =  (1/(6*pi*eta*R_c(i)))*F_c_x(i);
    v_y(i) =  (1/(6*pi*eta*R_c(i)))*F_c_y(i);
    v_z(i) =  (1/(6*pi*eta*R_c(i)))*F_c_z(i);
 end 

 
 
%-------------------------------------------------------------------------    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~ANIMATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%-------------------------------------------------------------------------

t_total_h   = fix(t_total/60);

if t_total_h >= h_counter 
    
    elap_h =  fix(toc/3600);
    elap_m =  fix(toc/60)-60*fix(toc/3600);
    elap_s =  toc-60*elap_m-3600*elap_h;
    
    disp("--------------------------------------")
    disp("   Elapsed time: " + elap_h + ' h ' +elap_m+ ' min ' +elap_s+ ' s');
    disp("   Simulated time:  "+ h_counter/24 + ' days'); 
    disp("   Number of cells:  " + n_c); 
    disp("--------------------------------------")
    h_counter = h_counter+24; 
end



%---------------------------------ANIMATION--------------------------------

% if t_total >= t_plot 
%                      
%   t_plot = t_plot + Dt_plot; 
%  
%   t_total_h   = fix(t_total/60);
%   t_total_min = t_total-60*t_total_h;
% 
%   
% 
% %...ALPHA-SHAPE
% theta = 0:0.5:2*pi; 
% gamma = 0:0.5:2*pi;
% 
% clear x_ca 
% clear y_ca 
% clear z_ca 
% 
% counter_p = 1; 
% for c=1:n_c 
% %      if z_c(c)<= 0 
%     for i=1:length(theta) 
%         for j=1:length(gamma) 
%                 
%                 x_ca(counter_p) = x_c(c) + R_c(c)*sin(theta(i))*cos(gamma(j));
%                 y_ca(counter_p) = y_c(c) + R_c(c)*sin(theta(i))*sin(gamma(j));
%                 z_ca(counter_p) = z_c(c) + R_c(c)*cos(theta(i));
%                 
%                 counter_p = counter_p+1; 
%         end  
%     end
% %      end
% end
% 
% 
% 
% shp1 = alphaShape(x_ca',y_ca',z_ca',20); 
% h1=plot(shp1);
% 
% set(h1,'FaceColor',[44/255 164/255 0/255], ...
%   'FaceAlpha',0.7,'EdgeColor','none','FaceLighting','gouraud')
% hold on
% 
% 
% 
% 
% 
% %...CELL CENTERS
% 
% for i = 1:n_c 
% 
%         %     if z_c(i)<= 0
%         %Real magnitude    
%         [xs,ys,zs] = sphere;     
%         xs = xs*R_c(i);
%         ys = ys*R_c(i); 
%         zs = zs*R_c(i); 
% 
%         hSurface=surf(xs+x_c(i),ys+y_c(i),zs+z_c(i));
%           set(hSurface,'FaceColor',[255/255 71/255 5/255], ...
%           'FaceAlpha',1,'FaceLighting','gouraud','EdgeColor','none')
%         daspect([1 1 1]);
%         camlight
%         hold on
%         material metal
%         %    end
% end
% 
% 
% 
% 
% 
% % set(gca,'facelighting','none');
% 
%  z_scale_bar  = [-100 -100]; 
%  y_scale_bar  = [-100 -100];
%  x_scale_bar  = [-100 -100]; 
%  
% %  plot3(x_scale_bar, y_scale_bar, z_scale_bar,'Color','w','LineWidth',5)
% %  text(-97,-93,0, '30 $\mu m$','Interpreter','Latex','FontSize',10,'Color','w')
%  
%  
% %  z_scale_bar  = [0 0]; 
% %  y_scale_bar  = [-100 -70];
% %  x_scale_bar  = [100 100]; 
% %  
% % plot3(x_scale_bar, y_scale_bar, z_scale_bar,'Color','w','LineWidth',5)
% %  text(93.,-97,0, '30 $\mu m$','Interpreter','Latex','FontSize',10,'Color','w')
% %  
%  
% 
% 
% % xlabel('($\mu m$)','Interpreter','Latex','FontSize',15)
% % ylabel('($\mu m$)','Interpreter','Latex','FontSize',15)
% % zlabel('($\mu m$)','Interpreter','Latex','FontSize',15)
% 
% 
% 
% %  set(gca,'Color','k')
%  set(gca,'Color',[38/255 38/255 38/255])
% %  yticks([])
% %  xticks([])
% %  zticks([])
% 
% 
% 
% xlim([x_dom(1) x_dom(2)])
% ylim([y_dom(1) y_dom(2)]) 
% zlim([z_dom(1) z_dom(2)])
% xlabel('$X\;(\mu m)$','Interpreter','Latex','FontSize',15)
% ylabel('$Y\;(\mu m)$','Interpreter','Latex','FontSize',15)
% zlabel('$Z\;(\mu m)$','Interpreter','Latex','FontSize',15)
% 
% % title("Time =  " + t_total_h + " h " + t_total_min + " min")
% title("Time =  " + t_total_h + " h ")
% set(gca,'FontSize',10,'FontName','Times New Roman')
% 
% Wfinal(counter) = getframe(gcf);
% counter=counter+1; 
% pause(0.1)
% hold off
% 
% end


%--------------------------------------------------------------------------




% if t_total >= t_plot 
%                      
%   t_plot = t_plot + Dt_plot; 
%  
%   t_total_h   = fix(t_total/60);
%   t_total_min = t_total-60*t_total_h;
% 
%   
%   
%   
%     for i = 1:n_c 
%         
%      %Real magnitude    
%        [xs,ys,zs] = sphere;     
%        xs = xs*R_c(i);
%        ys = ys*R_c(i); 
%        zs = zs*R_c(i); 
%     
%         if state(i) == 1 % inactivity 
%             
%             hSurface = surf(xs+x_c(i),ys+y_c(i),zs+z_c(i));
%                   set(hSurface,'FaceColor',[0 1 0], ...
%                   'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none')
%                   daspect([1 1 1]);
%             camlight
%             hold on 
%             
%         else  %proliferation
% 
%             hSurface = surf(xs+x_c(i),ys+y_c(i),zs+z_c(i));
%                   set(hSurface,'FaceColor',[0 0 1], ...
%                   'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none')
%                   daspect([1 1 1]);
%             camlight
%             hold on 
%             
%         end
%     end
% 
% 
% xlim([x_dom(1) x_dom(2)])
% ylim([y_dom(1) y_dom(2)]) 
% zlim([z_dom(1) z_dom(2)])
% 
% xlabel('$X\;(\mu m)$','Interpreter','Latex','FontSize',15)
% ylabel('$Y\;(\mu m)$','Interpreter','Latex','FontSize',15)
% zlabel('$Z\;(\mu m)$','Interpreter','Latex','FontSize',15)
% 
% title("Time =  " + t_total_h + " h " + t_total_min + " min")
% set(gca,'FontSize',10,'FontName','Times New Roman')
% 
% 
% Wfinal(counter) = getframe(gcf);
% counter=counter+1; 
% pause(0.1)
% hold off
% 
% end



%-------------------------------------------------------------------------   
%~~~~~~~~~~~~~~~~~~~~~~~~~REALLOCATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%-------------------------------------------------------------------------

%...Reallocation of variables for next time step
    x_c = x_c_n; 
    y_c = y_c_n; 
    z_c = z_c_n; 
    
    
%Update the target values
t_total = t_total + Dt_mech;   


    if  t_total >= t_t(evo_t)*24*60
        
        n_c_evo(evo_t)   = n_c; 
        t_n_c_evo(evo_t) = (t_total)/(24*60);

     
% %----------------------------SAVE FIGURE-------------------------
%         
% theta = 0:0.5:2*pi; 
% gamma = 0:0.5:2*pi;
% 
% clear x_ca 
% clear y_ca 
% clear z_ca 
% 
% counter_p = 1; 
% for c=1:n_c 
% %      if z_c(c)<= 0 
%     for i=1:length(theta) 
%         for j=1:length(gamma) 
%                 
%                 x_ca(counter_p) = x_c(c) + R_c(c)*sin(theta(i))*cos(gamma(j));
%                 y_ca(counter_p) = y_c(c) + R_c(c)*sin(theta(i))*sin(gamma(j));
%                 z_ca(counter_p) = z_c(c) + R_c(c)*cos(theta(i));
%                 
%                 counter_p = counter_p+1; 
%         end  
%     end
% %      end
% end
% 
% 
% 
% shp1 = alphaShape(x_ca',y_ca',z_ca',20); 
% h1=plot(shp1);
% 
% set(h1,'FaceColor',[44/255 164/255 0/255], ...
%   'FaceAlpha',0.7,'EdgeColor','none','FaceLighting','gouraud')
% hold on
% 
% 
% 
% 
% %...CELL CENTERS
% 
% for i = 1:n_c 
%         %     if z_c(i)<= 0
%             %Real magnitude    
%             [xs,ys,zs] = sphere;     
%             xs = xs*R_c(i);
%             ys = ys*R_c(i); 
%             zs = zs*R_c(i); 
% 
%             hSurface=surf(xs+x_c(i),ys+y_c(i),zs+z_c(i));
%               set(hSurface,'FaceColor',[255/255 71/255 5/255], ...
%               'FaceAlpha',1,'FaceLighting','gouraud','EdgeColor','none')
%             daspect([1 1 1]);
%             camlight
%             hold on
%             material metal
%         %    end
% end
% 
% 
% 
% 
% % set(gca,'facelighting','none');
% 
%  z_scale_bar  = [-100 -100]; 
%  y_scale_bar  = [-100 -100];
%  x_scale_bar  = [-100 -100]; 
%  
%  
% set(gca,'Color',[38/255 38/255 38/255])
% 
% 
% xlim([x_dom(1) x_dom(2)])
% ylim([y_dom(1) y_dom(2)]) 
% zlim([z_dom(1) z_dom(2)])
% xlabel('$X\;(\mu m)$','Interpreter','Latex','FontSize',15)
% ylabel('$Y\;(\mu m)$','Interpreter','Latex','FontSize',15)
% zlabel('$Z\;(\mu m)$','Interpreter','Latex','FontSize',15)
% 
% set(gca,'FontSize',10,'FontName','Times New Roman')        
%         
% 
%  fname = sprintf('Solid_org_evo_c3_t%d.fig', t_t(evo_t));
%  saveas(gca,fname);
%  hold off
%  
%  
%  
%  
% %------------------------------------------------------------------------        
        

        evo_t            = evo_t + 1;        
    end
    
    
%...store temporal data     
n_c_time(counter_steps)    = n_c; 
for i=1:n_c
    state_evo(i,counter_steps) = state(i);
end


counter_steps = counter_steps + 1; 
end %of simulation; while



%-------------------------------------------------------------------------
%~~~~~~~~~~~~~~~~~~~~~~~END OF SIMULATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%-------------------------------------------------------------------------

t_total_h   = fix(t_total/60);

if t_total_h >= h_counter 
    
    elap_h =  fix(toc/3600);
    elap_m =  fix(toc/60)-60*fix(toc/3600);
    elap_s =  toc-60*elap_m-3600*elap_h;
    
    disp("--------------------------------------")
    disp("   Elapsed time: " + elap_h + ' h ' +elap_m+ ' min ' +elap_s+ ' s');
    disp("   Simulated time:  "+ h_counter/24 + ' days'); 
    disp("   Number of cells:  " + n_c); 
    disp("--------------------------------------")
    h_counter = h_counter+24; 
end



%...Save the simulation 
savedir = fullfile(pwd,'/Results');
fname = sprintf('AISolidOrganoid.mat');
save(fullfile(savedir,fname));

end




%% REPRESENTATION 

%     
%   t_plot = t_plot + Dt_plot; 
%  
%   t_total_h   = fix(t_total/60);
%   t_total_min = t_total-60*t_total_h;
% 
%   figure
%     for i = 1:n_c 
%         
%      %Real magnitude    
%        [xs,ys,zs] = sphere;     
%        xs = xs*R_c(i);
%        ys = ys*R_c(i); 
%        zs = zs*R_c(i); 
%     
%         if state(i) == 1 % inactivity 
%             
%             hSurface = surf(xs+x_c(i),ys+y_c(i),zs+z_c(i));
%                   set(hSurface,'FaceColor',[0 1 0], ...
%                   'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none')
%                   daspect([1 1 1]);
%             camlight
%             hold on 
%             
%         else  if state(i) == 2 %proliferation
% 
%             hSurface = surf(xs+x_c(i),ys+y_c(i),zs+z_c(i));
%                   set(hSurface,'FaceColor',[0 0 1], ...
%                   'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none')
%                   daspect([1 1 1]);
%             camlight
%             hold on 
%             
%         else 
% 
%             hSurface = surf(xs+x_c(i),ys+y_c(i),zs+z_c(i));
%                     set(hSurface,'FaceColor',[1 0 0], ...
%                     'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none')
%                     daspect([1 1 1]);
%             camlight
%             hold on   
%             
%             end
%         end
%     end
% 
% 
% xlim([x_dom(1) x_dom(2)])
% ylim([y_dom(1) y_dom(2)]) 
% zlim([z_dom(1) z_dom(2)])
% 
% xlabel('$X\;(\mu m)$','Interpreter','Latex','FontSize',15)
% ylabel('$Y\;(\mu m)$','Interpreter','Latex','FontSize',15)
% zlabel('$Z\;(\mu m)$','Interpreter','Latex','FontSize',15)
% 
% title("Time =  " + t_total_h + " h " + t_total_min + " min")
% set(gca,'FontSize',10,'FontName','Times New Roman')




%% ALPHA-SHAPES
% 
%             %alpha-shapes for volume estimation
%             clear x_ca 
%             clear y_ca 
%             clear z_ca 
%             counter_p = 1; 
%             
%             for c=1:n_c 
%                 
%                 for k=1:length(theta) 
%                     for l=1:length(gamma) 
% 
%                     x_ca(counter_p) = x_c(c) + R_c(c)*sin(theta(k))*cos(gamma(l));
%                     y_ca(counter_p) = y_c(c) + R_c(c)*sin(theta(k))*sin(gamma(l));
%                     z_ca(counter_p) = z_c(c) + R_c(c)*cos(theta(k));
% 
%                     counter_p = counter_p+1; 
%                 
%                     end  
%                 end
% 
%             end
% 
%             figure
%             shp1 = alphaShape(x_ca',y_ca',z_ca');   
%             v_l  = volume(shp1);
%             
%             
% h1 = plot(shp1);
% 
%  set(h1,'FaceColor',[0 0 1], ...
%   'FaceAlpha',0.7,'EdgeColor','none','FaceLighting','none')
% 
% xlim([x_dom(1) x_dom(2)])
% ylim([y_dom(1) y_dom(2)]) 
% zlim([z_dom(1) z_dom(2)])
% 
% xlabel('$X\;(\mu m)$','Interpreter','Latex','FontSize',15)
% ylabel('$Y\;(\mu m)$','Interpreter','Latex','FontSize',15)
% zlabel('$Z\;(\mu m)$','Interpreter','Latex','FontSize',15)
% 
% title("Time =  " + t_total_h + " h " + t_total_min + " min")
% set(gca,'FontSize',10,'FontName','Times New Roman')

% surfaceArea(shp1)


%% ANIMATION 

% % % create the video writer with 1 fps
%   writerObj = VideoWriter('Solid_org_evo_c3.avi');
%   writerObj.FrameRate = 15;
%   
% % set the seconds per image
% % open the video writer
%   open(writerObj);
%   
% % write the frames to the video
% 
% for i=1:length(Wfinal)
%     % convert the image to a frame
%      frame = Wfinal(i) ;    
%      writeVideo(writerObj, frame);
% end
% 
% % close the writer object
% close(writerObj);


%% Postprocesado 

% %...ALPHA-SHAPE
% theta = 0:0.5:2*pi; 
% gamma = 0:0.5:2*pi;
% 
% clear x_ca 
% clear y_ca 
% clear z_ca 
% 
% counter_p = 1; 
% for c=1:n_c 
% %      if z_c(c)<= 0 
%     for i=1:length(theta) 
%         for j=1:length(gamma) 
%                 
%                 x_ca(counter_p) = x_c(c) + R_c(c)*sin(theta(i))*cos(gamma(j));
%                 y_ca(counter_p) = y_c(c) + R_c(c)*sin(theta(i))*sin(gamma(j));
%                 z_ca(counter_p) = z_c(c) + R_c(c)*cos(theta(i));
%                 
%                 counter_p = counter_p+1; 
%         end  
%     end
% %      end
% end
% 
% 
% 
% shp1 = alphaShape(x_ca',y_ca',z_ca',20); 
% h1=plot(shp1);
% 
% set(h1,'FaceColor',[44/255 164/255 0/255], ...
%   'FaceAlpha',0.7,'EdgeColor','none','FaceLighting','gouraud')
% hold on
% 
% 
% 
% 
% 
% %...CELL CENTERS
% 
% for i = 1:n_c 
% 
%         %     if z_c(i)<= 0
%         %Real magnitude    
%         [xs,ys,zs] = sphere;     
%         xs = xs*R_c(i);
%         ys = ys*R_c(i); 
%         zs = zs*R_c(i); 
% 
%         hSurface=surf(xs+x_c(i),ys+y_c(i),zs+z_c(i));
%           set(hSurface,'FaceColor',[255/255 71/255 5/255], ...
%           'FaceAlpha',1,'FaceLighting','gouraud','EdgeColor','none')
%         daspect([1 1 1]);
%         camlight
%         hold on
%         material metal
%         %    end
% end
% 
% 
% 
% 
% 
% % set(gca,'facelighting','none');
% 
%  z_scale_bar  = [-100 -100]; 
%  y_scale_bar  = [-100 -100];
%  x_scale_bar  = [-100 -100]; 
%  
% %  plot3(x_scale_bar, y_scale_bar, z_scale_bar,'Color','w','LineWidth',5)
% %  text(-97,-93,0, '30 $\mu m$','Interpreter','Latex','FontSize',10,'Color','w')
%  
%  
% %  z_scale_bar  = [0 0]; 
% %  y_scale_bar  = [-100 -70];
% %  x_scale_bar  = [100 100]; 
% %  
% % plot3(x_scale_bar, y_scale_bar, z_scale_bar,'Color','w','LineWidth',5)
% %  text(93.,-97,0, '30 $\mu m$','Interpreter','Latex','FontSize',10,'Color','w')
% %  
%  
% 
% 
% % xlabel('($\mu m$)','Interpreter','Latex','FontSize',15)
% % ylabel('($\mu m$)','Interpreter','Latex','FontSize',15)
% % zlabel('($\mu m$)','Interpreter','Latex','FontSize',15)
% 
% 
% 
% %  set(gca,'Color','k')
%  set(gca,'Color',[38/255 38/255 38/255])
% %  yticks([])
% %  xticks([])
% %  zticks([])
% 
% 
% 
% xlim([x_dom(1) x_dom(2)])
% ylim([y_dom(1) y_dom(2)]) 
% zlim([z_dom(1) z_dom(2)])
% xlabel('$X\;(\mu m)$','Interpreter','Latex','FontSize',15)
% ylabel('$Y\;(\mu m)$','Interpreter','Latex','FontSize',15)
% zlabel('$Z\;(\mu m)$','Interpreter','Latex','FontSize',15)
% 
% % title("Time =  " + t_total_h + " h " + t_total_min + " min")
% set(gca,'FontSize',10,'FontName','Times New Roman')


%% POSTPROCESS
% set(gca,'FontSize',20,'FontName','Times New Roman')
% xlabel('$X\;(\mu m)$','Interpreter','Latex','FontSize',20)
% ylabel('$Y\;(\mu m)$','Interpreter','Latex','FontSize',20)
% zlabel('$Z\;(\mu m)$','Interpreter','Latex','FontSize',20)

%% EVOLUTIONS 

% plot(t_t,n_ct,'ro','LineWidth',2)
% hold on 
% plot((1:length(n_c_time)).*(Dt_mech/(60*24)),n_c_time,'k','LineWidth',2) 
% 
% xlabel('$t\;(days)$','Interpreter','Latex','FontSize',20)
% ylabel('$N_{c}$','Interpreter','Latex','FontSize',20)
% 
% xticks([0 3 5 7])
% ylim([0 70])
% xlim([0 7])
% % legend('target values','Total number of cells')
% 
% set(gca,'FontSize',20,'FontName','Times New Roman')


%% STATE 

%    state_evo(n_c,evo) 
% 1=inactive  2=active
% 
% 
% for i = 1:counter_steps-1
% 
% %  n_evo_u = sum(state_evo(:,))
%  n_evo_q(i)    = sum(state_evo(:,i)==1); % 1=inactive  2=active
%  n_evo_p(i)    = sum(state_evo(:,i)==2); % 1=inactive  2=active
% 
% end
% 
% 
% plot((1:length(n_c_time)).*(Dt_mech/(60*24)),n_evo_p,'b','LineWidth',1) 
% hold on
% plot((1:length(n_c_time)).*(Dt_mech/(60*24)),n_evo_q,'g','LineWidth',1) 
% 
% 
% 
% 
% legend('Target values','Total number of cells','Number of proliferative cells',...
%        'Number of quiescent cells')


%% STATE EVO

% counter_evo = counter_steps; 
% evo_nc   = [1:n_c]'; 
% evo_time = [1:counter_evo-1]'.*Dt_mech/60 ; 
%  
% 
% 
% 
% figure
% imagesc(state_evo); 
% hs=imagesc(state_evo); 
% 
% %convierte xaxis de length a tiempo
% set(hs, 'XData', [1, counter_evo-1].*Dt_mech/60);
% xlim([0 (counter_evo-1)*Dt_mech/60])
% 
% set(gca,'YDir','normal') 
% set(gca,'FontSize',10.5,'FontName','Times New Roman')
% xlabel('$t\;(h)$','Interpreter','Latex','FontSize',15)
% ylabel('$N_{c}$','Interpreter','Latex','FontSize',15)
% 
% % cc = colorbar('Ticks',[0,1,2,3]);
% % % cc.YTickLabel = {'Unborn cell', 'Quiescence', 'Proliferation', 'Secretion'};
% % cc.YTickLabel = {'U', 'Q', 'P', 'S'};
% % caxis([0 3]) % sets colorbar limits 
% 
% map = [1 1 1
%        0 1 0
%        0 0 1
%       ];
% colormap(map)


%% 


% n_lines = (0:1:n_c)-0.5;
% 
% yline(n_lines)




%% 


% set(gca,'FontSize',25,'FontName','Times New Roman')
% % xlabel('$t\;(h)$','Interpreter','Latex','FontSize',40)
% 
% ylabel('$N_{c}$','Interpreter','Latex','FontSize',40) 
% 
% 
% xlabel('')
% set(gca,'xTicklabel',[])
% 
% 
% % ylim([0.5 60])
% % yticks([10 20 30 40 50 59])
% 
% 
% yticks([5 15 25 35 45 53 60])


%% 

% xlabel('')
% ylabel('')
% zlabel('')
% 
% xticks('')
% yticks('')
% zticks('')
% 


%% POST WITH NUCLEUS 


% %-------------------------------------------------------------------

% theta = 0:0.5:2*pi; 
% gamma = 0:0.5:2*pi;
% 
% clear x_ca 
% clear y_ca 
% clear z_ca 
% 
% figure
% counter_p = 1; 
% for c=1:n_c 
% %     if z_c(c)<=0
%     for i=1:length(theta) 
%         for j=1:length(gamma) 
%                 
%                 x_ca(counter_p) = x_c(c) + R_c(c)*sin(theta(i))*cos(gamma(j));
%                 y_ca(counter_p) = y_c(c) + R_c(c)*sin(theta(i))*sin(gamma(j));
%                 z_ca(counter_p) = z_c(c) + R_c(c)*cos(theta(i));
%                 
%                 counter_p = counter_p+1; 
%         end  
%     end
% %     end
% end
% 
% 
% shp1 = alphaShape(x_ca',y_ca',z_ca',20); 
% h1=plot(shp1);
% 
% set(h1,'FaceColor',[44/255 164/255 0/255], ...
%   'FaceAlpha',0.55,'EdgeColor','none','FaceLighting','gouraud')
% 
% 
%     x_dom = [-120 120];    %um 
%     y_dom = [-120 120];    %um 
%     z_dom = [-120 120];    %um
%     
% xlim([x_dom(1) x_dom(2)])
% ylim([y_dom(1) y_dom(2)]) 
% zlim([z_dom(1) z_dom(2)])
% 
% 
% 
% hold on
% for i = 1:n_c 
% 
% %     if z_c(i)<= 0
% 
% %Real magnitude    
% [xs,ys,zs] = sphere;     
% xs = xs*R_c(i)/2;
% ys = ys*R_c(i)/2; 
% zs = zs*R_c(i)/2; 
% 
% 
% hSurface=surf(xs+x_c(i),ys+y_c(i),zs+z_c(i));
%   set(hSurface,'FaceColor',[156/255 70/255 0/255], ...
%   'FaceAlpha',1,'FaceLighting','gouraud','EdgeColor','none')
% 
% 
% daspect([1 1 1]);
% 
% camlight
% hold on
% material metal
% 
% end
% 
%  
% 
% 
% 
%  z_scale_bar  = [0 0]; 
%  y_scale_bar  = [-100 -100];
%  x_scale_bar  = [-100 -70]; 
%  
%  plot3(x_scale_bar, y_scale_bar, z_scale_bar,'Color','w','LineWidth',5)
%  text(-97,-93,0, '30 $\mu m$','Interpreter','Latex','FontSize',10,'Color','w')
%  
%  
%  
% xlabel('','Interpreter','Latex','FontSize',15)
% ylabel('','Interpreter','Latex','FontSize',15)
% zlabel('','Interpreter','Latex','FontSize',15)
% 
%  set(gca,'Color','k')
%  yticks([])
%  xticks([])
%  zticks([])
%  




