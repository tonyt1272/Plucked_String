clear
clc
close all
close hidden
%This simulation is an ideal string with damping.
%
j=1i;
h=.05;     %Initial displacement of the string
L=1.5;   %Length of the string
T=100;    %Tension in the string (Newtons)
rho_L=.1; % length density, (Kg/meter)
Beta=10;
c=sqrt(T/rho_L); %wave speed on the string (meters/second)
prompt='Hit return when ready';
dx=.0005;
x=0:dx:L+dx;% 

% y1=2*h/L*x(1:length(x)/2);      %Setting initial displacement
% y2=h-2*h/L*x(1:length(x)/2);    %Plucked dead center

y1=8*h/L*x(1:(length(x)/8));          %Setting initial displacement 
y2=h-8/7*h/L*x(1:(length(x)*7/8));    %Plucked 1/4 from the left

y_x_0=[y1 y2 0];
y_x_total=[];

n_modes=30; %sets how high the mode number goes
for n = 1:n_modes        %sets how high the mode number goes
        mode(n,:)=sin(pi*n*x/L);
        A(n)=2/L*sum(y_x_0.*mode(n,:))*dx;  %Mode amplitude
        omega_d(n)=sqrt(((pi*n/L)*c)^2-Beta^2)+j*Beta;
        %damped frequency with decay applied to modes individually
end


duration=.25; %time in seconds
dt=1/2400;%.0002;
str_wait=waitbar(0,'data progress');
for t=0:dt:duration        %end value is the total time in seconds
    for n = 1:n_modes        
            if t==0
            modef_init(n,:)=A(n)*mode(n,:)*exp(j*omega_d(n)*t); 
            y_x_init=sum(modef_init); %series approximation of initial 
                                      %condition
            end
          modef(n,:)=A(n)*mode(n,:)*exp(j*omega_d(n)*t);    %each mode's  
                                                            %individual 
                                                            %response
    waitbar(t/duration,str_wait,'data progress');
    end 
    
    y_x=sum(modef);  %the state of the string at time t
    y_x_total=[y_x_total; y_x]; %collecting the entire time response in a 
                                %matrix
end
close(str_wait)

for i =1:t/dt+1                         %Real part of the complex 
y_x_total(i,:)=real(y_x_total(i,:));    %displacement is the real solution
end                                     %for y(x,t), the transverse
                                        %displacement of the string.
  movie_progress=waitbar(0,'video progress');                                      
for i =1:t/dt+1                         
  if i==1
    figure(1)
    hFig1 = figure(1);
    set(hFig1, 'Position', [200 300 740*2 455*.5])
    plot(x,(y_x_total(i,:)))
    title(['Plucked String, Damped, \beta=',num2str(Beta)])
    ylabel('Transverse displacement (y)')
    xlabel('position(x)')
    ylim([-1.3*h,1.3*h])
    xlim([0,L])
    grid on
    str = input(prompt,'s');
    tic
  end
    plot(x,(y_x_total(i,:)))
    ylim([-1.8*h,1.3*h])
    xlim([0,L])
    pause(.01)
    waitbar(i/(t/dt+1),movie_progress,'movie progress'); 
end
close(movie_progress)
toc   
figure(2);
hFig2 = figure(2);
set(hFig2, 'Position', [300 300 740*2 455*.5])
plot(x,y_x_0,x,y_x_init)
title(['Initial String Displacement, Fourier Series (',num2str(n_modes),' modes) vs Actual'])
legend('y_0','Fourier y_0')
ylabel('Transverse displacement (y)')
xlabel('position(x)')
ylim([-1.3*h,1.3*h])
xlim([0,L])
grid on
