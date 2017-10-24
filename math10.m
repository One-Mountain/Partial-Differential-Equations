%% Title 
% Heat and wave equations
%%  Initializations 
%{
Solving the differential equation with 
Utt+cUt=Uxx
U(0,t)=0, Ux(1,t)=0, U(x,0)=f(x), Ut(x,0)=0
where f(x)= 4x for 0< x < 1/4
          = 1  for 1/4< x < 1
          = dc  o/w
%}
clc; clear
m=800;                        %approximation
n=0:1:m;
Bet=(2*n+1)*pi/2;

c=[0, 1/2, 1, 2];             %arbitrary choices for c in: Utt+ cUt=Uxx
Lc=length(c);                 %length of c
x=0:.01:1;                    %Let that be our range for x-domain
Lx=length(x);
%% Utt+cUt=Uxx
for k=1:1:Lc                        %for each choice of c, plot it
    subplot(2,2,k)
    C1=c(k);
    a=sqrt(abs(C1^2-4*Bet.^2))/2;   %stuff inside the sin/cos;
    X=sin(Bet'*x);                  %X(x)
    Cn=8*sin(Bet/4)./(Bet.^2);      %Cn
    
    t=0:.1:2*pi/a(1);      %time domain
    Lt=length(t);          %length of t (aka pages of flip book)
    U1=zeros(1,Lx);        %initialize as zero.
    for r=1:1:length(n)
        U1(r,:)=Cn(r)*X(r,:);
    end
    U=sum(U1);             %plot f(x)=sum(Cn*X(x)) for t=1;
    plot(x,U);
    axis([0,1 -1,1])
    subplot(2,2,k)
    h=plot(x,U,'YDataSource','U');
    
    for w=1:1:Lt                    %animation
        tt=t(w);
        T=exp(-tt*C1/2)*(cos(a*tt)+C1./(2*a).*sin(a*tt)); %T(t)
        for r=1:1:length(n)
            U1(r,:)=Cn(r)*X(r,:)*T(r);
        end
        U=sum(U1);
        refreshdata(h,'caller');    %Update plot
        subplot(2,2,k)
        axis([0,1 -1,1])
        title(['Utt+',num2str(C1),'Ut=Uxx'])
        xlabel('x');
        drawnow; pause(.1)
    end
 end

%% Ut=Uxx
%Using the same initial conditions excep T'(0)=0
%then only T(t) changes
figure
m1=800;
n=1:1:m1;

x1=0:.01:1;
t1=0:.05:2;
Lx1=length(x1);
Lt1=length(t1);

U2=zeros(1,Lx1);
for r1=1:1:length(n)
    U2(r1,:)=Cn(r1)*X(r1,:);
end
Un=sum(U2);
plot(x1,Un);
axis([0,1 -.1,1])

 h=plot(x1,Un,'YDataSource','Un');
  for ww=1:1:Lt1
        tn=t1(ww);
        T1=exp(-(Bet).^2*tn);
        for r1=1:1:length(n)
            U2(r1,:)=Cn(r1)*X(r1,:)*T1(r1);
        end
        Un=sum(U2);
        refreshdata(h,'caller');
        axis([0,1 -.1,1])
        title('Ut=Uxx')
        xlabel('x')
        drawnow; pause(.5)
  end
  