%%% Produces the 9 graphs of Figure 3 & 4 and all country figures in the Appnedixfor each country
%%% Countrry codes are given in the file Country_Codes excel file 
clear all; 
warning('off','all')

%time_data=xlsread('testdata.xlsx');
time_data=xlsread('Incubation_death_6_20.xlsx');
time_data2=xlsread('Incubation_inf_6_20.xlsx');

fnameF = 'F:\COVID19_Manuscript_Repository - Copy\Appendix';

for n=1:127
    
%n=113
    
%    n
%n =input('Enter a number: ');
sn=string(n);

tvals=time_data(ismember(time_data(:,1), n), :);
tvals2=time_data2(ismember(time_data2(:,1), n), :);

%Incub=tvals(4);
%if Incub==6.5
%   Incub=6.0;
%end1
%Infec_duration=tvals(5);

%%%%%%
delIe=tvals2(2);
delIs=delIe+tvals2(3);
%delhE=Incub;
%delhS=delhE+Infec_duration;

%%%%%%%
delhE=tvals(2);%10;8;%10;8;10;%7;    
delhS=delhE+tvals(3);%20;22;18;%7;
%%%%%%%
%delIe=2.5;6;%%%3.5;3.5; 4;4.5; 4.76; %5;
%delIs=delIe+10;%8; %10;%%%9;10;% 8; 10;   %8;
%%%%%%%%%%%%%%%%%%%%%
sh=105; 

[num country N beta theta mu k eps]=CaseNo_3(n);
global T1 data Hdata h S00 del S00 O00 D00 N sh delhS delhE delIs delIe fac m  V  
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 3.2]);
x_width=3 ;y_width=3.2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
sn=string(n);
%axx = gca;
%exportgraphics(ax,'barchartaxes.png','Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu=0;  
LB=sqrt([0.0001 0.0001 0.15  0.0001 0.0001  0.0001  0.4]);
UB=sqrt([6        1    0.6   1       4     10000   0.8]);
LB=([0.0001 0.0001 sqrt(0.15)  0.0001 0.0001  0.0001  sqrt(0.4)]);
UB=([6        1    sqrt(0.6)   1       4     10000    sqrt(0.9)]);
nu=0;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pI=100;
fac=1;
nu1=0.5;
m=1;
fac=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Table=array2table(num, 'VariableNames',{'Day'	'Outd'	'Home'	'Inf'	'Death' 'Vacc' 'FVacc'});
%%%%%%%%%%%%%%%%
t=Table.Day(1:end-sh+1);
theta_H=max(Table.Home(sh:end),0.0001); 
Id=max(Table.Inf(sh:end),0);
Dd=max(Table.Death(sh:end),0);
data(:,1) = t;
data(:,2) = Id;
data(:,3) = Dd;
data(:,4) =max((Table.Home(sh:end)),0.0001);
H=data(:,4);
Hdata=max((Table.Home),0);%,0.001);
V=max(Table.FVacc,0);
%%%%%%%%%%
T1=length(H);
h=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O(1)=mean(Table.Inf(sh:sh+4)); %
pI=2*max(O(1),1);
I(1)=O(1)+pI;
D(1)=Dd(1); % death
S(1)=N-sum(Table.Inf(1:sh))-sum(Table.Death(1:sh)) -nu1*sum(V(1:sh-14));
%%%%%%%%%%%%
S00=S(1);
O00=O(1);
D00=D(1);
%%%%%%%%
x0=([beta, theta, eps, mu, k, pI nu1]);%
beta0 = sqrt(x0);
y1=(data(:,2));
y2=(data(:,3));
xD=(1:T1);
yD2(:,1)=y1;
yD2(:,2)=y2;
options = optimset('Display','off','MaxFunEvals',10000000000000);
[x0,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@modelR10,beta0,xD,yD2,LB,UB, options);
CIX=x0.^2;
conf = nlparci(CIX,residual,'jacobian',jacobian);
RSS=sum(sum((residual.^2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=x0(1)^2;
gamma=beta;
theta=x0(2)^2;
eps=x0(3)^2;
mu=x0(4)^2;
k=x0(5)^2;
pI=x0(6)^2;
nu1=x0(7)^2;
nu=nu1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df=5;
betaL=max(conf(1,1),0);
betaU=min(conf(1,2),UB(1)^2);
thetaL=max(conf(2,1),0);
thetaU=min(conf(2,2),UB(2)^2);
epsL=max(conf(3,1),0);
epsU=min(conf(3,2),1);
muL=max(conf(4,1),0);
muU=min(conf(4,2),1);
kL=max(conf(5,1),0);
kU=min(conf(5,2),UB(5)^2);
nu1L=max(conf(7,1),LB(7));
nu1U=min(conf(7,2),UB(7));



clear S O I D Ef2 R0 Ir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES
S(1)=S00;
O(1)=O00;
I(1)=O00+pI;
D(1)=D00;
for i=1:length(xD)-1 
Ef2(i)=max(fac*(mean(Hdata(sh+i))/100)^k,0);
Ef2(i)=min(fac*(mean(Hdata(sh+i))/100)^k,1);
R0(i)=max((beta*S(i)*max((1-theta*Ef2(i)),0)/N -eps),0);
Ir(i)=max((mean(I(max(1,i-delIs):max(1,i-delIe)))*max((1-theta*Ef2(i)),0)),1);
S(i+1)=S(i) -beta*S(i)*Ir(i)*h/N -nu1*V(max(sh-14+i,1));
I(i+1)=beta*S(i)*Ir(i)*h/N -eps*mean(I(max(1,i-delIs):max(1,i-delIe)))*h; 
O(i+1)=eps*mean(I(max(1,i-delIs):max(1,i-delIe)))*h;
D(i+1)=mu*mean(O(max(1,i-delhS):max(1,i-delhE)))^m;
end

H0=min((100/theta)*(1-(1+eps)/beta)^(1/k),100); 
Vsum=sum(V);
Vsum_per=min((sum(V)/N)*100,100);
ss=length(xD)-1;
NLL=(ss/2)*log(RSS/ss)+(ss/2)*log(2*pi)+(ss/2);
R0_cur=R0(length(xD)-1);
S_Over_N=S./N;
R0_H(:,1)=beta.*S_Over_N'.*(1-theta.*(Hdata(sh:sh+T1-1)./100).^k)-eps;
R0_H_U(:,1)=betaU.*S_Over_N'.*(1-thetaU.*(Hdata(sh:sh+T1-1)./100).^kU)-epsL;
R0_H_L(:,1)=betaL.*S_Over_N'.*(1-thetaL.*(Hdata(sh:sh+T1-1)./100).^kL)-epsU;
R0_MN=mean(R0_H(end-7:end,1));



c1=[0.4940 0.1840 0.5560]; %purple
c2=[0.4660 0.6740 0.1880]; %green
c3=[0.8500 0.3250 0.0980]; %red
c4=[0 0.4470 0.7410]; %blue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% PLOTS
figure(1); clf; hold on
u1=plot(((O)),'-','LineWidth',2, 'color', c3); 
u2=plot(((Id)),'o-','MarkerSize',4, 'color', c1) ; 
title(country);
b1=legend([u1,u2],'Model-3 (Home&Vac)','New Cases-Data','FontSize', 10,'Location','northwest');
set(b1,'Box','off')
axis([0, length(H), 0, inf]);
box on
ylabel('Number of New Cases (daily)', 'FontSize', 12)
xlabel('Day since April 15, 2021', 'FontSize', 12)

%%%%%%%%%%%%%%%%%
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 3.2]);
x_width=3 ;y_width=3.2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    
filename1 = sprintf('Fig_%s%s','1_',sn);
saveas(gcf, fullfile(fnameF, filename1), 'fig');
saveas(gcf, fullfile(fnameF, filename1), 'pdf');
%saveas(gcf, fullfile(fnameF, filename1), 'jpg');
fn=fullfile(fnameF, filename1);
print(gcf,fn,'-djpeg','-r300')
%saveas(gcf, pr, 'jpg')
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
figure(2); clf; hold on
u1=plot([1:length(D)], D(1:length(D)),'-','LineWidth',2, 'color', c3); 
u2=plot(((Dd)),'o-','MarkerSize',4, 'color', c1); 
title(country);
b3=legend([u1,u2],'Model-3 (Home&Vac)','Deaths-Data','FontSize', 10,'Location','northwest');
set(b3,'Box','off')
axis([0, length(H), 0, Inf]);
box on
ylabel('Number of Deaths (daily)', 'FontSize', 12)
xlabel('Day since April 15, 2021', 'FontSize', 12)


%%%%%%%%%%%%%%%%%
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 3.2]);
x_width=3 ;y_width=3.2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    
filename1 = sprintf('Fig_%s%s','2_',sn);
saveas(gcf, fullfile(fnameF, filename1), 'fig');
saveas(gcf, fullfile(fnameF, filename1), 'pdf');

fn=fullfile(fnameF, filename1);
print(gcf,fn,'-djpeg','-r300')

%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
figure(3); clf; hold on
yyaxis right
p1=plot(H,'square','MarkerSize',4, 'color', c2); 
plot(H,'-','LineWidth',0.5, 'color', c2);
ylabel('Home-stay (% diff. from pre-Covid Normal)', 'FontSize', 12,'color', c2)
xlabel('Day since April 15, 2021', 'FontSize', 12)
ax = gca;
ax.YColor = c2;
axis([0, length(H), 0, 50]);
yyaxis left
cVc=(cumsum(V)+1);
p2=plot(cVc(sh+1:end)/N*100,'<','MarkerSize',4, 'color', c3); 
plot(cVc(sh+1:end)/N*100,'-','LineWidth',0.5, 'color', c3); 
ylabel('Cumulative Numb. Vacc (% of Populaton)', 'color', c3)
xlabel('Day since April 15, 2021', 'FontSize', 12)
ax = gca;
ax.YColor = c3;
b1=legend([p1, p2],'Home-Stay Data', 'Vaccination Data','FontSize', 10,'Location','northwest');
set(b1,'Box','off')
axis([0, length(H), 0, 50]);
box on
title(country);

%%%%%%%%%%%%%%%%%
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 3.2]);
x_width=3 ;y_width=3.2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    
filename1 = sprintf('Fig_%s%s','3_',sn);
saveas(gcf, fullfile(fnameF, filename1), 'fig');
saveas(gcf, fullfile(fnameF, filename1), 'pdf');

fn=fullfile(fnameF, filename1);
print(gcf,fn,'-djpeg','-r300')

%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%
figure(4); clf; hold on
title(country);
%a1=plot(beta.*(1-theta.*Ef2),'-','LineWidth',2, 'color', c3); 
a2=plot(R0,'-','LineWidth',2, 'color', c4); 
%%%%%
xl = linspace(0, length(H))';
yl = [0*(xl)+1];
a3=line(xl,yl,'Color','black','LineStyle','--');

%if betaL>=betaU
%a22=plot(R0_H_U(1:length(R0)),':','LineWidth',2, 'color', c1); 
%a32=plot(max(-R0_H_U(1:length(R0)),0),':','LineWidth',2, 'color', c1); 
a22=plot(R0_H_U,':','LineWidth',2, 'color', c1); 
a32=plot(max(R0_H_L,0),':','LineWidth',2, 'color', c1); 

%a32=plot(max(R0_H_L(1:length(R0)),0),':','LineWidth',2, 'color', c1); 
%else
%a22=plot(R0_H_U(1:length(R0))-R0',':','LineWidth',2, 'color', c1); 
%a32=plot(max(R0_H_L(1:length(R0))+R0',0),':','LineWidth',2, 'color', c1);     
%end

%%%%
ylabel('R0(H,V)','FontSize', 12)
xlabel('Day since April 15, 2021', 'FontSize', 12)
b1=legend([a2,a22,a3],'R0(t)(Model 3)','R0(t) 95% CI', 'FontSize', 10,'Location','northwest');
set(b1,'Box','off')
axis([0, length(H), 0, 8]);
box on

%%%%%%%%%%%%%%%%%
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 3.2]);
x_width=3 ;y_width=3.2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    
filename1 = sprintf('Fig_%s%s','4_',sn);
saveas(gcf, fullfile(fnameF, filename1), 'fig');
saveas(gcf, fullfile(fnameF, filename1), 'pdf');

fn=fullfile(fnameF, filename1);
print(gcf,fn,'-djpeg','-r300')

%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%
%figure(5); clf; hold on
%title(country);
%a1=plot(H(2:end), beta.*(1-theta.*Ef2), 'o','MarkerSize',6, 'LineWidth',2, 'color', c3); ; 
%a2=plot(H(2:end), R0, 'o','MarkerSize',6,  'LineWidth',2,'color', c4);
%VumV=cumsum(V(1:end));
%sz =round(max((100*(VumV(sh+1-14:end-14)./N)),1));% 25;
%c = linspace(1,sz(end),length(H(2:end)));
%scatter(H(2:end), R0',sz, c,'filled');
%hcb = colorbar;
%hcb.Title;sz;
%hcb.Title.String = "Vc(-14d)%";
%xl = linspace(0,126)';
%yl =[0*(xl)+1];
%line(xl,yl,'Color','black','LineStyle','--');
%a3=line(xl,yl,'Color','black','LineStyle','--');
%%%%
%ylabel('R0(H,V)','FontSize', 12)
%xlabel('Home-stay (H%) (daily)', 'FontSize', 12)
%%%%%
%xl2 = [H(2:end)]';
%yl2 = [beta.*(1-theta.*Ef2)];
%pv = polyfit(xl2,yl2,2);
%x1 = linspace(0,50);
%y1 = polyval(pv,x1);
%axis([0, 50, 0, 4]);
%hold on
%xl2 = [H(2:end)]';
%yl2 = [R0];
%pv = polyfit(xl2,yl2,2);
%x1 = linspace(0,50);
%y1 = polyval(pv,x1);
%b1=legend([a2,a3],'R0 (Model 3)', 'R0=1','FontSize', 10,'Location','northeast');
%set(b1,'Box','off')
%axis([0, 50, 0, 4]);
%box on

%%%%%%%%%%%%%%%%%
%set(gcf, 'PaperUnits', 'inches');
%set(gcf, 'PaperSize', [3 3.2]);
%x_width=3 ;y_width=3.2;
%set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    
%filename1 = sprintf('Fig_%s%s','5_',sn);
%saveas(gcf, fullfile(fnameF, filename1), 'fig');
%saveas(gcf, fullfile(fnameF, filename1), 'pdf');

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
figure(5); clf; hold on
title(country);
a1=plot(H(2:end), beta.*(1-theta.*Ef2), 'o','MarkerSize',6, 'LineWidth',2, 'color', c3); ; 
a2=plot(H(2:end), R0, 'o','MarkerSize',6,  'LineWidth',2,'color', c4);
VumV=cumsum(V(1:end));
sz =round(max((100*(VumV(sh+1-14:end-14)./N)),1));% 25;
c = linspace(1,sz(end),length(H(2:end)));
scatter(H(2:end), R0',sz, c,'filled');
hcb = colorbar;
hcb.Title;sz;
hcb.Title.String = "Vc(-14d)%";
xl = linspace(0,126)';
yl =[0*(xl)+1];
line(xl,yl,'Color','black','LineStyle','--');
a3=line(xl,yl,'Color','black','LineStyle','--');
%%%%
ylabel('Infection rate \beta_a(H,V) & R0(H,V)','FontSize', 12)
xlabel('Home-stay (H%) (daily)', 'FontSize', 12)
%%%%
xl2 = [H(2:end)]';
yl2 = [beta.*(1-theta.*Ef2)];
pv = polyfit(xl2,yl2,2);
x1 = linspace(0,50);
y1 = polyval(pv,x1);
axis([0, 50, 0, 4]);
hold on
xl2 = [H(2:end)]';
yl2 = [R0];
pv = polyfit(xl2,yl2,2);
x1 = linspace(0,50);
y1 = polyval(pv,x1);
b1=legend([a1,a2,a3],'\beta_a=\gamma(1-Ef)', 'R0 (Model 3)', '\beta,R0=1','FontSize', 10,'Location','northeast');
set(b1,'Box','off')
axis([0, 50, 0, 4]);
box on

%%%%%%%%%%%%%%%%%
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 3.2]);
x_width=3 ;y_width=3.2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    
filename1 = sprintf('Fig_%s%s','5_',sn);
saveas(gcf, fullfile(fnameF, filename1), 'fig');
saveas(gcf, fullfile(fnameF, filename1), 'pdf');

fn=fullfile(fnameF, filename1);
print(gcf,fn,'-djpeg','-r300')

%%%%%%%%%%%%%
figure(8); clf; hold on
title(country);
VumV=cumsum(V(1:end));
zV1=100.*VumV(sh+1-14:end-14)./N;
zH1=H(2:end);
zR1=R0;
surf([zV1(:) zV1(:)], [zR1(:) zR1(:)], [zH1(:) zH1(:)], ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', 8);            % Make a thicker line
view(2);  
colorbar; 
hcb = colorbar;
hcb.Title;sz;
hcb.Title.String = "H%";
xl = linspace(0,length(xD)-1)';
yl =[0*(xl)+1];
line(xl,yl,'Color','black','LineStyle','--');
a3=line(xl,yl,'Color','black','LineStyle','--');
%%%%
ylabel('R0(H,V)','FontSize', 12)
xlabel('Vc%', 'FontSize', 12)
%%%%
xl2 = [H(2:end)]';
yl2 = [beta.*(1-theta.*Ef2)];
pv = polyfit(xl2,yl2,2);
x1 = linspace(0,50);
y1 = polyval(pv,x1);
hold on
box on
axis([0, 60, 0, 4]);


%%%%%%%%%%%%%%%%%
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 3.2]);
x_width=3 ;y_width=3.2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    
filename1 = sprintf('Fig_%s%s','8_',sn);
saveas(gcf, fullfile(fnameF, filename1), 'fig');
saveas(gcf, fullfile(fnameF, filename1), 'pdf');

fn=fullfile(fnameF, filename1);
print(gcf,fn,'-djpeg','-r300')

%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%
figure(6); clf;
Vcc=[0, 50000];
Hc=[0,10,20];  
for j=1:1
    Vc=Vcc(j);
    for i=1:3
    Hx=Hc(i);    
for i=T1:T1+220 ; %   
Ef2(i)=fac*(Hx/100)^k;
R0(i)=max((beta*S(i)*max((1-theta*Ef2(i)),0)/N -eps),0);
Ir(i)=max((mean(I(max(1,i-delIs):max(1,i-delIe)))*max((1-theta*Ef2(i)),0)),0);
S(i+1)=S(i) -beta*S(i)*Ir(i)*h/N -nu1*Vc;%
I(i+1)=beta*S(i)*Ir(i)*h/N -eps*mean(I(max(1,i-delIs):max(1,i-delIe)))*h;
O(i+1)=eps*mean(I(max(1,i-delIs):max(1,i-delIe)))*h; 
%if (i-del)>0
%D(i+1)=mu*mean(O(i-del:i))^m; 
%else
%D(i+1)=mu*mean(O(1:i))^m; 
%end
%len=abs(1:min(i,delhS));
%y = lognpdf(len,pm,pvv) ;%plot([0:40],y)
%D(i+1)=mu*mean(O(max(1,i-delhS+1):i).*flip(y ,2)); 
D(i+1)=mu*mean(O(max(1,i-delhS):max(1,i-delhE)))^m;
end;
S=max(S,1);
I=max(I,1);
O=max(O,1);
D=max(D,1);
%%%%%%%%%%%%%%
figure(6); hold on
plot(log(cumsum(I)),'LineWidth',2); hold on
xlabel('Day since 15th April 2021', 'FontSize', 12)
ylabel('Cumulative Number Infected (Log)', 'FontSize', 12)
b1=legend('H0%,V=0','H10%,V=0','H20%,V=0','H0%,V=5E+5','H10%,V=5E+5','H20%,V=5E+5');%,'Ef2') 
set(b1,'Box','off')
axis([0, 365, 0, inf]);
box on
    end
end

for j=2:2
    Vc=Vcc(j);
 for i=1:3
    Hx=Hc(i);    
for i=T1:T1+220   
Ef2(i)=fac*(Hx/100)^k;
R0(i)=max((beta*S(i)*max((1-theta*Ef2(i)),0)/N -eps),0);
Ir(i)=max((mean(I(max(1,i-delIs):max(1,i-delIe)))*max((1-theta*Ef2(i)),0)),0);
if S(i)>(beta*S(i)*Ir(i)*h/N +nu1*Vc);
S(i+1)=S(i) -beta*S(i)*Ir(i)*h/N -nu1*Vc;%
else
S(i+1)=S(i) -beta*S(i)*Ir(i)*h/N;%    
end
I(i+1)=beta*S(i)*Ir(i)*h/N -eps*mean(I(max(1,i-delIs):max(1,i-delIe)))*h; 
O(i+1)=eps*mean(I(max(1,i-delIs):max(1,i-delIe)))*h; 
%if (i-del)>0
%D(i+1)=mu*mean(O(i-del:i))^m; 
%else
%D(i+1)=mu*mean(O(1:i))^m; 
%end
%len=abs(1:min(i,delhS));
%y = lognpdf(len,pm,pvv) ;%plot([0:40],y)
%D(i+1)=mu*mean(O(max(1,i-delhS+1):i).*flip(y ,2)); 
D(i+1)=mu*mean(O(max(1,i-delhS):max(1,i-delhE)))^m;
end;

S=max(S,1);
I=max(I,1);
O=max(O,1);
D=max(D,1);

%%%%%%%%%%%%%%
figure(6); hold on
plot(log(cumsum(I)),'--','LineWidth',2); hold on
xlabel('Day since 15th April 2021', 'FontSize', 12)
ylabel('Cumulative Number Infected (Log)', 'FontSize', 12)
b1=legend('H0%,V=0','H10%,V=0','H20%,V=0','H0%,V=5E+5','H10%,V=5E+5','H20%,V=5E+5');%,'Ef2') 
set(b1,'Box','off')
axis([0, 365, 0, inf]);
box on
    end
end

%%%%%%%%%%%%%%%
figure(6); hold on
yl = linspace(0,20)';
xl =[0*(yl)+T1];
a3=line(xl,yl,'Color','black','LineStyle','--');
legend('H0%,    V=0','H10%,   V=0','H20%,   V=0','H0%,   V=5E+5','H10%, V=5E+5','H20%, V=5E+5', 'Simu starting day','Location','southeast');%,'Ef2') 
title(country);

%%%%%%%%%%%%%%%%%
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 3.2]);
x_width=3 ;y_width=3.2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    

filename1 = sprintf('Fig_%s%s','6_',sn);
saveas(gcf, fullfile(fnameF, filename1), 'fig');
saveas(gcf, fullfile(fnameF, filename1), 'pdf');

fn=fullfile(fnameF, filename1);
print(gcf,fn,'-djpeg','-r300')

%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%
figure(7);clf; hold on
con=1;   pr=1;
for i=1:4
    Hy=0;
for j=1:11 
yR(j)=beta*(pr)*(1-theta*(Hy/100)^k)-eps;
Hy=Hy+10;
end
pr=pr-0.25;
figure(7); hold on
plot([0:10:100],yR','LineWidth',1.5); hold on
xlabel('Home-stay H%', 'FontSize', 12)
ylabel('R0', 'FontSize', 12)
b1=legend('Vp=0%','Vp=25%','Vp=50%','Vp=75%');%,'Ef2') 
set(b1,'Box','off')
axis([0, 100, 0, 2]);
box on
end
xl = linspace(0,126)';
yl =[0*(xl)+1];
line(xl,yl,'Color','black','LineStyle','--');
b1=legend('Vp=0%','Vp=25%','Vp=50%','Vp=75%', 'R0=1');%,'Ef2') 
title(country);


%%%%%%%%%%%%%%%%%
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 3.2]);
x_width=3 ;y_width=3.2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

filename1 = sprintf('Fig_%s%s','7_',sn);
saveas(gcf, fullfile(fnameF, filename1), 'fig');
saveas(gcf, fullfile(fnameF, filename1), 'pdf');

fn=fullfile(fnameF, filename1);
print(gcf,fn,'-djpeg','-r300')

%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%
figure(9); clf; hold on
title(country);
VumV=cumsum(V(1:end));
zV=100.*VumV(sh+1-14:end-14)./N;
zH=H(2:end);
zR=R0(1:length(H)-1);
hs=surf([zV(:) zV(:)], [zH(:) zH(:)], [zR(:) zR(:)], ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', 5);            % Make a thicker line
hcb = colorbar;
hcb.Title;sz;
hcb.Title.String = "R0(H,V)";
hold on
for i=1:length(zR)
if zR(i)<1.2 & zR(i)>0.98 
    dd=plot3(zV(i), zH(i), zR(i),'o','MarkerSize', 8, 'color', 'red', 'LineWidth',2); hold on;
end
end
xlabel('Vc%', 'FontSize', 12)
ylabel('H%','FontSize', 12)
set(b1,'Box','off')
text( 10 ,  47 , 1 , 'R0(H,V)=1','Color', 'k', 'FontSize',10);
text( 5 , 47 , 1 , '0','Color', 'r', 'FontSize',10);
title(country);
box on
axis([0, 60, 0, 50]);
%end

%%%%%%%%%%%%%%%%%
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [3 3.2]);
x_width=3 ;y_width=3.2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

filename1 = sprintf('Fig_%s%s','9_',sn);
saveas(gcf, fullfile(fnameF, filename1), 'fig');
saveas(gcf, fullfile(fnameF, filename1), 'pdf');

fn=fullfile(fnameF, filename1);
print(gcf,fn,'-djpeg','-r300')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%


end
