% SEIaIhAHRD-IR model
% Variables observadas: Ia, Ih, A, H, Rc, Dc, Cc
% Variables no observadas: ...
% Parametros a estimar: 12
%%
clear model data params options

% Leer datos de curva epidemica
%

load datos_20200930.mat data

%Tiempo(inicia en 1)
% InfecciososAmb, InfecciososHosp, AisladosAmb, AisladosHosp, 
% RecuperadosAcumulados, MuertesAcumuladas, ContagiosAcumulados
% DatosEpi=readmatrix('../DatosCurvaEpidemicaColima.csv');
% DatosEpiVdeA=readmatrix('../DatosCurvaEpidemicaVilladeAlvarez.csv');
% 
% data.ydata=DatosEpi(:,2:8)+DatosEpiVdeA(:,2:8);
% data.ylabels={'Prevalencia','Prevalencia','Prevalencia','Prevalencia',...
%     'Acumulados','Acumulados','Acumulados'};
% data.titles={'Infecciosos Amb','Infecciosos Hosp','Aislados Amb','Aislados Hosp',...
%     'Recuperados','Defunciones','Contagios acumulados'};
Ttot=length(data.ydata(:,1));

% data.xdata=DatosEpi(:,1)-1;
% data.xlabels={'dias (a partir del 24/03)'};

diasNum=[8,22,38,52,69,83,99,113,130,144,161,175,191,205,222,236];
diasEtiq={'1/abr','15/abr','1/may','15/may','1/jun','15/jun','1/jul','15/jul',...
    '1/ago','15/ago','1/sep','15/sep','1/oct','15/oct','1/nov','15/nov'};

% NColima=172513; %Colima
% NVdeA=150479; %Villa de Alvarez
% N=NColima+NVdeA;
% Ic0=data.ydata(1,2); %contagiosos confirmados
% I0=30; %contagiosos sin confirmar -- parametro por estimar
% % S,E,Iac,Ihc,Ac_(1-4),Hc_(1-4),HDc_(1-4),Rc,Dc,Cc,I,R,C
% data.y0=[N-Ic0-I0,0,Ic0,0,zeros(1,4),zeros(1,4),zeros(1,4),0,0,Ic0,I0,0,I0]; 

% Change working directory
cd('../../../mcmctoolbox')

%%
% The model sum of squares in function |fermentass| is
% given in the model structure.
model.ssfun = @funSS;

%% Prior distributions of parameters
% 'name', start, [min,max], N(mu,s^2) 
params = {
% positive parameters
{'beta',1,0,2}
{'factorConf',40,20,70}
{'pA',0.5,0,1}
{'sigma',0.21,0.1,0.3}
{'alphaAC',0.3,0,1}
{'alphaHC',0.3,0,1}
{'gammaAC',0.1,0,1}
{'gammaHC',0.1,0,1}
{'muC',0.07,0.01,0.15}
{'gamma',0.5,0,1}
{'I0',30,20,60}
{'pD',0.5,0,1}
};

%%
% We assume having at least some prior information on the
% repeatability of the observation and assign rather non informational
% prior for the residual variances of the observed states. The default
% prior distribution is sigma2 ~ invchisq(S20,N0), the inverse chi
% squared distribution (see for example Gelman et al.). Each
% component may have separate variances.
model.S20 = [1 1 1 1 1 1 1];
model.N0  = [4 4 4 4 4 4 4];

%%
% First generate an initial chain.
options.nsimu = 15000; %VARY ACCORDINGLY TO EACH EXAMPLE
[results, chain, s2chain]= mcmcrun(model,data,params,options);
figure
mcmcplot(chain,[],results,'chainpanel')

%%
% Then re-run starting from the results of the previous run,
% this will take several minutes.
options.nsimu = 20000; %VARY ACCORDINGLY TO EACH EXAMPLE
[results, chain, s2chain] = mcmcrun(model,data,params,options, results);

%%
% Chain plots should reveal that the chain has converged and we can
% use the results for estimation and predictive inference.
figure
mcmcplot(chain,[],results,'pairs');
%%%%%%%%%%%%%%%%%%%%%%%%
figure
mcmcplot(chain,[],results,'denspanel',3); % change mcmcpredplot.m (lines 55-56) if big number of parameters
%%%%%%%%%%%%%%%%%%%%%%%%
figure
mcmcplot(chain,[],results,'chainpanel')

%%
% Function |chainstats| calculates mean and std from the chain and
% estimates the Monte Carlo error of the estimates. Number |tau| is
% the integrated autocorrelation time and |geweke| is a simple test
% for a null hypothesis that the chain has converged.
stats=chainstats(chain,results);

%%
% Plot results with mean of estimates
theta=stats(:,1)';
%disp(theta) % show parameter values (mean of distribution)
ss = funSS(theta,data);
rss=sum(ss);
%disp(rss) % show sum of squares (value of error) for such parameters values

% % Akaike information criterion
% n=prod(size(data.ydata))
% k=length(theta)
% AIC=n*log(rss/n)+2*k+2*k*(k+1)/(n-k-1)

timeorig=data.xdata;
time=linspace(data.xdata(1),data.xdata(end),100)'; % create time vector (more detail than in data)
%time=linspace(data.xdata(1),240,500)'; % create time vector (more detail than in data)
y0=data.y0;
y0(20)=theta(11);
y0(22)=theta(11);
%[t,y] = ode45(@modelSystem,time,y0,[],theta); % numerical solution of model with such parameters
ymod=modelSolution(time,theta,y0);
figure
for i=1:7 
    subplot(3,3,i)
    plot(data.xdata,data.ydata(:,i),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3); 
    hold on
    plot(time,ymod(:,i),'k') % plot numerical solution of model
    xlabel('dias (a partir del 24/03)');
    ylabel(data.ylabels(i)); 
    title(data.titles(i));    
    hold off
end

% subplot(2,3,6)
% plot(t,y(:,9),'--b') % plot numerical solution of model
% hold on
% plot(t,y(:,10),'b')
% plot(t,y(:,9)+y(:,10),'k')
% xlabel('dias (a partir del 24/03)');
% ylabel('Acumulados')
% title('Contagios totales')
% hold off


%% Predictive plots

% We need to augment the time in ydata to get more plot points
% |datamerge| in the toolbox does this.
data.xdata = datamerge(data.xdata,time);

% In order to use the |mcmcpred| function we need
% function |modelfun| with input arguments given as
% |modelfun(xdata,theta)|. We construct this as an anonymous function.
modelfun = @(d,th) modelSolution(d,th,[data.y0(1:19),th(11),data.y0(21),th(11)]);

% We sample 500 parameter realizations from |chain| and |s2chain|
% and calculate the predictive plots.
nsample = 1000;
%%% If s2chain is given to mcmcpred, then the plot shows 95% probability limits 
%%% for new observations and for model parameter uncertainty.
% out = mcmcpred(results,chain,s2chain,data.xdata,modelfun,nsample); %<--NOT SELECTED
%%% If s2chain is not used then the plot contains 50%, 90%, 95%, and 99% predictive 
%%% probability limits due parameter uncertainty.
out = mcmcpred(results,chain,[],data.xdata,modelfun,nsample); % prediction plots
figure
mcmcpredplot(out);
% add the data observations to the plot
hold on
for i=1:7 % depends on number of variables
  subplot(3,3,i)
  hold on
  plot(timeorig,data.ydata(:,i),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3);
  ylabel(data.ylabels(i)); 
  title(data.titles(i));    
  %xlabel('dias (a partir del 24/03)');
  %xlim([0,DatosEpi(end,1)])
  xlim([data.xdata(1),data.xdata(end)])
  set(gca,'Xtick',diasNum,'XTickLabel',diasEtiq)
  xtickangle(90)
  hold off
end

% % save data results in a file
% datares.chain = chain;
% datares.s2chain = s2chain;
% datares.results = results;
% save ../modelo_A1/datos2020-09-30/ModeloColima-VdeA/resultadosv1/dataresults.mat datares

%%
% PROYECCIONES A LARGO PLAZO
%
%timeorig=data.xdata;
%time=linspace(data.xdata(1),data.xdata(end),100)'; % create time vector (more detail than in data)
time=linspace(timeorig(1),240,500)'; % create time vector (more detail than in data)
y0=data.y0;
y0(20)=theta(11);
y0(22)=theta(11);
%[t,y] = ode45(@modelSystem,time,y0,[],theta); % numerical solution of model with such parameters
ymod=modelSolution(time,theta,y0);
% figure
% % Infecciosos
% subplot(2,3,1)
% plot(t,y(:,3),'m') % Iac
% hold on
% plot(t,y(:,4),'r') % Ihc
% ylabel('Prevalencia')
% title('Infecciosos')
% legend('amb','hosp','Location','NW')
% % Aislados
% subplot(2,3,2)
% plot(t,y(:,5),'m') % Ac
% hold on
% plot(t,y(:,6),'r') % Hc
% ylabel('Prevalencia')
% title('Aislados')
% legend('amb','hosp','Location','NW')
% % Recuperados
% subplot(2,3,3)
% plot(t,y(:,7),'k') % Rc
% hold on
% ylabel('Acumulados')
% title('Recuperados')
% plot(timeorig,data.ydata(:,1),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',2);
% % Defunciones
% subplot(2,3,4)
% plot(t,y(:,8),'k') % Dc
% hold on
% ylabel('Acumulados')
% title('Defunciones')
% plot(timeorig,data.ydata(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',2);
% % Contagios acumulados confirmados
% subplot(2,3,5)
% plot(t,y(:,9),'k') % Cc
% hold on
% ylabel('Acumulados')
% title('Contagios')
% plot(timeorig,data.ydata(:,3),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',2);
% % Contagios acumulados totales
% subplot(2,3,6)
% plot(t,y(:,9)+y(:,12),'b') % Cc+C
% hold on
% ylabel('Acumulados')
% title('Contagios totales')
% 
% for i=1:6 
%     subplot(2,3,i)
%     xlim([0,time(end)])
%     set(gca,'Xtick',diasNum,'XTickLabel',diasEtiq)
%     xtickangle(90)
%     grid on
%     hold off
% end

% subplot(2,3,6)
% plot(t,y(:,9),'--b') % plot numerical solution of model
% hold on
% plot(t,y(:,10),'b')
% plot(t,y(:,9)+y(:,10),'k')
% xlabel('dias (a partir del 24/03)');
% ylabel('Acumulados')
% title('Contagios totales')
% hold off

format long
%defunciones, proporcion de infectados, IFR
%disp([ymod(end,18),(ymod(end,19)+ymod(end,22))/N,100*(ymod(end,18)/(ymod(end,19)+ymod(end,22)))])

%% Predictive plots

% We need to augment the time in ydata to get more plot points
% |datamerge| in the toolbox does this.
data.xdata = datamerge(timeorig,time);

% In order to use the |mcmcpred| function we need
% function |modelfun| with input arguments given as
% |modelfun(xdata,theta)|. We construct this as an anonymous function.
modelfun = @(d,th) modelSolution(d,th,[data.y0(1:19),th(11),data.y0(21),th(11)]);

% We sample 500 parameter realizations from |chain| and |s2chain|
% and calculate the predictive plots.
nsample = 1000;
%%% If s2chain is given to mcmcpred, then the plot shows 95% probability limits 
%%% for new observations and for model parameter uncertainty.
% out = mcmcpred(results,chain,s2chain,data.xdata,modelfun,nsample); %<--NOT SELECTED
%%% If s2chain is not used then the plot contains 50%, 90%, 95%, and 99% predictive 
%%% probability limits due parameter uncertainty.
out = mcmcpred(results,chain,[],data.xdata,modelfun,nsample); % prediction plots
figure
mcmcpredplot(out);
% add the data observations to the plot
hold on
for i=1:7 % depends on number of variables
  subplot(3,3,i)
  hold on
  plot(timeorig,data.ydata(:,i),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',2);
  ylabel(data.ylabels(i)); 
  title(data.titles(i));    
  %xlabel('dias (a partir del 24/03)');
  xlim([0,time(end)])
  set(gca,'Xtick',diasNum,'XTickLabel',diasEtiq)
  xtickangle(90)
  grid on
  hold off
end

%%

function ss = funSS(theta,data)
    % sum-of-squares function
    time   = data.xdata;
    ymodel = modelSolution(time,theta,data.y0);
    ss = sum((ymodel - data.ydata).^2);
end

function ymod=modelSolution(time,theta,y0)
    % model function
    [t,y] = ode45(@modelSystem,time,y0,[],theta);
    ymod=[y(:,3:4),sum(y(:,5:8),2),sum(y(:,9:16),2),y(:,17:19)];
end

function dy = modelSystem(t,y,theta)
    % ode system function 

    beta=theta(1);
    epsilon=1/theta(2);
    pA=theta(3);
    sigma=theta(4);
    alphaA=theta(5);
    alphaH=theta(6);
    gammaAC=theta(7);
    gammaHC=theta(8);
    muC=theta(9);
    gamma=theta(10);
    q=theta(12);

    % S,E,Iac,Ihc,Ac_(1-4),Hc_(1-4),HDc_(1-4),Rc,Dc,Cc,I,R,C
    S=y(1);
    E=y(2);
    Iac=y(3);
    Ihc=y(4);
    
    Ac1=y(5);
    Ac2=y(6);
    Ac3=y(7);
    Ac4=y(8);
    
    Hc1=y(9);
    Hc2=y(10);
    Hc3=y(11);
    Hc4=y(12);

    HDc1=y(13);
    HDc2=y(14);
    HDc3=y(15);
    HDc4=y(16);

    
    Rc=y(17);
    Dc=y(18);
    Cc=y(19);
    I=y(20);
    R=y(21);
    C=y(22);
    
    n=4;
    
    N=S+E+Iac+Ihc+Ac1+Ac2+Ac3+Ac4+Hc1+Hc2+Hc3+Hc4+HDc1+HDc2+HDc3+HDc4+Rc+I+R;
    dS=-beta*S*(Iac+Ihc+I)/N;
    dE=beta*S*(Iac+Ihc+I)/N-sigma*E;
    dIac=pA*epsilon*sigma*E-alphaA*Iac;
    dIhc=(1-pA)*epsilon*sigma*E-alphaH*Ihc;
    
    dAc1=alphaA*Iac-n*gammaAC*Ac1;
    dAc2=n*gammaAC*Ac1-n*gammaAC*Ac2;
    dAc3=n*gammaAC*Ac2-n*gammaAC*Ac3;
    dAc4=n*gammaAC*Ac3-n*gammaAC*Ac4;
    
    dHc1=alphaH*Ihc-n*gammaHC*Hc1;
    dHc2=n*gammaHC*Hc1-n*gammaHC*Hc2;
    dHc3=n*gammaHC*Hc2-n*gammaHC*Hc3;
    dHc4=n*gammaHC*Hc3-n*gammaHC*Hc4;

    dHDc1=q*n*gammaHC*Hc4-n*muC*HDc1;
    dHDc2=n*muC*HDc1-n*muC*HDc2;
    dHDc3=n*muC*HDc2-n*muC*HDc3;
    dHDc4=n*muC*HDc3-n*muC*HDc4;
    
    dRc=n*gammaAC*Ac4+(1-q)*n*gammaHC*Hc4;
    dDc=n*muC*HDc4;
    dCc=epsilon*sigma*E;

    dI=(1-epsilon)*sigma*E-gamma*I;
    dR=gamma*I;
    dC=(1-epsilon)*sigma*E;

    dy = [dS;dE;dIac;dIhc;...
          dAc1;dAc2;dAc3;dAc4;...
          dHc1;dHc2;dHc3;dHc4;...
          dHDc1;dHDc2;dHDc3;dHDc4;...
          dRc;dDc;dCc;dI;dR;dC];

end