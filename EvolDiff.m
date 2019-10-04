function [BParab1, NIter] = EvolDiff(Px, Py, th, Gk) 
% Differential Evolution 
% gavina - Implementation
% avina@ugto.mx

nPix    = length(Px);
sF      = 0.125;        % Factor to speed up, reducing useful points [0,1]
D       = 6;            % Number of Population to optimize
G       = 250;          % Number of generations or generationSize
N       = round(sF*nPix);    % Population Size (Portion of Total Population)
nPE     = 4;            % nPoints to obtain parabola equation
F       = 0.5;             % [0, 1]   %0.5
CR      = 0.8;            % [0.1, 0.3] % 0.8
F2      = F/3;  % F/2;
Error   = 1e-3; % 5e-2

    Population  = zeros(N, D);
    Donors      = zeros(N, D);
    Trials      = zeros(N, D);
    cc          = zeros(N, 1);
    BParab1     = zeros(1, D);
    BEntropy    = 0;
    IndCount    = 1;
    fsort       = ones(N, 1);

    Stats = zeros(G,1);
    
rng('shuffle');

%% Initialization
    Set4P  = randi(nPix, N, nPE);
    for i = 1 : N % Calcule of Population for population        
        Pij= Set4P(i,:);           
        [Parab1, Parab2] = ObliqueParabola([Px(Pij) Py(Pij)], Gk, 0);
        while isempty(Parab1)%% Choising only valid parabolas
              Pij= randi(nPix, 1, nPE);
             [Parab1, Parab2] = ObliqueParabola([Px(Pij) Py(Pij)], Gk, 0);
        end
          
        Entropy1 = EntropyParabola(Parab1, Px, Py, th, 0);
        Entropy2 = EntropyParabola(Parab2, Px, Py, th, 0);
           
        if (Entropy1 > Entropy2) Population(i,:) = Parab1;            
        else                     Population(i,:) = Parab2; end
    end        
    P = [Px(Pij) Py(Pij)];
    
for gen = 1 : G  

    %% Full mutation
    rij = randi(N, N, 3);
    for i = 1 : N
        ri = rij(i,:); 
        %% Searching repetition at least once
        if ((ri(1)==ri(2))||(ri(1)==ri(3))||(ri(3)==ri(2))|| ...
            (ri(1)== i   )||(ri(2)== i   )||(ri(3)== i  )) ...
            rij(i,:)= randi(N, 1, 3); end
    end
        
    Donors  =    Population(rij(:,1),:) + ...
              F *(Population(rij(:,2),:) - Population(rij(:,3),:)) +...
              F2*(repmat(BParab1,N,1)    - Population) ;

    %% Crossover or Recombinations 

    ijRand = rand (N, D);
     Jrand = randi(D, N, D);
        
     %if ((ijRand(i,j)<=CR)||(j==Jrand)) Trials(i,j) = Donors(i,j);     end  
     %if ((ijRand(i,j)> CR)&&(j~=Jrand)) Trials(i,j)=Population(i,j); end
           
     check = (ijRand <= CR)|(repmat(1:D, N, 1)==Jrand);
     Trials = (check).*Donors + (~check).*Population;
     
    %% Selection

    for i = 1 : N
        Entropy1 = EntropyParabola(    Trials(i,:), Px, Py, th, 0);
        Entropy2 = EntropyParabola(Population(i,:), Px, Py, th, 0);
       
        if (Entropy1 > Entropy2) 
             cc(i) = Entropy1; Population(i,:) = Trials(i,:); 
        else, cc(i) = Entropy2; end
    end
 
    % Normalize first coefficients
%     check = abs(Population(:,1)) > abs(Population(:,3));
%     Cnorm = repmat(check.*Population(:,1) + ~check.*Population(:,3), 1, D);
%     Population = Population./Cnorm;
    
    %Population = Population./repmat(max(abs(Population'))',1,D);
    
    %% Stopping conditions    
    [maxCC, maxId1] = max(cc);        fsort(maxId1) =  0; %Detecting 1 Max Entropy
    [~,     maxId2] = max(fsort.*cc); fsort(maxId1) =  1; %Detecting 2 Max Entropy   
%     
    Parab1 = Population(maxId1,:);
    Parab2 = Population(maxId2,:); 
    %  Div = (abs(Parab1))/1; Div(find(Div == 0)) = 0.001;
    Div = (0.80*abs(Parab1)+0.2*abs(Parab2));  Div(find(Div == 0)) = eps;
    
    NormErr= (Parab1-Parab2)./Div; %mean(Population)
    CurrentErr = sum(NormErr.^2);

    Stats(gen) = (CurrentErr);
    if (CurrentErr<Error && gen > 25)
        fprintf('Gen: %d, Error: %.4e\n', gen, CurrentErr);
        break; 
    end %2.5
       
    if (BEntropy < maxCC) 
        BEntropy = EntropyParabola(Parab1, Px, Py, th, 0);%P
        BParab1  = Parab1;
    end
end
     NIter  = gen;
end

function [Average, StdDev]=plot2Ev(fname, Title, XLabel, YLabel, Data, NoFig)
Average = mean(Data);
StdDev  = std(Data);

set(0,'defaultLineMarkerSize',  1); % set the default line marker size to msz
set(0,'defaulttextInterpreter','latex')
set(0,'defaultLineLineWidth',   0.75);   % set the default line width to lw
 
width  = 2.5;    % Width in inches
height = 3/4*width;    % Height in inches

set(gcf, 'Units','inches', 'PaperUnits', 'inches',...
         'Position', [0 0 width height], 'PaperPositionMode', 'auto'); 
     
axis tight
 
set(gca,'TickLabelInterpreter','LaTeX','LooseInset', get(gca,'TightInset'),...
        'FontUnits', 'points', 'FontWeight', 'normal','FontSize', 14, ...
        'FontName', 'Times', 'LineWidth', 1.25)

ax = gca;
%  
 outerpos = ax.OuterPosition;
 ti = ax.TightInset; 
 left = outerpos(1) + ti(1);
 bottom = outerpos(2) + ti(2);
 ax_width = outerpos(3) - ti(1) - ti(3);
 ax_height = outerpos(4) - ti(2) - ti(4);
 ax.Position = [left bottom ax_width ax_height];


if length(Data) > 10
    ax.XLim  = [0 length(Data)];
    ax.XTick = [0:length(Data)/10:length(Data)]; 
end

MinD = 0; %min(Data)
MaxD = max(Data);
ax.YLim  = [MinD MaxD];
ax.YTick = [MinD : (MaxD-MinD)/5:MaxD];

fig = figure(NoFig); 
hold on, plot(Data,'b'), plot(Data, '*r'), grid, hold off

title(Title);
ylabel(XLabel);
xlabel(YLabel);

%legend({'MSE Error'}, 'FontUnits','points',...
%       'Interpreter','latex','FontName','Times', 'Location', 'SouthEast');

% print(gcf, fname, '-dpdf' );
% system(sprintf('pdfcrop --margins 1 %s %s', fname, fname))
%print(gcf, 'Result01.eps','-opengl','-depsc', '-tiff','-r300');
end


