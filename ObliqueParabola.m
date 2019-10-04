function [Parab1, Parab2] = ObliqueParabola(P, Gk, PlotOn)
    % Oblique Parabola Equation obtained from 4 points P
    % PlotON > 2 Numero de figure for plotting
    % Using only 4 points
    Parab1 = [];
    Parab2 = [];

    P2 = P(2,:) - P(1,:);
    P3 = P(3,:) - P(1,:);
    P4 = P(4,:) - P(1,:);
    % It is Posible to use an additional P5 with some error
    P5 = P2;

    % Extracting coordinates
    x2 = P2(1); x3 = P3(1); x4 = P4(1); x5 = P5(1);
    y2 = P2(2); y3 = P3(2); y4 = P4(2); y5 = P5(2);

    % Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0;
    % Considering a refence point P11 = P1-P1 = (0,0) Origing F should be 0
    
    % Ax^2 + Bxy + C*y^2 + Dx + Ey = 0; %Using 2 set of equation (P21,P32)and 
    % (P41, P51) considering P5 = P2 egaling E from Krammer%
    
    D45 = x5*y4 - x4*y5;
    D23 = x2*y3 - x3*y2;
    
    %% AH(1) + B*H(2) + C*H(3) = 0 + B^2-4Ac + Delta ;
    if ((D45==0) || (D23 == 0)) return ; end; % Not solution at all
        H = [(x2 * x3 * (x2-x3)) -(x4 *x5 * (x5-x4));
             (x2 * x3 * (y2-y3)) -(x4 *x5 * (y5-y4));
             (x3*y2^2 - x2*y3^2) -(x4*y5^2 - x5*y4^2)]*[D45; D23];
    
    if ((H(1)==0)&&(H(3)==0))
       % Possible a circle
        Parab1 = GetParab_ABCDEF(1, 0, 1, P(1,:), P2, P3);
        Parab2 = GetParab_ABCDEF(0.01+rand(1), 0, 0.01+rand(1), P(1,:), P2, P3);
        %disp('Circle 1: Proba');
    % H1*A + H2*B + H3*C=0 & Using Condition of Conics B^2-4AC = Gk
    % For stability of conics (Designs)
    else
        if (abs(H(1)) > abs(H(3))) % C = 1 &   (B^2/4)*H1 + H2*B + (H3 + Gk*H1/4)=0
            C  = 1;  H(3) = H(3) - Gk*H(1)/4; % disp('Correction to conics diversity 1');
            DC = H(2)^2 - H(1)*H(3);
            if (DC < 0) return; end
            B1 = 2*(-H(2) + sqrt(DC))/H(1); %1st solution
            B2 = 2*(-H(2) - sqrt(DC))/H(1); %2nd Solution          
            if (B1 == 0) A = 1; else A = (B1^2 - Gk)/(4*C); end %Circle
            Parab1 = GetParab_ABCDEF(A, B1, C, P(1,:), P2, P3);
            if (B2 == 0) A = 1; else A = (B2^2 - Gk)/(4*C); end %Circle
            Parab2 = GetParab_ABCDEF(A, B2, C, P(1,:), P2, P3);
        else                         % A = 1 &   (H1 + Gk*H3/4) + H2*B + (B^2/4)*H3=0
            A  = 1;  H(1) = H(1) - Gk*H(3)/4; % disp('Correction to conics diversity 2'); 
            DC = H(2)^2 - H(1)*H(3);
            if (DC < 0) return; end
            B1 = 2*(-H(2) + sqrt(DC))/H(3); %1st solution
            B2 = 2*(-H(2) - sqrt(DC))/H(3); %2nd Solution
            if (B1 == 0) C = 1; else C = (B1^2 - Gk)/(4*A); end %Circle
            Parab1 = GetParab_ABCDEF(A, B1, C, P(1,:), P2, P3);
            if (B2 == 0) C = 1; else C = (B2^2 - Gk)/(4*A); end %Circle
            Parab2 = GetParab_ABCDEF(A, B2, C, P(1,:), P2, P3);       
        end
    end

    % Ploting results
    if (PlotOn > 0)
        figure(PlotOn)
        subplot(2,2,1), plotParabolaXY(Parab1, P, 5);
        subplot(2,2,2), plotParabolaXY(Parab2, P, 5);
        % Rotating parabolas
        R =@(th) [cos(th) -sin(th); sin(th) cos(th)];
        [ParRot1, th1] = RotateParabola(Parab1);
        [ParRot2, th2] = RotateParabola(Parab2);
        subplot(2,2,3), plotParabolaXY(ParRot1, P*R(th1), 5);
        subplot(2,2,4), plotParabolaXY(ParRot2, P*R(th2), 5);
    end
end

function [ObliquePar] = GetParab_ABCDEF(A, B, C, P1, P2, P3)
   % Translating to original coordinates (to P1)
    h = P1(1);  k = P1(2); %Original Reference point
    x2= P2(1);  y2= P2(2); % Traslated second point
    x3= P3(1);  y3= P3(2); % Traslated third  point
    
    F = 0;
    H = [x2  x3; -y2  -y3]*...
        [(A*x3^2+B*x3*y3+C*y3^2); -(A*x2^2+B*x2*y2+C*y2^2)]/...
        (x3*y2-x2*y3); 
    D = H(2); E = H(1);     
    %[A B C D E F] Parabola passing by origin
    
    % Translating to original coordinate P
    F = A*h^2 + B*h*k + C*k^2 - D*h - E*k;
    D = D - 2*A*h - B*k;
    E = E - 2*C*k - B*h;
    ObliquePar = [A B C D E F]; %Parabola passing by the original points
end

function [ParRotated,th]=RotateParabola(Parab)
    %Rotating  x^2 + B*x*y + C*y^2 + D*x +E*y + F;  C=B^/4
    %Output    y = A1 *x^2 + D1*x + F1  (Rotation angle)
    ParRotated = [];
    A = Parab(1); B = Parab(2); C = Parab(3); 
    D = Parab(4); E = Parab(5); F = Parab(6);

    % Rotation angle
    th = atan2(B, A-C )/2;

    A2 = A*(cos(th))^2 + B*sin(th)*cos(th) + C*(sin(th))^2;
    B2 = (B/2)*cos(2*th)-(A-C)*sin(2*th)/2;
    C2 = A*(sin(th))^2 - B*sin(th)*cos(th) + C*(cos(th))^2;
    D2 = (+D*cos(th) + E*sin(th));
    E2 = (-D*sin(th) + E*cos(th));
    F2 = F; 
    ParRotated = [A2 B2 C2 D2 E2 F2];
end

