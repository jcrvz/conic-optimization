function [X1, Y1, X2, Y2] = plotParabolaXY(Parab, Po, nVar)
    % Ploting Oblique parabola XY from some points Po
    if isempty(Parab), return; end
    
    Prn = true;
    
    A = Parab(1); B = Parab(2); C = Parab(3); 
    D = Parab(4); E = Parab(5); F = Parab(6);

    plot(Po(:,1), Po(:,2), '.r','LineWidth', 1, 'MarkerSize', 1), 
    
    if (nVar < 0), nVar = 1; Prn = false; end
    if (Prn), hold on, end % If negative, we get only the vector points

    if (abs(A) > abs(C))
        Y = -nVar*std(Po(:,2)) + min(Po(:,2)): 0.01: ...
            +nVar*std(Po(:,2)) + max(Po(:,2));
        X1 = (-(D+B*Y)+sqrt((D+B*Y).^2-4*A*(C*Y.^2+E*Y+F)))/(2*A);
        X2 = (-(D+B*Y)-sqrt((D+B*Y).^2-4*A*(C*Y.^2+E*Y+F)))/(2*A);  
        % To avoid complex numbers
        Idy = find (X1 == real(X1));
        Y1 = Y(Idy); Y2= Y1; X1 = X1(Idy); X2 = X2(Idy); 
        if (Prn)
            plot(X1, Y1,'.b', X2, Y2,'.b','LineWidth', 1, 'MarkerSize', 1); 
            hold off,
        end
    else
        X = -nVar*std(Po(:,1)) + min(Po(:,1)): 0.01: ...
            +nVar*std(Po(:,1)) + max(Po(:,1));
        Y1 = (-(E+B*X)+sqrt((E+B*X).^2-4*C*(A*X.^2+D*X+F)))/(2*C);
        Y2 = (-(E+B*X)-sqrt((E+B*X).^2-4*C*(A*X.^2+D*X+F)))/(2*C);
        % To avoid complex numbers
        Idx = find (Y1 == real(Y1));
        X1 = X(Idx); X2 = X1; Y1 = Y1(Idx); Y2 = Y2(Idx); 
        if (Prn)
            plot(X1, Y1,'.b', X2, Y2,'.b','LineWidth', 1, 'MarkerSize', 1); 
            hold off,
        end
    end    
end

% function [X, Y] = plotParabola(Parab, Po)
%     % Ploting parabola rotated from some points Po
%     A1 = Parab(1);  D1 = Parab(2); F1 = Parab(3); 
% 
%     plot(Po(:,1), Po(:,2), 'r*'), hold on
% 
%     X = -2*var(Po(:,1))+min(Po(:,1)): 0.01: ...
%         +2*var(Po(:,1))+max(Po(:,1));
%     Y = (A1*X.^2 + D1*X + F1);
% 
%     plot(X, Y,'b'), hold off
% end