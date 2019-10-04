function [MAE, RMSE, MP, X1, Y1, X2, Y2] = ParabolaErrors(Parab, Px, Py, th)
    % testing fitting on parabola PxY from some points Po
        infty= 1e30;
        MAE  = infty;
        RMSE = infty; 
        MP   = 1;
        X1   = [];  X2 = [];
        Y1   = [];  Y2 = [];
      
        if isempty(Parab), return; end
                       
        A = Parab(1); B = Parab(2); C = Parab(3); 
        D = Parab(4); E = Parab(5); F = Parab(6);   

        %% Root Mean Square Error RMSE & MAE Mean Absolute Error

           if A ~= 0
               X1 = (-(D+B*Py)+sqrt((D+B*Py).^2-4*A*(C*Py.^2+E*Py+F)))/(2*A);
               X2 = (-(D+B*Py)-sqrt((D+B*Py).^2-4*A*(C*Py.^2+E*Py+F)))/(2*A);
               I1 = find (X1 == real(X1));
               Y1 = Py(I1); Y2 = Y1;  X1 = X1(I1); X2 = X2(I1); 
               Err1 = min(abs(Px(I1)-X1), abs(Px(I1)-X2));
               Rn1 = max(Px) - min(Px);  
               N1   = length(I1)+1;
           else
               N1 = 1; 
               Rn1 = 1; 
               Err1 = infty; 
           end
          if C ~= 0
              Y1 = (-(E+B*Px)+sqrt((E+B*Px).^2-4*C*(A*Px.^2+D*Px+F)))/(2*C);
              Y2 = (-(E+B*Px)-sqrt((E+B*Px).^2-4*C*(A*Px.^2+D*Px+F)))/(2*C);
              I2 = find (Y1 == real(Y1));
              X1 = Px(I2); X2 = X1; Y1 = Y1(I2); Y2 = Y2(I2);
              Err2 = min(abs(Py(I2)-Y1), abs(Py(I2)-Y2));
              Rn2 = max(Py) - min(Py);
              N2   = length(I2)+1;
          else
              N2 = 1; 
              Rn2 = 1; 
              Err2 = infty; 
          end
          
        chk1 = (N1~=1); chk2 = (N2~=1);
        MAE  = (chk1)*sum(Err1)/N1/Rn1 + (chk2)*sum(Err2)/N2/Rn2;
        RMSE = (chk1)*sqrt(sum(Err1.^2)/N1)/Rn1 +...
               (chk2)*sqrt(sum(Err2.^2)/N2)/Rn2;
        %  chi-2;
        % Normalized: Matching points threshold (in pixels)
        MP   = 1 + sum(Err1 < th)+sum(Err2 < th);  % th < 10; %INLIERS
        % MAE  = 100*MAE /MP; % MAE
        % RMSE = 100*RMSE/MP; % RMSE
end

