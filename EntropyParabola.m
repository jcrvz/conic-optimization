function [Entropy] = EntropyParabola(Parab, Px, Py, th, P)
%%% For this implementation entropy is not measured.
%%% MSE is used instedad
         tiny = 1e-30;
         Entropy = 0;
         % Coef B
         if isempty(Parab) return; end;
        % if ~isreal(Parab(2))  return; end  

    %% Root Mean Square Error RMSE & MAE Mean Absolute Error
    %Translaped Binary histogram between image and proposed parabola
         [MAE, RMSE, MP, X1, Y1, X2, Y2] = ParabolaErrors(Parab, Px, Py, th);
         if (RMSE ==0 || isnan(RMSE)) Entropy = tiny;
         else Entropy = MP/RMSE;  end 
         
         % (1) Firts  option MP or MP^2/RMSE
         
%          Prob = MP/(length(Px)^2); % To tunning for application
%          if (Prob > 0)  % (2) Second option
%               Entropy = -((1-Prob)*log2(1-Prob) + Prob*log2(Prob));
%          end
         
         % In the image scope for videos
         if (P>0)
                     figure(10)
                     plot (Px, Py, '.r'), hold on
                     % set(par,'AlphaData',0.75)
                     plot(X1, Y1, '.b', X2, Y2 ,'.b');
                     plot(P(:,1), P(:,2), '*g');
                     axis([min(Px) max(Px) min(Py) max(Py)]);
                     title(sprintf('Entropy = %2.8f   %2.8f', Entropy, MP));
                     hold off
                     pause(0.20);
         end
end