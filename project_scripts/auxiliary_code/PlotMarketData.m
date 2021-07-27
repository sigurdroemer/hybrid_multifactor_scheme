function PlotMarketData(S_contracts,VIX_contracts,arrangePlots)
%
%    Plots the fit of a set of models to a set of market prices for options on the underlying
%    asset, written S, and the VIX index. Also plots yield curves.
%
%    VIX options are currently not supported.
%
% -----------------------------------------------------------------------------------------------
%  Main parameters  
% -----------------------------------------------------------------------------------------------
%   S_contracts:     [1 x 1 struct]     Information on S(t) options. The following fields apply:
%       o k:         [Ns x 1 real]      Log-moneyness values.
%       o T:         [Ns x 1 real]      Maturities.
%       o bid_iv:    [Ns x 1 real]      Implied volatilities (bid). 
%       o ask_iv:    [Ns x 1 real]      Implied volatilities (ask). 
%       o mid_iv:    [Ns x 1 real]      Implied volatilities (mid). 
%       o q:         [1x1 CurveClass]   Dividend yield curve.
%       o y:         [1x1 CurveClass]   Interest rate yield curve.
%
%   VIX_contracts:     [1 x 1 struct]     Information on VIX options. The following fields apply:
%       o K:           [Nv x 1 real]      Strikes.
%       o T:           [Nv x 1 real]      Maturities.
%       o bid_iv:      [Nv x 1 real]      Implied volatilities (bid). 
%       o ask_iv:      [Nv x 1 real]      Implied volatilities (ask). 
%       o mid_iv:      [Nv x 1 real]      Implied volatilities (mid). 
%       o F:           [Nv x 1 real]      VIX futures.
%       o vix_index:   [1 x 1 real]       VIX index.
%
%   Either of the parameters 'S_contracts' or 'VIX_contracts' can be left empty.
%
%   arrangePlots:    [1x1 logical]  Arrange plot windows for easier viewing (default = true). 
%
% -----------------------------------------------------------------------------------------------

    % Parse variable inputs:
    if ~exist('arrangePlots','var') || isempty(arrangePlots)
        arrangePlots = true;
    end
    
    if ~exist('VIX_contracts','var') || isempty(VIX_contracts)
        VIX_contracts = [];
    end
    
    % Plot options on S:
    if ~isempty(S_contracts)
        uniqT = unique(S_contracts.T);
        for i=1:size(uniqT,1)
            idx = S_contracts.T == uniqT(i);
            if any(idx)
                % Plot market prices:
                figure;
                plot(S_contracts.k(idx),S_contracts.bid_iv(idx),'.','DisplayName','Bid',...
                     'Color','Blue');hold on;
                plot(S_contracts.k(idx),S_contracts.ask_iv(idx),'.','DisplayName','Ask',...
                     'Color','Red');hold on;
                plot(S_contracts.k(idx),S_contracts.mid_iv(idx),'--','DisplayName','mid',...
                     'Color','Green');hold on;                 
                title(['T = ',num2str(uniqT(i))]);
                xlabel('Log-moneyness');
                ylabel('Implied volatility');

                hold off;

                legend();

            end
        end
        S_contracts.q.plot();title('Dividend yields');
        ylim([min(S_contracts.q.values)-0.01,max(S_contracts.q.values)+0.01])
        S_contracts.y.plot();title('Interest rate yields');
    end
    
    % Plot options on VIX: 
    if ~isempty(VIX_contracts)
        uniqT = unique(VIX_contracts.T);
        for i=1:size(uniqT,1)
            idx = VIX_contracts.T == uniqT(i);
            if any(idx)
                % Plot market prices:
                figure;
                plot(VIX_contracts.k(idx),VIX_contracts.bid_iv(idx),'.','DisplayName','Bid',...
                     'Color','Blue');hold on;
                plot(VIX_contracts.k(idx),VIX_contracts.ask_iv(idx),'.','DisplayName','Ask',...
                     'Color','Red');hold on;
                plot(VIX_contracts.k(idx),VIX_contracts.mid_iv(idx),'--','DisplayName','mid',...
                     'Color','Green');hold on;    
                title(['T = ',num2str(uniqT(i))]);
                xlabel('Log-moneyness');
                ylabel('Implied volatility');
                hold off;
                legend();
            end
        end
        [uniqT,idx,~] = unique(VIX_contracts.T);
        figure;
        plot(uniqT,VIX_contracts.F(idx),'x-');hold on;
        plot(1/12,VIX_contracts.vix_index,'x-');
        legend('VIX futures','VIX index');xlabel('Maturity');ylabel('Volatility');
        
    end
    
    % Rearrange plot windows:
    if arrangePlots
        autoArrangeFigures();
    end

end

