function PlotFit(models,S_contracts,VIX_contracts,varargin)
%
%    Plots the fit of a set of models to a set of market prices for options on the underlying
%    asset, written S, and the VIX index.
%
% -----------------------------------------------------------------------------------------------
%  Main parameters  
% -----------------------------------------------------------------------------------------------
%   models:          [1xN cell]       Set of pricing models. Must be of type PricingModelClass.
%
%   S_contracts:     [1 x 1 struct]   Information on S(t) options. The following fields apply:
%       o k:         [Ns x 1 real]    Log-moneyness values.
%       o ttm:       [Ns x 1 real]    Maturities.
%       o bid:       [Ns x 1 real]    Implied volatilities (bid). 
%       o ask:       [Ns x 1 real]    Implied volatilities (ask). 
%
%   VIX_contracts:   [1 x 1 struct]   Information on VIX options. The following fields apply:
%       o K:         |Nv x 1 real]    Strikes.
%       o ttm:       [Nv x 1 real]    Maturities.
%       o bid:       [Nv x 1 real]    Implied volatilities (bid). 
%       o ask:       [Nv x 1 real]    Implied volatilities (ask). 
%
%   Either of the parameters 'S_contracts' or 'VIX_contracts' can be left empty.
%
% -----------------------------------------------------------------------------------------------
%  Optional parameters  
% -----------------------------------------------------------------------------------------------
%  * Specify as name-value pairs
%
%   modelNames:      [1xN cell]     Model names. 
%   modelColors:     [Nx3 real]     RGB colors for each model. 
%   marketColors:    [2x3 real]     RGB colors for observed prices. Order is: bid, ask.
%   arrangePlots:    [1x1 logical]  Arrange plot windows for easier viewing (default = false). 
%
% -----------------------------------------------------------------------------------------------

    % Default color scheme:
    nModels = size(models,1);
    cm = hsv(nModels + 2);
    
    % Default model names:
    modelNames = cell(nModels,1);
    for i=1:size(nModels,1);modelNames{i}=num2str(i);end

    % Parse variable inputs:
    p = inputParser;
    addParameter(p,'modelNames',modelNames);
    addParameter(p,'modelColors',cm(1:nModels,:));
    addParameter(p,'marketColors',cm(nModels+1:end,:));
    addParameter(p,'arrangePlots',false);
    parse(p,varargin{:});
    v2struct(p.Results);
    
    if ~isempty(S_contracts) && ~isStructFieldsSameSize(S_contracts)
        error('PlotFit: ''S_contracts'' fields must be of the same size.');
    end
    if ~isempty(VIX_contracts) && ~isStructFieldsSameSize(VIX_contracts)
        error('PlotFit: ''VIX_contracts'' fields must be of the same size.');
    end
    
    if ~isempty(S_contracts)
        S_contracts_price = struct;
        S_contracts_price.k = S_contracts.k;
        S_contracts_price.ttm = S_contracts.ttm;
    else
        S_contracts_price = [];
    end
    if ~isempty(VIX_contracts)
        VIX_contracts_price = struct;
        VIX_contracts_price.K = VIX_contracts.K;
        VIX_contracts_price.ttm = VIX_contracts.ttm;
    else
        VIX_contracts_price = [];
    end
    
    % Compute model prices:
    if ~isempty(S_contracts);S_contracts.priceFormat = 'implied_volatility';end
    if ~isempty(VIX_contracts);VIX_contracts.priceFormat = 'implied_volatility';end
    [S_modelPrices,VIX_modelPrices] = deal(cell(nModels,1));
    for i=1:size(models,1)
        [S_modelPrices{i},VIX_modelPrices{i}] = models{i}.GetPrices(S_contracts_price,...
                                                                    VIX_contracts_price);
    end
    
    % Plot options on S:
    if ~isempty(S_contracts)
        uniqT = unique(S_contracts.ttm);
        for i=1:size(uniqT,1)
            idx = S_contracts.ttm == uniqT(i);
            if any(idx)
                % Plot market prices:
                figure;
                plot(S_contracts.k(idx),S_contracts.bid(idx),'.','DisplayName','Bid',...
                     'Color',marketColors(1,:));hold on;
                plot(S_contracts.k(idx),S_contracts.ask(idx),'.','DisplayName','Ask',...
                     'Color',marketColors(2,:));hold on;
                title(['T = ',num2str(uniqT(i))]);
                xlabel('Log-moneyness');
                ylabel('Implied volatility');

                % Plot model prices:
                for j=1:size(nModels,1)
                    plot(S_modelPrices{j}.k(idx),S_modelPrices{j}.p(idx),'-',...
                         'DisplayName',modelNames{j},'Color',modelColors(j,:));hold on;            
                end
                
                if nModels == 1
                    % Plot confidence intervals too:
                    plot(S_modelPrices{1}.k(idx),...
                        S_modelPrices{1}.p(idx) - 1.96*S_modelPrices{1}.se(idx),'--',...
                         'DisplayName','95 % CI','Color',modelColors(1,:));hold on;   
                    plot(S_modelPrices{1}.k(idx),...
                        S_modelPrices{1}.p(idx) + 1.96*S_modelPrices{1}.se(idx),'--',...
                         'DisplayName','95 % CI','Color',modelColors(1,:));hold on;                        
                end
                
                hold off;

                legend();

            end
        end
    end

    % Plot options on VIX: 
    if ~isempty(VIX_contracts)
        uniqT = unique(VIX_contracts.ttm);
        for i=1:size(uniqT,1)
            idx = VIX_contracts.ttm == uniqT(i);
            if any(idx)
                % Plot market prices:
                figure;
                plot(VIX_contracts.K(idx),VIX_contracts.bid(idx),'.','DisplayName','Bid',...
                     'Color',marketColors(1,:));hold on;
                plot(VIX_contracts.K(idx),VIX_contracts.ask(idx),'.','DisplayName','Ask',...
                     'Color',marketColors(2,:));hold on;
                title(['T = ',num2str(uniqT(i))]);
                xlabel('Strike');
                ylabel('Implied volatility');
                
                % Plot model prices:
                for j=1:size(nModels,1)
                    plot(VIX_modelPrices{j}.K(idx),VIX_modelPrices{j}.p(idx),'-',...
                         'DisplayName',modelNames{j},'Color',modelColors(j,:));hold on;            
                end
                
                if nModels == 1
                    % Plot confidence intervals too:
                    plot(VIX_modelPrices{1}.K(idx),...
                        VIX_modelPrices{1}.p(idx) - 1.96*VIX_modelPrices{1}.se(idx),'--',...
                         'DisplayName','95 % CI','Color',modelColors(1,:));hold on;   
                    plot(VIX_modelPrices{1}.K(idx),...
                        VIX_modelPrices{1}.p(idx) + 1.96*VIX_modelPrices{1}.se(idx),'--',...
                         'DisplayName','95 % CI','Color',modelColors(1,:));hold on;                        
                end
                
                hold off;

                legend();
            end
        end
    end
    
    % Rearrange plot windows:
    if arrangePlots
        autoArrangeFigures();
    end

end

