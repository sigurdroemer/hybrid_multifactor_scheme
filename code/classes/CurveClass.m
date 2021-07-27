classdef CurveClass < GenericSuperClass
%
% 	Defines a general curve object. A curve is here defined by specifying a set of fixed points 
%   (t_i,y_i), i=1,...,N, and then placing an interpolation and extrapolation method on top.
%
% -----------------------------------------------------------------------------------------------
%  Properties
% -----------------------------------------------------------------------------------------------
%    gridpoints:    [Nx1 real]   Grid points t_i, i = 1,...,N. Must be sorted in ascending order.
%    values:        [Nx1 real]   Curve values y_i, i = 1,...,N.
%    interpolation: [1x1 string] Interpolation method, options are 'flat', 'linear', 'pchip'.
%    extrapolation: [1x1 string] Extrapolation method, currently the only option is 'flat.
%
% -----------------------------------------------------------------------------------------------

properties
    gridpoints    
    values        
    interpolation 
    extrapolation 
end
    
methods
function obj = CurveClass(varargin)
%
%   Constructor.
%
% -----------------------------------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------------------------------
%   varargin: Inputs must come in name-value pairs corresponding to each property of the object. 
%             The properties of the object will then be set to the values. Default values may be 
%             chosen if some properties are not set.
%
% -----------------------------------------------------------------------------------------------
% Output
% -----------------------------------------------------------------------------------------------
%   [1x1 CurveClass] The object.
%
% -----------------------------------------------------------------------------------------------
%

    % Parse input name-value pairs to object properties:
    obj.ParseConstructorInputs(varargin{:});

    % Validation checks:
    if ~isequal(size(obj.values),size(obj.gridpoints)) ...
            || size(obj.values,2) > 1 || size(obj.gridpoints,2) > 1
        error(['CurveClass: Grid points and values',...
               ' must be column vectors of the same size.']);
    end

    % Default values:
    if isempty(obj.interpolation)
        obj.interpolation = 'flat';
    else
        if ~(strcmpi(obj.interpolation,'flat') || ...
             strcmpi(obj.interpolation,'linear') || strcmpi(obj.interpolation,'pchip'))
          error('CurveClass: Invalid choice of interpolation method.');
        end
    end
    if isempty(obj.extrapolation)
        obj.extrapolation = 'flat';
    else
        if ~strcmpi(obj.extrapolation,'flat')
          error('CurveClass: Invalid choice of extrapolation method.');
        end        
    end

end
function vals = Eval(obj,pts)
%
%   Evaluates the curve using interpolation and extrapolation as necessary.
%
% ---------------------------------------------------------------------------------------------
% Parameter
% ---------------------------------------------------------------------------------------------
%   pts: [NxM real] Evaluation points.
%
% ---------------------------------------------------------------------------------------------
% Output
% ---------------------------------------------------------------------------------------------
%   vals: [NxM real] Value of curve at evaluation points.
%
% ---------------------------------------------------------------------------------------------
%

    N = size(pts,1);M = size(pts,2);
    
    if M > 1
        % Compute columns recursively:
        vals = zeros(N,M);
        for i=1:size(vals,2)
            vals(:,i) = obj.Eval(pts(:,i));
        end
        return;
    end

    if ~strcmpi(obj.extrapolation,'flat')
        error(['CurveClass:GetValues: Extrapolation method', ...
              ' is not supported.']);
    end

    switch obj.interpolation
        case {'linear','pchip'}
            [maxPt, idxMax] = max(obj.gridpoints);
            [minPt, idxMin] = min(obj.gridpoints);
            idxHigh = pts >= maxPt;
            idxLow = pts <= minPt;
            idxMiddle = ~idxHigh & ~idxLow;
            vals = NaN(size(pts));
            if strcmpi(obj.interpolation,'linear')
                vals(idxMiddle) = interp1(obj.gridpoints,obj.values,pts(idxMiddle));
            else
                vals(idxMiddle) = pchip(obj.gridpoints,obj.values,pts(idxMiddle));
            end
            if any(idxHigh)
                vals(idxHigh) = obj.values(idxMax);
            end
            if any(idxLow)
                vals(idxLow) = obj.values(idxMin);
            end                
        case 'flat'
            n = size(obj.gridpoints,1);
            vals = reshape(obj.values(min(...
                            sum(bsxfun(@(a,b)(a < b) ,...
                            obj.gridpoints,pts(:)'), 1) + 1, n)), ...
                            size(pts, 1), size(pts, 2));                            
        otherwise
            error(['CurveClass:GetValues: Interpolation method', ...
                  ' is not supported.']);
    end

end
function vals = Integrate(obj,a,b)
%
%   Returns the integrated curve.
%
% ---------------------------------------------------------------------------------------------
% Parameters 
% ---------------------------------------------------------------------------------------------
%   a: [Nx1 real] Lower end-points.
%   b: [Nx1 real] Upper end-points.
%
% ---------------------------------------------------------------------------------------------
% Output
% ---------------------------------------------------------------------------------------------
%   vals: [Nx1 real] Curve integrated from a to b.
%
% ---------------------------------------------------------------------------------------------
%

    if size(a,1) ~= size(b,1)
        error('CurveClass:Integrate: The input vectors must have the same length');
    end

    if strcmpi(obj.interpolation,'flat') && strcmpi(obj.extrapolation,'flat') ...
            && size(obj.values,1) == 1
        vals = obj.values(1)*(b - a);

    else
        vals = NaN(size(a));
        for i=1:size(a,1)
            vals(i) = integral(@(t)(obj.Eval(t)),a(i),b(i));
        end

    end

end
function plot(obj,curves,a,b,nPts)
%
%   Creates a simple plot showing the curve.
%
% ---------------------------------------------------------------------------------------------
% Parameters 
% ---------------------------------------------------------------------------------------------
%   curves:  [1xN cell]      Each element is a [1x1 CurveClass] that will then be plotted. If
%                            empty we default to only plotting the current object.
%   a:       [1x1 real]      Lower end-point. Default is 0.
%   b:       [1x1 real]      Upper end-point. Default is max(obj.gridpoints).
%   nPts:    [1x1 integer]   Number of points to sample between a and b. Default is 1000.
%
% ---------------------------------------------------------------------------------------------
%    
    
    if ~exist('a','var') || isempty(a)
        a = 0;
    end
    if ~exist('b','var') || isempty(b)
        b = max(obj.gridpoints);
    end
    if ~exist('nPts','var') || isempty(nPts)
        nPts = 1000;
    end    
    
    dx = (b-a)/(nPts-1);
    x = (a:dx:b)';
    
    if ~exist('curves','var') || isempty(curves)
        curves = {obj};
    end    
        
    figure;
    for i=1:size(curves,2)
        y = curves{i}.Eval(x);
        plot(x,y,'DisplayName',['Interpolated curve (',num2str(i),')']);hold on;
    end
    xlabel('t');
    ylabel('Values');
    legend();
    xlim([a,b]);
    hold off;    
    
end
end

end

