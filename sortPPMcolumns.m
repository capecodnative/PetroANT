% Sort *_ppm columns to the end of table T
% Ordering is based on completeness and abundance, with flexible user options.
%
% OUTPUT:
%   Tsorted          reordered table
%   ppmStats         table of stats for each *_ppm column
%
% USER OPTIONS:
%   opts.completenessMode : 'count' or 'fraction'
%   opts.abundanceMode    : 'median' or 'logmedian'
%   opts.sortDirection    : 'descend' or 'ascend'
%   opts.weightComplete   : weight on completeness metric
%   opts.weightAbundance  : weight on abundance metric
%   opts.normalizeMetrics : true/false
%   opts.zeroAsMissing    : true/false
%
% NOTES:
% - All non-*_ppm columns stay in their original order at the front.
% - All *_ppm columns are moved to the end.
% - By default, columns with more populated values and higher median abundance
%   sort earlier within the *_ppm block.

%% ---------------- USER OPTIONS ----------------
opts.completenessMode = 'count';     % 'count' or 'fraction'
opts.abundanceMode    = 'median';    % 'median' or 'logmedian'
opts.sortDirection    = 'descend';   % 'descend' or 'ascend'
opts.weightComplete   = 0.9;         % weight on completeness metric
opts.weightAbundance  = 0.1;         % weight on abundance metric
opts.normalizeMetrics = true;        % normalize metrics to 0-1 before combining
opts.zeroAsMissing    = false;       % if true, treat 0 like missing
%% ------------------------------------------------

varNames = PetroANT_DCO.Properties.VariableNames;
isMol = endsWith(varNames, '_mol');

molVars    = varNames(isMol);
otherVars  = varNames(~isMol);

nMol = numel(molVars);

nonNanCount = nan(nMol,1);
nonNanFrac  = nan(nMol,1);
medVal      = nan(nMol,1);

for i = 1:nMol
    x = PetroANT_DCO.(molVars{i});

    % Skip nonnumeric columns safely
    if ~isnumeric(x)
        continue
    end

    x = x(:);

    if opts.zeroAsMissing
        valid = ~isnan(x) & x ~= 0;
    else
        valid = ~isnan(x);
    end

    xv = x(valid);

    nonNanCount(i) = sum(valid);
    nonNanFrac(i)  = sum(valid) / numel(x);

    if ~isempty(xv)
        switch lower(opts.abundanceMode)
            case 'median'
                medVal(i) = median(xv, 'omitnan');
            case 'logmedian'
                % Only positive values allowed for log metric
                xv = xv(xv > 0);
                if ~isempty(xv)
                    medVal(i) = median(log10(xv), 'omitnan');
                end
            otherwise
                error('Unknown opts.abundanceMode: %s', opts.abundanceMode);
        end
    end
end

% Choose completeness metric
switch lower(opts.completenessMode)
    case 'count'
        completeMetric = nonNanCount;
    case 'fraction'
        completeMetric = nonNanFrac;
    otherwise
        error('Unknown opts.completenessMode: %s', opts.completenessMode);
end

abundanceMetric = medVal;

% Optional normalization to comparable scales
if opts.normalizeMetrics
    cMetric = normalizeToUnitRange(completeMetric);
    aMetric = normalizeToUnitRange(abundanceMetric);
else
    cMetric = completeMetric;
    aMetric = abundanceMetric;
end

% Combined sort score
sortScore = opts.weightComplete .* cMetric + opts.weightAbundance .* aMetric;

% Build stats table for inspection
ppmStats = table( ...
    molVars(:), ...
    nonNanCount, ...
    nonNanFrac, ...
    medVal, ...
    cMetric, ...
    aMetric, ...
    sortScore, ...
    'VariableNames', {'Variable','NonNanCount','NonNanFraction','MedianValue', ...
                      'CompleteMetric','AbundanceMetric','SortScore'});

% Sort ppm stats
switch lower(opts.sortDirection)
    case 'descend'
        ppmStats = sortrows(ppmStats, {'SortScore','NonNanCount','MedianValue'}, {'descend','descend','descend'});
    case 'ascend'
        ppmStats = sortrows(ppmStats, {'SortScore','NonNanCount','MedianValue'}, {'ascend','ascend','ascend'});
    otherwise
        error('Unknown opts.sortDirection: %s', opts.sortDirection);
end

% Reorder table columns: non-ppm first (original order), ppm at end (sorted)
newOrder = [otherVars, ppmStats.Variable(:)'];
Tsorted = PetroANT_DCO(:, newOrder);

%% -------- helper function --------
function y = normalizeToUnitRange(x)
    y = nan(size(x));
    valid = ~isnan(x);
    if ~any(valid)
        return
    end

    xmin = min(x(valid));
    xmax = max(x(valid));

    if xmax == xmin
        y(valid) = 1;
    else
        y(valid) = (x(valid) - xmin) ./ (xmax - xmin);
    end
end

clear opts varNames isMol molVars otherVars nMol ...
      nonNanCount nonNanFrac medVal ...
      x valid xv ...
      completeMetric abundanceMetric ...
      cMetric aMetric sortScore ...
      ppmStats newOrder i