function interfaceMatrix = buildInterfaceMatrix(cellValues, cellScaling, gridSpacing, method)
% buildTransportMatrix Creates a matrix for transport in battery simulations.

% Get values and spacings for neighboring cells
leftValues = cellValues(1:end-1);
rightValues = cellValues(2:end);
leftSpacing = gridSpacing(1:end-1);
rightSpacing = gridSpacing(2:end);

% Calculate property at cell interfaces
switch method
    case "harmonic" % Harmonic mean: Good for diffusion
        interfaceValues = (leftValues .* rightValues .* (leftSpacing + rightSpacing)) ./ ...
                          (rightValues .* leftSpacing + leftValues .* rightSpacing);

    case "linear" % Linear mean: Smooth blending
        interfaceValues = (rightValues .* leftSpacing + leftValues .* rightSpacing) ./ ...
                          (leftSpacing + rightSpacing);
    case "weighted" % Weighted mean: Based on cell size
        interfaceValues = (rightValues .* rightSpacing + leftValues .* leftSpacing) ./ ...
                          (leftSpacing + rightSpacing);
end

% Compute coefficients for the matrix
coeffs = interfaceValues ./ (leftSpacing + rightSpacing);

% Build tridiagonal matrix
mainDiagonal = [-coeffs(1); -coeffs(2:end)-coeffs(1:end-1); -coeffs(end)];
matrixCore = diag(coeffs, -1) + diag(mainDiagonal) + diag(coeffs, 1);

% Apply scaling
scalingMatrix = diag(2 ./ cellScaling);
interfaceMatrix = scalingMatrix * matrixCore;
end