classdef FVMelectrode
    % FVM class for computing Single Particle Model (SPM) and Multi-Particle Model (MPM) matrices
    % using Finite Volume Method for lithium-ion battery diffusion
    % Reference: http://dx.doi.org/10.1109/CDC.2015.7402829

    methods(Static)

        function [numRadialNodes, radialGridSpacing, particleRadius, diffusion_coefficient_reference] = compute_initialVectors(GeometricParams, TransportParams, electrode)

            % Prefix for electrodes ('neg', 'pos')
            prefix = char(extractBefore(electrode,4));

            % numRadialNodes - Number of radial nodes (nrn or nrp)
            numRadialNodes = GeometricParams.electrode.(electrode).particles.(['number_radial_nodes_' prefix]);
            diffusion_coefficient_reference = TransportParams.electrode.(electrode).(['diffusion_coefficient_reference_' prefix]);

            radialGridSpacing = GeometricParams.electrode.(electrode).particles.(['radialGridSpacing_' prefix]);
            particleRadius = GeometricParams.electrode.(electrode).particles.(['particle_radius_' prefix]);

        end
        %% ------------------------ Base Diffusion Matrix ---------------------- %%
        function Acs_j = Acs_j(numRadialNodes)
            % Computes base diffusion matrix Acs_jn or Acs_jp for a single particle
            % Equation: Tridiagonal matrix for finite volume diffusion (FVM) discretization
            % Inputs: electrode
            % Outputs: Acs_j - Tridiagonal matrix (size nrn * nrn or nrp x* nrp)
            % Symbol: Acs,jn (negative electrode) or Acs,jp (positive electrode) - Base diffusion operator
                        
            subDiagonal = (0:numRadialNodes-2) ./ (1:numRadialNodes-1);
            superDiagonal = [2, (2:numRadialNodes-1) ./ (1:numRadialNodes-2)];

            Acs_j = diag(subDiagonal, -1) - 2 * eye(numRadialNodes) + diag(superDiagonal, 1);
            Acs_j(numRadialNodes, numRadialNodes-1) = 2;
        end

        %% ------------------------ Discrete-Time Diffusion Matrix ---------------------- %%
        function Acs_hat = compute_Acs_hat(numRadialNodes, radialGridSpacing, diffusion_coefficient_reference, timeStep)

            % Computes discrete-time diffusion matrix Acs_hat_n,i or Acs_hat_p,j for a single particle
            % Equation: Acs_hat = Delta t * Acs - I_nrn --> Acs (Scaled Diffusion Matrix)
            % Inputs: Same as Acs
            % params - Structure: numRadialNodes_neg/pos, radialGridSpacing_neg/pos, diffusionCoeff_neg/pos
            % stoichiometry_neg/pos - Scalar (SPM) or vector (MPM): stoichiometry values
            % temperature - Scalar: temperature
            % electrode - String: 'neg' or 'pos'
            % particleIdx - Integer: particle index (1 for SPM)

            % Outputs: Acs_hat - Discrete-time diffusion matrix
            % Symbol: Acs_hat_n,i (negative) or Acs_hat_p,j (positive) - Discrete-time diffusion for a particle         

            Acs_j = FVMelectrode.Acs_j(numRadialNodes);

            particleId = 1;
            numParticles = 1;

            Acs(1+numRadialNodes*(particleId-1):numRadialNodes*particleId,1+numRadialNodes*(particleId-1):numRadialNodes*particleId) = ...
                sparse(kron(eye(1),Acs_j) * ...
                diag((1/radialGridSpacing(1, particleId) ^ 2) * ...
                diffusion_coefficient_reference' * kron(eye(1), ones(1,1))));

            totalNodes = numRadialNodes * numParticles;
           
            Acs_hat = timeStep * Acs - speye(totalNodes);
        end

        %% ----------------- Discrete-Time Diffusion Matrix For Diffusion Function -------------- %%
        function Acs_hat = compute_Acs_hat_diff_func(numRadialNodes, radialGridSpacing, TransportParams, timeStep, electrode, ThermalParams, stoichiometry, temperature)
            % Computes scaled diffusion matrix Acs,n,i or Acs,p,j for a single particle
            % Equation: Acs,n,i = (D_s,n,i / Delta r_n,i^2) * A_cs,jn (base diffusion)
            % Inputs:
            %   params - Structure: numRadialNodes_neg/pos, radialGridSpacing_neg/pos, diffusionCoeff_neg/pos
            %   stoichiometry_neg/pos - Scalar (SPM) or vector (MPM): stoichiometry values
            %   temperature - Scalar: temperature
            %   electrode - String: 'neg' or 'pos'
            %   particleIdx - Integer: particle index (1 for SPM)
            % Outputs: A_cs - Scaled diffusion matrix (size n_r,n x n_r,n or n_r,p x n_r,p)
            % Symbol: A_cs,n,i (neg) or A_cs,p,j (pos) - Continuous diffusion for a particle

            if nargin < 7
                temperature = ThermalParams.cell.('temperature_reference_cell');
            end

            calTransport = TransportParametersCalculator();

            diffusion_coeff = calTransport.compute_electrode_diffusion_coefficient(ThermalParams, TransportParams, electrode, stoichiometry, temperature);

            Acs_j = FVMelectrode.Acs_j(numRadialNodes);

            particleId = 1;
            numParticles = 1;

            Acs(1+numRadialNodes*(particleId-1):numRadialNodes*particleId,1+numRadialNodes*(particleId-1):numRadialNodes*particleId) = sparse(kron(eye(1),Acs_j) * ...
                diag((1/radialGridSpacing(1, particleId) ^ 2) * ...
                diffusion_coeff' * kron(eye(1), ones(1,1))));

            totalNodes = numRadialNodes * numParticles;
           
            Acs_hat = timeStep * Acs - speye(totalNodes);
        end
                
        %% ------------------------ Surface Concentrations matrix ---------------------- %%
        function Acs_bar = Acs_bar(numRadialNodes)     % Acs_bar
            % Computes output matrix Acs_bar for surface concentrations

            numParticles = 1;
            Acs_bar = sparse(kron(eye(numParticles), ...
                                  [zeros(1, numRadialNodes - 1) 1]));
        end
        
        %% ------------------------ Boundary Flux Matrix ---------------------- %%
        function boundaryConcentrationMatrix = compute_boundaryConcentrationMatrix(numRadialNodes, particleRadius, radialGridSpacing)    % Bcs
            % Computes boundary flux matrix for both electrodes (SPM or MPM)
            % It is Bcs
            % Inputs:
            %   params - Structure with fields:
            % Outputs:
            %   boundaryFluxTotal - Block-diagonal matrix of boundary flux vectors

            numParticles = 1;
            particleId = 1;
            Bcs_blocks = cell(1, numParticles);
            Bcs_blocks{particleId} = sparse(-2 * (radialGridSpacing(1, particleId) + ...
                particleRadius(1, particleId)) / (radialGridSpacing(1, particleId) * ...
                particleRadius(1, particleId)) * ...
                [zeros(numRadialNodes-1, 1); 1]);
            boundaryConcentrationMatrix = blkdiag(Bcs_blocks{:});
        end

        %% ------------------------ Discrete-Time Boundary Flux Matrix ---------------------- %%        
        function boundaryConcentrationDiscreteTimeMatrix = compute_boundaryConcentrationDiscreteTimeMatrix(numRadialNodes, particleRadius, radialGridSpacing, timeStep)
            % Computes discrete-time boundary flux matrix boundary Concentration Discrete-Time Matrix
            boundaryConcentrationMatrix = FVMelectrode.compute_boundaryConcentrationMatrix(numRadialNodes, particleRadius, radialGridSpacing);
            boundaryConcentrationDiscreteTimeMatrix = boundaryConcentrationMatrix * timeStep;
        end
        
        %% ------------------------ total discrete-time diffusion matrix ---------------------- %%
        function Acs_hat_total = A_cs_hat_total(GeometricParams, transportParams, stoichiometry_neg, stoichiometry_pos, temperature)
            % Computes total discrete-time diffusion matrix Acs_hat_total for all particles

            Acs_hat_neg = FVMelectrode.Acs_hat(GeometricParams, transportParams, GeometricParams, 'negative', stoichiometry_neg, temperature);
            Acs_hat_pos = FVMelectrode.Acs_hat(GeometricParams, transportParams, GeometricParams, 'positive', stoichiometry_pos, temperature);

            Acs_hat_total = blkdiag(Acs_hat_neg, Acs_hat_pos);
        end

        %% ------------------------ Electrode FVM matirces ---------------------- %%
        function preComputedMatrix = compute_electrodeFVMMatirces(GeometricParams, TransportParams, timeStep, electrode)

            [numRadialNodes, radialGridSpacing, particleRadius, diffusion_coefficient_reference] = ...
                FVMelectrode.compute_initialVectors(GeometricParams, TransportParams, electrode);

            Acs_bar = FVMelectrode.Acs_bar(numRadialNodes);

            boundaryConcentrationDiscreteTimeMatrix = ...
                FVMelectrode.compute_boundaryConcentrationDiscreteTimeMatrix(numRadialNodes, ...
                particleRadius, radialGridSpacing, timeStep);

            Acs_hat = FVMelectrode.compute_Acs_hat(numRadialNodes, ...
            radialGridSpacing, diffusion_coefficient_reference, timeStep);                            

            preComputedMatrix.Acs_bar = Acs_bar;
            preComputedMatrix.boundaryConcentrationDiscreteTimeMatrix = boundaryConcentrationDiscreteTimeMatrix;
            preComputedMatrix.Acs_hat = Acs_hat;

        end
    end
end

