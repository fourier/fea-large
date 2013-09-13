function T = residual(el0,el)
% function T = construct_t(el0,el)
% Function for constuction of internal force vector T
% el - element in current configuration
% q - displacements from initial configuration(for deformation gradient)
% T = int(B^T * Sigma)
    
    T = integrate5nodes(el,@tvector,el0);
    
