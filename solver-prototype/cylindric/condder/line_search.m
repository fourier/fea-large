function [rtu,eta,eta0] = line_search(R,X,eta0,eta,rtu0,rtu)
% function [rtu,eta,eta0] = line_search(R,X,eta0,eta,rtu0,rtu)
	rtu = 0;
	rtu = X'*R;
	if rtu == 0.0
		return
	end	
	rtu1 = (rtu-rtu0*(1-eta))/(eta*eta);
	Alpha = rtu0/rtu1;
	eta0 = eta;
	
	if Alpha < 0 
		q = (Alpha - sqrt ( Alpha * (Alpha - 4) ) )/ 2.0;
		eta = Alpha/q;
	elseif Alpha < 2
		eta = Alpha/2.0;
	else	
		eta = 1.0;
		rtu = 0;
	end

