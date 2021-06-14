% function [K,M] = assemble_matrix(p,tri)
% assemble the stiffness matrix K and the mass matrix M
% function [K,M] = assemble_matrix(p,tri,varargin)
% parameter_name = {'diff_coef'};
% default_value = {1};
% diff_coef can be a constant or a row vector at the triangle centers of
% mass, but it cannot be a column vector. 

% Copyright: Shaoying Lu and Yingxiao Wang 2012
function [K,M] = assemble_matrix(p,tri,varargin)
parameter_name = {'diff_coef'};
default_value = {1};
diff_coef = parse_parameter(parameter_name, default_value, varargin);


c = diff_coef; a = 1; f=0;
% c the diffusion coefficient, which can be a row vector of the same
% length as tri
% c can also contain 3 rows representing a 2x2 symmetric matrix
% at the centers of mass for the trianglulation
% c = [   ... c(1,1) ...   ]
%            ... c(1,2) ...   
%            ... c(2,2) ...   
[K,M,~] = assema(p,tri,c,a,f);
%
return;
% mass lumping
% M = diag(sum(M));

