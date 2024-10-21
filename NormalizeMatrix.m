%* ---------------------------------------------------------------------------------------------------------------------
%*   Function description: 
%    Normalize each non-zero column in the input matrix and retain the zero columns.
%* ---------------------------------------------------------------------------------------------------------------------
%%
function Normalized_U = NormalizeMatrix(U)
[M,~] = size(U);
U_norm = sqrt(sum(abs(U).^2,1));
Normalized_U = U;
Normalized_U(:,U_norm>0) = U(:,U_norm>0)./repmat(U_norm(U_norm>0),M,1);
end
