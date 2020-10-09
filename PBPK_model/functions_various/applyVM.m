function [ Mout ] = applyVM( M, V, in )

% Nicola Melillo 20/01/2017

% apply vector to matrix
% multiply each row (in=1) or column (in=2) of the matrix M for the vector
% V.
% * M: matrix to be multiplied.
% * V: vector, must be of the length of M's columns (in=2) or rows (in=1). V
% must be column.
% * in: if in=1 I multiply the vector for the rows and if in=2 for the
% columns.

%%

nrows = size(M,1);
ncol = size(M,2);
Mout = zeros(nrows,ncol);

if in==1 % multiply for the rows
    
    for i=1:nrows
        Mout(i,:) = M(i,:).*V'; 
    end
    
elseif in==2 % multiplu for the columns
    
    for i=1:ncol
        Mout(:,i) = M(:,i).*V; 
    end
    
end


end

