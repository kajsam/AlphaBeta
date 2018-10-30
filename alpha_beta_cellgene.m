function [pi] = alpha_beta_cellgene(X, W, H, notW, notH, fig_nr)

% Kajsa Mollersen (kajsa.mollersen@uit.no) 29th October 2018

% Block by block

% Maximum likelihood estimation of alpha and beta in the regression model
% with Bernoulli distribution and logit link function. For a thorough
% presentation, see alphabeta.tex.

% The Bernoulli pi consists of a probability of success, called 'effect',
% which is p for all entries where A_ij = 1, and (1 - p) for all entries
% where A_ij = 0, and in addition a cell effect, called 'cell_effect',
% which is added to the 'effect'. And equivalently a gene effect. 

% Input:        X (n x d)   - observed matrix
%               A (n x d)   - class membership (0 or 1)

% The alpha's and beta's are estimated as sums, and the four possible terms
% of the sums are the combinations of 0 and 1 for X_ij and A_ij

if ~(islogical(X))
  disp('Logical, please')
  return
end
     
[n,d] = size(X);

% The overall effect, common for each entry in X
p = 1 - sum(sum(xor(X,W*H)))/(n*d)

K = size(W,2)
Pi = zeros(n,d);
for k = 1: K
  A = logical(W(:,k)*H(k,:));
  n11 = sum(X(A));
  denom = sum(A(:));
  pi1 = n11/denom;
  Pi(A) = pi1;
end

for k = 1: K+1
  notA = logical(notW(:,k)*notH(k,:));
  n10 = sum(X(notA));
  denom = sum(notA(:));
  pi0 = n10/denom;
  Pi(notA) = pi0;
end

if fig_nr
  figure(fig_nr), imagesc(Pi), colormap(gray), title('\Pi')
end
if any(find(Pi == 0)) || any(find(Pi == 1))
  disp('/pi = 0 or 1')
  % return
end

pi0_g = zeros(1,d);
pi1_g = zeros(1,d);

pj = zeros(1,d);
for j = 1: d
  
  pi0_g(j) = sum(Pi(~A(:,j),j))/sum(~A(:,j));
  pi1_g(j) = sum(Pi(A(:,j),j))/sum(A(:,j));
  if sum(~A(:,j)) == 0
    pi0_g(j) = max(pi1_g(j) - 2*p + 1,1/n);
  elseif sum(A(:,j)) == 0
    pi1_g(j) = min(2*p - 1 + pi0_g(j),(n-1)/n);
  end
  pj(j) = (pi1_g(j) - pi0_g(j)+1)/2;
end
p = mean(pj(pj>0));
for j = 1: d
  if sum(~A(:,j)) == 0
    pi0_g(j) = max(pi1_g(j) - 2*p + 1,1/n);
  elseif sum(A(:,j)) == 0
    pi1_g(j) = min(2*p - 1 + pi0_g(j),(n-1)/n);
  end
end



pi = [pi0_g; pi1_g];


% The d gene effects
p_gene = zeros(1,d);
for j = 1: d
  if sum(A(:,j)) > sum(~A(:,j))    
    p_gene(j) = pi1_g(j) - p;
  else
    p_gene(j) = pi0_g(j) - (1 - p);
  end
end

if fig_nr
  figure, imagesc(repmat(p_gene,10,1)), colormap(gray), 
  title('Estimated gene effect'), drawnow
end
  