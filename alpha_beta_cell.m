function [alphabeta_c, cell_effect, Pi] = alpha_beta_cell(X, A, fig_nr)

% Kajsa Mollersen (kajsa.mollersen@uit.no) 15th October 2018

% Similar to the alpha_beta function, but here ignoring the gene effect

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

if ~(islogical(X) && islogical(A))
  disp('Logical, please')
  return
end
     
[n,d] = size(X);


% Cell effect: There are n cells, each with different alpha and beta. Once
% alpha and beta are known, the Bernoulli pi can be extracted.

% We assume that 0 < pi < 1, and therefore add the machine epsilon.
% That did not have a good effect. Let's set the minimum count to 1. 

alpha_c = zeros(1,n);
beta_c = zeros(1,n);
pi0_c = zeros(1,n);
pi1_c = zeros(1,n);

t = zeros(2,n);
for i = 1 : n
  n00 = max(sum(~X(i,:) & ~ A(i,:)),1); % Let it be minimum 1 to avoid 1/0
  n01 = max(sum(~X(i,:) & A(i,:)),1); % Let it be minimum 1 to avoid 1/0
  n10 = max(sum(X(i,:) & ~A(i,:)),1); % Let it be minimum 1 to avoid log 0
  n11 = max(sum(X(i,:) & A(i,:)),1);  % Let it be minimum 1 to avoid log 0
  alpha_c(i) = log(n10/(n00+eps) + eps);
  beta_c(i) = log(n11/(n01+eps) + eps) - log(n10/(n00+eps) + eps);
  pi0_c(i) = exp(alpha_c(i))/(1 + exp(alpha_c(i))); 
  pi1_c(i) = exp(alpha_c(i) + beta_c(i))/(1 + exp(alpha_c(i)+ beta_c(i))); 
  t(1,i) = log((n00+n10)/(n00+eps) + eps);
  t(2,i) = log((n01+n11)/(n01+eps) + eps);
end

alphabeta_c = [alpha_c; beta_c; t];

% The (n x d) probability matrix can be calculated
% pi0 if A_ij = 0, pi1 if A_ij = 1
Pi = zeros(n,d);
for i = 1 : n
  for j = 1 : d
    Pi(i,j) = pi0_c(i) + (pi1_c(i) - pi0_c(i))*A(i,j);
  end
end

% Maybe I should check for 0's and 1's. Don't know what to do about them,
% though.
if any(find(Pi == 0)) || any(find(Pi == 1))
  disp('/pi = 0 or 1')
  return
end

% The overall effect, common for each entry in X
effect = 1 - sum(sum(xor(X,A)))/(n*d); 

% The n cell effects
cell_effect = zeros(n,1);
for i = 1: n
  pi0 = 0;
  pi1 = 0;
  if sum(A(i,:))    
    pi1 = sum(Pi(i,A(i,:)))/sum(A(i,:));
  end
  if sum(~A(i,:))
    pi0 = sum(Pi(i,~A(i,:)))/sum(~A(i,:));
  end
  cell_effect(i) = (pi0 + pi1 - 1)/2;
end

if fig_nr
  figure(fig_nr), imagesc(repmat(cell_effect,1,10)), colormap(gray), title('Cell effect')
end








    
      

  