

####@# 5.1

+ research paper
    + idea
        + probability of losses given a credit portfolio L
        + L: set of bonds, with interest
        + L is probabilistic model 
    + L(Z, \Sigma): 2 normal(0, 1) distribution 
        + first 2 chapters
        + biggest loss: everything defaults
        + also loss: reduce credit rating (AA -> A)
        + gain: increase credit rating
    + try to estimate with monte carlo 
    ```
    P(L(Z, \Sigma) > e) = E(1_{L(Z, \Sigma) > \epsilon})
    \approx \frac{1}{MN} \sum_{n=1}^N \sum_{n=1}^M 1_{L(Z, \Sigma) > \epsilon}
    ```
    + problem, MC runs too slowly
+ importance sampling 
    + how fast MC runs depends on variance
        + i..e Var_f(g(x))
    ```
    E(f(x)) = \int g(x) f(x) dx
    = \int \frac{g(x) f(x)}{\hat{f}(x)} \hat{f}(x) dx
    = E_{\hat{f}}{ \frac{g(x) f(x)}{\hat{f}(x)} }
    ```
    + potential source of error 
        + 2 MC on lhs and rhs, they should be approximately equal
        + but there seems to be some consistent bias, not as close as wed expect
    + importance sampling reduce variance 
        + var_{\hat{f}}{\frac{g(x) f(x)}{\hat{f}(x)}}
+ variance reduction method
    + why variance reduction for glasserman li doesnt work
+ 2-state

##### 5.21

+ stock return
    + log normal model
    + not that accurate, underestimate large losses
    + want thicker tails
+ motivation for H
    + P(y \leq H_{c(n)}^n)
+ why not model the covariance matrix and in this way sampling uses only 1 layer
    + 2 level useful
+ numerical analysis computing resources 
    + loginto appS0, appS1, ... then login to
    + compS0, compS1, ...
        + `ssh kjr@compS0.cs.toronto.edu`
+ TODOS:
    + check 2 simple MC algorithm 2.1 2.2 yield same result 
        + the result are consistently different, but try 
            + run on larger sample sizes 
            + try run 2 algorithm on the same saved parameters, maybe they are different and somehow, however unlikely, effect the result
        + also some MC estimate in 1 run yield 0 ??!!

##### 5.28

+ TODO:
    + verify 2.1 2.2 again with CI, on large sample sizes
    + re-read glasserman lee
+ since model is simulated, can compute the real statistics not relying on sampling
+ try different sample size -> observe varied variance ...
+ instead storing mean(...)


``` 
% naive bernoulli 2 lvl mc as estimator for true value ...

% 500 sample z 10000 sample e
m =
0.0108498  % log(m) =  -4.5236 

s =
0.08176819905509370367000561334108
mci =
0.0108498 - 0.0032130476786333435318806632420398*5^(1/2)    % 0.003665206976
0.0032130476786333435318806632420398*5^(1/2) + 0.0108498    % 0.01803439302
sci =
0.076994856480868758566102903459682
0.087177301847900966319842520379664

% 740 sample z 10000 sample e
m =
0.0096022972972972972972972972972973    % log(m) = -4.6458
s =
0.074992749366207626482320314398412

mci =
 0.0096022972972972972972972972972973 - 0.00039790327266436918932749219745937*185^(1/2)    % 0.004190227669
 0.00039790327266436918932749219745937*185^(1/2) + 0.0096022972972972972972972972972973    % 0.01501436693
sci =
 0.071356932703114390005142823954217
 0.079021889730105408411683231859265

 % 1000 sample z 10000 sample e

m =
0.0086053

s =
0.07264915664941935079140065088676

mci =
0.0086053 - 0.0014256245220952668699852953629357*10^(1/2)   % 0.004097079422
0.0014256245220952668699852953629357*10^(1/2) + 0.0086053   % 0.01311352058

sci =
0.069598789223672154152100560846431
0.075981236183767535161490983828342
```
