

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
    + check 2.1 2.2 the same algorithm