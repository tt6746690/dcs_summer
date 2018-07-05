
#### TODO
 

+ use beamer and sharelatex 
+ 10 slides     
    + 1 intro 
    + 1 ending 
    + 1 overview 
    + 1 finance terminology + definition stuff
    + 1 MC
    + 1 IS
    + ~2 problem setup 
        + copula model
        + what is the problem
            + how to simulate rare event better 
        + motivation for MC IS
            + rare event simulation 
    + ~2 steps and solutions 
        + introduce copula factor model
        + how glasserman and li works
    + experiments/tests/future direction



#### abstract 

In finance, credit risk is the risk of default on a debt that may arise when a borrower fails to make required payments. Monte Carlo simulation is often used to measure such risk, in other words the probability of losses from a large number of defaults. We will briefly introduce Monte Carlo simulation as a tool for such computation. In reality, default rarely occur but when they do, often incur significant losses. We will explain how importance sampling makes Monte Carlo simulation more effective at computing probabilities of rare events. We will then describe how Monte Carlo simulation and importance sampling is applied to the Gaussian copula factor model, a predominant model for evaluating credit risk. Finally, we will talk about current directions


Credit risk is one of the fundamental risks for financial entities to manage. It refers to the possible losses due to defaults or downgrades of ratings for obligors in a given portfolio. One outstanding question is to measure the probability that this possible loss will exceed a given amount for a specified portfolio. In real life, the measurement of credit risk is often a rare-event simulation problem, because risk management tends to focus on rare but significant losses which result from a large number of defaults. This nature motivates the application of importance sampling(IS) techniques on simulations. In this presentation, we will generally walk through the application of IS procedures under Gaussian Copula Factor Model, the most widely used model on credit risk estimations nowadays.
