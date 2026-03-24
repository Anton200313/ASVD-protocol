# ASVD-Protocol (Automatic Selection and Variance Decomposition)
ASVD is an advanced analytical framework designed for discovering emergent effects and high-order interactions within complex systems. The protocol integrates the feature-selection power of Elastic Net regularization with the statistical robustness of Non-parametric Bootstrap estimation and variance decomposition for final model interpretation.
### Core Features
1. **Automatic Selection**. Effortlessly extract non-linear patterns and interactions (2nd, 3rd order, and beyond) from large candidate predictor sets using penalized regression.
2. **Elastic Net Regularization**. Efficiently handles multicollinearity and overfitting by balancing Ridge and Lasso penalties (adjustable alpha and lambda).
3. **Robust Re-estimation**. Validates coefficient stability and 95% Confidence Intervals (CI) using a bootstrap resampling approach.
4. **Variance Decomposition**. Quantifies the relative importance of each predictor and interaction by calculating their unique contribution to the total R-squared.  
5. **Simple Slopes Diagnostics**. Automatic visualization and statistical testing of slopes to analyze moderation effects across different levels (Low, Mean, High).
## User Manual: the two-step workflow
The ASVD module implements a rigorous two-step algorithmic protocol to ensure only the most reliable effects are included in the final model.
### Step 1. Automatic Selection (Feature Screening)
In this phase, the algorithm "sifts" through all possible variable combinations, retaining only those with genuine predictive power.
1. **Variable Assignment**. Place your target variable in the Dependent Variable field and all potential predictors in the Independent Variables field.
2. **Define Complexity**. Set the Max Interaction Order.
*Example*. Setting this to 3 will test all main effects, two-way interactions and three-way interactions.
3. **Regularization Parameters**:
*Alpha*. Set to 1.0 for Lasso (maximum sparsity), 0.0 for Ridge, or 0.5 for a balanced Elastic Net search.
*Lambda Selection*. Choose min for maximum accuracy or 1se for the most parsimonious (simplest) model.
4. **Extract Formula.** The module generates a log titled "SELECTED CANDIDATES".
*Action*. Click the generated formula string in the text field to select it and Copy it.

### Step 2. Final Re-estimation & Variance Decomposition
Transition from automated discovery to classical statistical validation and effect decomposition.
1. **Switch Mode.** Change the Analysis Step option to "Step 2: Bootstrap & Variance decomposition".
2. **Input Formula.** Paste the copied string into the "STEP 2: Model Formula" field.
3. **Interpret Metrics:**
*Model Performance:* review R-squared, adjusted R-squared and cross-validated R-squared for overfit checking.
*Importance Table:* analyze the contribution of each term and standardized Beta coefficients with 95% Bootstrap CIs.
4. **Emergent Effect Diagnostics.**
*Enable "Show Interaction Plot"* to visualize moderation.
*Use the "Interaction Index" (1, 2, 3...)* to toggle between different detected interactions.
*Review the Simple Slopes Analysis table.* It tests the significance of X’s effect on Y at specific moderator levels (Low -1SD, Mean, High +1SD).

### Data Requirements & Best Practices
1. **Variable Types.** Predictors should be continuous or dummy-coded numeric variables.
2. **Sample Size.** A sample size of N > 100 is recommended for stable bootstrap intervals, especially when testing 3rd-order interactions.
3. **Interpreting Slopes.** If the 95% CI [Lower, Upper] for a slope does not include zero, the effect is statistically significant at that specific moderator level.
