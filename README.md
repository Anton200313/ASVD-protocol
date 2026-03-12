# ASVD-protocol
An analytical protocol for discovering emergent effects in complex systems. ASVD utilizes Elastic Net regularization to identify high-order interactions and non-linear patterns within large predictor sets. The protocol features automated screening, robust re-estimation via bootstrap, and simple slopes visualization.
## How to Use the ASVD Module

The ASVD module implements a two-step algorithmic protocol for discovering emergent effects and high-order interactions.

### Step 1: Automated Screening
1.  **Variable Selection**: Place your dependent variable in the `Dependent Variable` field and all candidate predictors in the `Independent Variables` field.
2.  **Complexity Setup**: Set the `Max Interaction Order` to define the maximum interaction degree to be tested (e.g., 2 for two-way, 3 for three-way interactions).
3.  **Run Selection**: Choose the regularization parameters (`Alpha` and `Lambda Selection`).
4.  **Get Results**: The output will display a log titled `SELECTED CANDIDATES`. 
    *   **Action**: Copy the list of identified predictors and interactions from this log.

### Step 2: Final Re-estimation and Visualization
1.  **Switch Mode**: Change the `Analysis Step` option from "Elastic Net Screening" to "Bootstrap & Variance decomposition".
2.  **Input Formula**: Paste the copied list of predictors into the **"STEP 2: Paste variables from Step 1 here"** field.
3.  **Analysis**: The module will perform a **Bootstrap** re-estimation, providing stable coefficients, importance metrics, and 95% confidence intervals.
4.  **Interaction Diagnostics**:
    *   Enable **"Show Interaction Plot (Step 2)"** to visualize detected emergent effects.
    *   Adjust the **"Interaction Index"** (1, 2, 3...) to toggle between plots of different identified interactions.
