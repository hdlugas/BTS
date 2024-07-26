# BTS: Basic Testing System

BTS is a web-based R Shiny application for teaching, learning, and performing statistical analyses common in the biological/healthcare sciences. Its functionality includes:

- Non-parametric and parametric tests for one, two, and paired sample comparisons, including t-tests and Wilcoxon rank sum tests. 
- One-Way and Two-Way ANOVA. Tumor growth analysis. 
- Tests for categorical data including the Chi-squared test, Fisher’s exact test, McNemar’s test, and Cochran-Armitage trend test. 
- Correlation calculation between variables including Pearson, Spearman, and Kendall correlations. 
- Regression analysis including linear and logistic regression. 
- Survival analysis, including Cox Proportional Hazards regression. 
- Multiple comparison correction adjustments to p-values, including Holm’s, Bonferroni, Hommel, Hochberg, Benjamini-Hochberg, and Benjamini-Yekutieli adjustments.

# Installation
To install R and RStudio, see [https://cran.r-project.org/](https://cran.r-project.org/) and [https://posit.co/products/open-source/rstudio/](https://posit.co/products/open-source/rstudio/), respectively.To run BTS, open the BTSv1.R file in RStudio and either (1) press Ctrl+Shift+Enter or (2) click "Run App" in the upper-right corner. This will open the Shiny application in a new window. Note that BTS requires the following R packages, which will be installed when running BTSv1.R the first time:

- shiny
- shinythemes
- Hmisc
- ggplot2
- psych
- car
- qqplotr
- ggpubr
- survival
- survminer
- plyr
- Rmisc
- nlme
- DT
- ISwR

By default, the "About BTS" tab will be displayed with a brief introduction upon launching the BTS app (Figure 1) as seen below. The app has thirteen tabs on the top of the right panel, including the About BTS tab, aligned next to each other. Statistical analyses can be conducted by choosing the corresponding tab for each analysis. Each tab primarily consists of left and right panels. The left panel is the input panel that facilitates data upload and various parameter selections for performing each statistical procedure. The corresponding outputs are displayed on the right output panel. 
![about_BTS](https://github.com/user-attachments/assets/1229c190-5c27-47c7-a799-9d5e42566336)


# Usage
The full BTS tutorial can be found in the PDF document BTS_tutorial.pdf

# Bugs/Questions
If you notice any bugs in this software or have any questions, feel free to reach out to fy7392@wayne.edu.

