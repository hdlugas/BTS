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
To download BTS, one can click the green "Code" box at the top-right of this repository, and then click "Download ZIP". Once this ZIP file is downloaded on your local computer, extract the contents.

To install R and RStudio, see [https://cran.r-project.org/](https://cran.r-project.org/) and [https://posit.co/products/open-source/rstudio/](https://posit.co/products/open-source/rstudio/), respectively. To run BTS, open the BTSv1.R file in RStudio and either (1) press Ctrl+Shift+Enter or (2) click "Run App" in the upper-right corner. This will open the Shiny application in a new window. Note that BTS requires the following R packages, which will be installed (if not already installed) when running BTSv1.R the first time:

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

By default, the "About BTS" tab will be displayed with a brief introduction upon launching the BTS app (Figure 1) as seen below. The app has thirteen tabs on the top of the right panel, including the "About BTS" tab. Statistical analyses can be conducted by choosing the corresponding tab for each analysis. Each tab primarily consists of left and right panels. The left panel is the input panel that facilitates data upload and various parameter selections for performing each statistical procedure. The corresponding outputs are displayed on the right output panel.  
![intro_window](https://github.com/user-attachments/assets/c4f6ebed-591a-463b-beb1-1852fc0f11fd)


# Usage
As an example of how to use BTS in the context of survival analysis for example, we can first upload our data by opening BTS, clicking on "Upload Data", clicking on Browse in the top of the left panel, and selecting the toydata_survivalanalysis_veteran.csv file in the example_data directory. We then click on the "Survival Curve" tab. Now if we want to construct a univariable Cox Proportional Hazards Regression model for example with the celltype group (i.e. squamous, smallcell, adeno, and large) as the predictor variable, we can choose diagtime as the survival time varaible (which indicates the amount of time to event), trt as the survival status variable (1 indicates no event, 2 indicates event), Input Unit as months, Output Unit as months, X-axis label as "Months After Diagnosis", Y-axis label as "Overall Survival", Time point for a survival rate as 1, Group comparison as yes, and group variable as celltype. The left panel of BTS will then appear as it does in the screenshot below. The distribution of the celltype groups, the Kaplan-Meier curve, and the hazard ratios and p-values of each indicator variable for the non-reference celltype groups are found in the right panel.

![side_panel_survival_analysis](https://github.com/user-attachments/assets/4b9538dd-2cd9-4f8e-ad2b-19ce304fe283)


The full BTS tutorial can be found in the PDF document BTSv1_user_tutorial.pdf.


