---
title: 'BTS: Basic Testing System - an R/Shiny-Based Web Tool for Statistical Analysis'
tags:
- Shiny
- R
- teaching statistics
- learning statistics
- statistical software for biology
date: "25 July 2024"
affiliations:
- name: Biostatistics and Bioinformatics Core, Karmanos Cancer Institute, Department of Oncology, Wayne
    State University School of Medicine, Detroit, MI, USA
  index: 1
authors:
- name: Hunter Dlugas
  orcid: 0000-0002-6819-0045
  affiliation: 1
- name: Janaka S. S. Liyanage
  orcid: 0000-0001-9882-0362
  affiliation: 1
- name: Seongho Kim
  orcid: null
  affiliation: 1
date: 26 July 2024
bibliography: paper.bib
editor_options: 
  markdown: 
    wrap: 72
---

# Summary

Statistical analysis is crucial for conducting rigorous and reproducible research across various fields, including cancer and biological research. Fundamental data-driven decisions can be effectively made using various statistical methods. For example, categorical variables can be summarized in a contingency table, and associations among them can be quantified with the Chi-Square Test, Fisher’s Exact Test, McNemar’s Test, and the Cochran-Armitage trend test. T-tests compare two groups to assess whether their means differ significantly, while ANOVA extends this analysis to multiple groups to determine significant differences among them. Regression analysis helps understand relationships between independent variables and continuous dependent variables, predicting outcomes based on predictor variables and quantifying predictor variable significance. Survival curves are essential for analyzing time-to-event data, such as patient survival times, offering vital information for medical research and practice. It is increasingly important for researchers and professionals from all backgrounds - including those in the biological sciences - to learn these common statistical concepts and methodologies and to apply them effectively. Despite the development of numerous courses and workshops to teach these critical skills, the availability of user-friendly statistical software poses a barrier to the effective learning and implementation of common statistical methods. Many existing statistical software packages are not accessible to those lacking a proper background in statistics or computer coding. To address this gap, we introduce ‘BTS: Basic Testing System’, an R/Shiny software package designed specifically for essential statistical analysis functions in biological research. BTS aims to provide an intuitive and accessible tool for both instructors teaching statistics and learners acquiring these essential skills. Our software emphasizes ease of use, making it suitable for individuals with limited statistical or programming experience, thereby supporting the broader adoption of sound statistical practices in scientific research. 

# Statement of need

Comprehensive education in statistical methods is crucial for trainees, including researchers in the biological sciences [@Weissgerber:2016, Colon:2011]. Mastery of methods such as contingency table analysis, T-tests, ANOVA, regression analysis, and time-to-event analysis is essential for making data-driven decisions and inferences in both laboratory and clinical settings. Researchers with this statistical expertise enhance data analysis quality, critical thinking, evidence-based decision-making, and the ability to conduct rigorous, reproducible research [@Kohrs:2023]. Statistical software is vital in this educational process, as it is the tool used to implement data analyses. Such software simplifies procedures, automates calculations, and generates visualizations which reduces the learning curve and promotes self-learning. Despite these benefits of high-quality statistics education and researcher fluency in statistical software, user-friendly, free, and efficient statistical software remains limited [@JSSv010i05]. Indeed, our experience teaching a statistics course for a cancer biology program revealed that many standard statistical software packages such as SAS (SAS Institute Inc., Cary, NC, USA), SPSS (IBM Corp, Armonk, NY, USA), STATA (StatCorp LLC, College Station, TX, USA), and R (R Core Team) require extensive training, posing challenges for beginners. Several R/Shiny packages for statistical data analysis have been introduced, each with unique features to enhance user experience and accessibility. Radiant [@radiant:2023] provides a comprehensive interface for business analytics, while vvdoctor [@vvdoctor:2024] offers a robust tool for verifying datasets. ShinyItemAnalysis [@Martinkova:2023] is tailored for psychometric analysis and test scoring, easystats [@Daniel:2022] simplifies statistical analysis with a user-friendly interface, and EZR [@EZR:2013] focuses on medical statistics, making complex analyses more approachable for health researchers. These packages collectively contribute to the expanding ecosystem of accessible statistical tools, aiding both novice and experienced users in performing advanced data analyses.

However, in our experience, all these existing packages require a steep learning curve and/or focus on more advanced data analysis methods which are beyond the needs for many researchers in the biological sciences. Thus, there is a pressing need for statistical software that combines robust analytical capabilities with user-friendly interfaces to enhance understanding and enthusiasm for statistics. Such advancements can empower researchers-in-training, foster innovation, and accelerate discovery in the biological sciences.

To fill this gap, we developed 'BTS: Basic Testing System,' an R/Shiny-based statistical software app for teaching, learning, and performing basic statistical analyses in the biological sciences. This Shiny application was built using R [@R:2020] and RStudio [@RStudio:2019] statistical software, along with Shiny [@Chang:2024]. The BTS was thoroughly tested as the supporting statistical software for teaching the aforementioned statistical course in biological studies. The BTS can be used by researchers with non-statistical backgrounds to execute essential statistical methodologies efficiently for conducting their research. The BTS application and a tutorial illustrating its utilization are available on GitHub ([https://github.com/hdlugas/BTS](https://github.com/hdlugas/BTS)).

# Functionality

The following statistical analyses can be conducted using BTS:
- Non-parametric and parametric tests for one, two, and paired sample comparisons, including t-tests and Wilcoxon rank sum tests.
- One-Way and Two-Way ANOVA. Tumor growth analysis.
- Tests for categorical data, including Chi-square test, Fisher’s exact test, McNemar’s test, and Cochran-Armitage trend test.
- Correlation calculation between variables, including Pearson, Spearman, and Kendall correlations.
- Regression analysis, including linear and logistic regression.
- Survival analysis, including Cox Proportional Hazards regression.
- Multiple comparison correction adjustments to p-values, including Holm’s, Bonferroni, Hommel, Hochberg, Benjamini-Hochberg, and Benjamini-Yekutieli adjustments.

# Acknowledgements

This research was partially funded by the National Institutes of Health (R21GM140352) and the National Cancer Institute (P30CA022453).

# References
