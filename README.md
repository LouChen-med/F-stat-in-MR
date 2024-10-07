
# F-stat-in-MR

**Arguing That the F-Value in MR Studies Is Meaningless**

![image](https://github.com/user-attachments/assets/ab395487-e3b8-4819-8684-741a8fcc5c77)

This is a paragraph from [Nature Human Behavior](https://doi.org/10.1038/s41562-024-01879-8), which discusses the formula for calculating the F-value of an SNP and also provides appropriate references. It is widely known in related research that $F > 10$ implies that the instrumental variable has sufficient strength. However, these widely applied formulas inevitably yield $F > 10$.

## Converting the F-calculation from the Nature Human Behavior paper into R code

```r

df <- read_excel(file_path, sheet = "Harmon_ex")


mr_F <- plyr::ddply(df, c("id.exposure", "id.outcome"), 
                      function(x1) {
                        x <- subset(x1, mr_keep)
                        if (nrow(x) == 0) {
                          message("No SNPs available for MR analysis of '", 
                                  x1$id.exposure[1], "' on '", x1$id.outcome[1], 
                                  "'")
                          return(NULL)
                        }else {
                          message("Analysing '", x1$id.exposure[1], "' on '", 
                                  x1$id.outcome[1], "'")
                        }
                        x<-dplyr::select(x,c(beta.exposure,se.exposure,samplesize.exposure))
                        k<-nrow(x)
                        #print(k)
                        # print(x)
                        # head(x)
                        x$R2<-((x$beta.exposure)^2)/(((x$beta.exposure)^2)+((x$se.exposure)^2)*x$samplesize.exposure)
                        x$Fstat<- ( x$R2/(1- x$R2))*((x$samplesize.exposure-k-1)/k)
                        x$`F-statistic`<-sum(x$Fstat)
                        x <- head(x, 1) 
                        x <- dplyr::select(x, `F-statistic`) 

                        return(x)
                      })

fwrite(mr_F,save_path,sep="\t")




```

## Single SNP R-squared Formula

The formula for calculating the proportion of variance explained by an SNP is:

$R^2 = \frac{\beta^2}{\beta^2 + SE^2 \times N}$

where $\beta$ is the effect size, $SE$ is the standard error, and $N$ is the sample size.

## F-statistic for a Single SNP

The formula for the F-statistic of a single SNP is:

$F = \frac{N - K - 1}{K} \cdot \frac{R^2}{1 - R^2} = \frac{N - K - 1}{K} \cdot \frac{\beta^2}{SE^2 \times N}$

where $K$ is the number of SNPs, and $N$ is the sample size.

## Conditions
In general, K and N have the following size ranges:

$K < 200$

$N > 10,000$

thus

$N - K - 1 \approx N$

## Cumulative F-statistic for All SNPs in a Cluster

The cumulative F-statistic for all SNPs in a cluster is:

$F_{sum} \approx \frac{1}{K} \left( \frac{\beta_1^2}{SE_1^2} + \frac{\beta_2^2}{SE_2^2} + \cdots + \frac{\beta_K^2}{SE_K^2} \right)$

This measures the cumulative strength of association for $K$ SNPs in one exposure.

## Significance Threshold

The significance threshold is calculated as:

$\frac{\beta}{SE} = Z$

A suggested threshold is $Z \geq 3.16$, corresponding to a significance level of $p = 0.001565$. For a typical SNP used as instrumental variable, $p$ should be less than $1 \times 10^{-5}$, which implies that $Z^2$ will be greater than 10. Thus, for any given SNP, the formula can be simplified to:

$$
\frac{\beta^2}{SE^2} > 10 \implies F_{sum} \geq 10
$$

This implies that for a given exposure, regardless of the number of SNPs, if the $p$-value for an SNP is less than $0.001565$, then $F_{sum}$ will be greater than 10.

In summary, the calculation of $F$ is largely a trick in many paper.

## Plot of P-Z relationship

![image](https://github.com/user-attachments/assets/1e5c3505-c5da-4245-a408-7508892a978b)

code for the plot
```r
# Load necessary packages
library(ggplot2)

# Generate a range of Z values, assuming Z values range from 0 to 10
z_vals <- seq(0, 10, by = 0.01)

# Calculate P-values using the cumulative distribution function of the normal distribution
p_vals <-stats::pnorm(z_vals, lower.tail = FALSE) * 2

# Put Z values and P values into a data frame
df <- data.frame(Z = z_vals, P = p_vals)

# Calculate the P-value corresponding to Z = 3.162278
z_specific <- 3.162278
p_specific <- stats::pnorm(z_specific, lower.tail = FALSE) * 2

# Use ggplot2 to plot the curve of P-values against Z values
ggplot(df, aes(x = Z, y = P)) +
  geom_line(color = "blue") +
  scale_y_log10() +  # P-values are typically plotted on a logarithmic scale
  labs(title = "P-value vs Z-score", 
       x = "Z-score", 
       y = "P-value") +
  theme_minimal() +
  
  # Add horizontal lines to mark the positions of P = 0.05 and P = 1e-3
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1e-3, linetype = "dashed", color = "green") +
  
  # Add a vertical line to mark the position of Z = 3.162278
  geom_vline(xintercept = z_specific, linetype = "dashed", color = "purple") +
  
  # Add text annotations to label the P-value on the vertical line
  annotate("text", x = z_specific, y = p_specific, 
           label = paste0("Z = 3.16\nP = ", round(p_specific, 6)), 
           vjust = -0.5, color = "purple") +
  
  # Add additional text annotations
  annotate("text", x = 9, y = 0.05, label = "P = 0.05", vjust = -1, color = "red") +
  annotate("text", x = 9, y = 1e-3, label = "P = 1e-3", vjust = -1, color = "green")

```






