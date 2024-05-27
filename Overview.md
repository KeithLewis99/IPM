---
title: "Overview of model updates from April 2024 onward"
output: 
   # html_document:
   bookdown::html_document2:
      toc: true
      number_sections: true
# output:
#   pdf_document:
#     citation_package: natbib
#   bookdown::pdf_book:
#     citation_package: biblatex
vertical_layout: scroll
bibliography: references.bib
author: Nicholas Gullage

---


  
# Version 1 {-}
   
## Version 1.0 (2024-04-25) {-}

I started the process model from scratch, but took much of the parameter priors and diagnostics from model "cap.v37" from Keith. Below are the Priors, Cohort model, and diagnostic equations used for the first version of the Capelin SS JAGS model.

This version details the "base" model to which modifications are be added. Unless otherwise specified, future versions will maintain every component of this version version to itself.
 
### Priors {-}

The priors for model variance defined as follows:

\begin{equation}
   \sigma_{proc} \sim dnorm(0.01, 20)
\end{equation}
\begin{equation}
   \sigma_{obs} \sim dnorm(0.01, 20)
\end{equation}
\begin{equation}
   \sigma_{mat} \sim dnorm(0.01, 1)
\end{equation}
\begin{equation}
   \sigma_{rec} \sim dnorm(0.01, 20)
\end{equation}

where $\sigma_{proc}$ is process error variance, $\sigma_{obs}$ is observation error variance, $\sigma_{mat}$ is variance in maturity values, and $\sigma_{rec}$ is variance in recruitment before larval density information is available.

Priors for abundances-at-age for the first year are defined as follows:

\begin{equation}
   rec_{y=1} \sim dnorm(15, 1/9)
\end{equation}
\begin{equation}
   N^{a=2}_{y=1} \sim dnorm(12.8, 1/9)
\end{equation}
\begin{equation}
   N^{a=3}_{y=1} \sim dnorm(11.3, 1/9)
\end{equation}
\begin{equation}
   N^{a=4}_{y=1} \sim dnorm(8.2, 1/9)
\end{equation}

Priors for age-dependent processes for survival and recruitment are defined as follows:

\begin{equation}
   \alpha \sim dnorm(0, 100^-2)
\end{equation}
\begin{equation}
   \beta \sim dnorm(0, 100^-2)
\end{equation}
\begin{equation}
   \delta^a \sim dnorm(0.01, 1)
\end{equation}
\begin{equation}
   \gamma^a \sim dnorm(1.29, 2.56)
\end{equation}
\begin{equation}
   \varepsilon^a \sim dnorm(0, 100^-2).
\end{equation}

Here, $\gamma^a$ represents the effect of ice width and $\delta^a$ represents the effect of timing/rate of ice coverage (see [Buren et al., 2014](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0087589)) on total mortality for each age, $\varepsilon^a$ represents the effect of fish condition on mortality for each age, and $\alpha$ and $\beta$ are parameters modeling a linear relationship between recruitment (abundance at age 1) once larval density information is available,

Recruitment for years prior to larval density estimates are defined as a random walk on previous years' recruitment estimates, where priors defined as follows:

\begin{equation}
   rec_y \sim dlnorm(ln(rec_{y-1}), \sigma_{rec}).
\end{equation}

Here, I use a lognormal relationship to bound the value of recruitment to $rec_y \in (0, \infty)$.
Priors for maturities-at-age are defined as follows:

\begin{equation}
   M^a_{est} \sim dbeta(100 M^a_{obs}, 100 (1-M^a_{obs})),
\end{equation}

where $M_{est}$ is the estimated maturity-at-age and $M_{obs}$ is the observed maturity-at-age.

A total mortality component was added to replace survival from previous models, and priors are defined as follows:

\begin{equation}
   \sigma_Z \sim dunif(0.01, 1),
\end{equation}
\begin{equation}
   Z^a_0 \sim dlnorm(0, \tau_Z),
\end{equation}

where $Z^a_0$ is used as an offset (intercept) for total mortality-at-age. I use a lognormal relatioship foo $Z^a_0$ for similar reasons as recruitment, i.e. we expect the *unknown* mortality component to be non-negative, since a negative value for mortality would increase stock size. We define values for natural mortality for ages 1 (i.e. recruitment), 2, and 3. Value for $Z_0$ for age 4 fish are irrelevant because we assume a total mortality at this age.


### Process Equations {-}

Total mortality is used in place of survival from previous model versions, where

\begin{equation}
   S^a_y = exp(-Z^a_y),\hspace{.5in} Z^a_y = -ln(S^a_y).
\end{equation}

We define total mortality as a function of ice coverage, $TI_y$, and Fall capelin condition, $CO_y$, when condition information becomes available (i.e, 1995),

$$
   Z^a_y = 
   \begin{cases}
   4 \gamma^a \frac{TI_y}{\delta^a} (1 - \frac{TI_y}{\delta^a}) + Z^a_0,\hspace{.5cm} y < 1995 \\
   4 \gamma^a \frac{TI_y}{\delta^a} (1 - \frac{TI_y}{\delta^a}) + \varepsilon^a CO_y + Z^a_0,\hspace{.5cm} y \geq 1995. \\
   \end{cases}
$$

Abundance-at-age is defined as in the cohort equation, such that

\begin{equation}
   N^a_y = \mu^{a-1}_{y-1} - Z^{a-1}_{y-1},
   (\#eq:process)
\end{equation}

where Abundances for proceeding ages depend on the *processed* abundances of previous ages. Note, we define *N* as the natural log of abundance, $N^a_y := ln(Abundance)$. Abundance across all ages for the first year are estimated parameters defined by priors (see above). 

Abundances at age 2 (i.e. $N^{a=2}_y$) are calculated similar to Eq. \@ref(eq:process), as

\begin{equation}
   N^{a=2}_y = \rho_{y-1} - Z^{rec}_{y-1}.
\end{equation}

Recruitment is defined by priors (see above) for years where larval density information is unavailable for their respective cohort. When larval density information is available (years 2004-present), values for recruitment, *rec*, are estimated as

\begin{equation}
   rec_y = \alpha + \beta LD_{y-1}.
\end{equation}

where *LD* is larval density.

Observed indices are fitted to process-corrected abundance estimates as

\begin{equation}
   I^a_y \sim dnorm(\mu^a_y, \sigma_{obs}).
\end{equation}

Additionally, priors are fitted against the projected cohort abundance (see below) to estimate process error in abundances-at-age,

\begin{equation}
   \mu^a_y \sim dnorm(N^a_y, \sigma_{proc}),
\end{equation}

where $\mu$ is the estimate of abundance accounting for process error. Priors are also fitted against recruit abundance,

\begin{equation}
   \rho_y \sim dnorm(rec_y, \sigma_{proc}).
\end{equation}



### Diagnostics {.unlisted .unnumbered .tabset}



#### Stock Trends {- .tabset}

##### Index Estimates {-}

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

##### Maturity Estimates {-}

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

##### Recruitment Estimates {-}

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

##### SPAY {-}

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

#### Summaries {- .tabset}

##### Parameters {-}

<div class="tabwid"><style>.cl-4e1eeec2{}.cl-4e02b554{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-4e02b55e{font-family:'Arial';font-size:6.6pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;position: relative;bottom:3.3pt;}.cl-4e02b55f{font-family:'Arial';font-size:6.6pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;position: relative;top:3.3pt;}.cl-4e1a1ee2{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4e1a1eec{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4e1a2f86{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4e1a2f90{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4e1a2f91{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4e1a2f9a{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4e1a2f9b{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4e1a2fa4{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.tabwid {
  font-size: initial;
  padding-bottom: 1em;
}

.tabwid table{
  border-spacing:0px !important;
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  border-color: transparent;
  caption-side: top;
}
.tabwid-caption-bottom table{
  caption-side: bottom;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td, .tabwid th {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
.katex-display {
    margin: 0 0 !important;
}</style><table data-quarto-disable-processing='true' class='cl-4e1eeec2'><thead><tr style="overflow-wrap:break-word;"><th class="cl-4e1a2f86"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">Parameter</span></p></th><th class="cl-4e1a2f90"><p class="cl-4e1a1eec"><span class="cl-4e02b554">Value</span></p></th><th class="cl-4e1a2f90"><p class="cl-4e1a1eec"><span class="cl-4e02b554">95% PI</span></p></th><th class="cl-4e1a2f90"><p class="cl-4e1a1eec"><span class="cl-4e02b554">05% PI</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">α</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">3.79</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">4.81</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">3.00</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">β</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.25</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.71</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">-0.19</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">δ</span><span class="cl-4e02b55e">2</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">1.90</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">2.42</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">1.42</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">δ</span><span class="cl-4e02b55e">3</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">2.02</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">2.46</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">1.50</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">δ</span><span class="cl-4e02b55e">4</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">1.83</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">2.41</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">1.39</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">ε</span><span class="cl-4e02b55e">2</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">-0.24</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">-0.03</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">-0.43</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">ε</span><span class="cl-4e02b55e">3</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">-0.06</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.18</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">-0.31</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">ε</span><span class="cl-4e02b55e">4</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.29</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.61</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">-0.03</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">γ</span><span class="cl-4e02b55e">2</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.50</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.88</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.11</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">γ</span><span class="cl-4e02b55e">3</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.44</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.81</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.10</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">γ</span><span class="cl-4e02b55e">4</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.62</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.94</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.14</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">Z</span><span class="cl-4e02b55e">1</span><span class="cl-4e02b55f">0</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">1.06</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">2.00</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.55</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">Z</span><span class="cl-4e02b55e">2</span><span class="cl-4e02b55f">0</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.83</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">1.17</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.46</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">Z</span><span class="cl-4e02b55e">3</span><span class="cl-4e02b55f">0</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">1.77</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">2.37</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">1.18</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">σ</span><span class="cl-4e02b55f">rec</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.18</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.30</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.11</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">σ</span><span class="cl-4e02b55f">Z</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.56</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.89</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.22</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f91"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">τ</span><span class="cl-4e02b55f">obs</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.94</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">1.26</span></p></td><td class="cl-4e1a2f9a"><p class="cl-4e1a1eec"><span class="cl-4e02b554">0.72</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e1a2f9b"><p class="cl-4e1a1ee2"><span class="cl-4e02b554">τ</span><span class="cl-4e02b55f">proc</span></p></td><td class="cl-4e1a2fa4"><p class="cl-4e1a1eec"><span class="cl-4e02b554">14.94</span></p></td><td class="cl-4e1a2fa4"><p class="cl-4e1a1eec"><span class="cl-4e02b554">128.42</span></p></td><td class="cl-4e1a2fa4"><p class="cl-4e1a1eec"><span class="cl-4e02b554">4.57</span></p></td></tr></tbody></table></div>

##### Statistics {-}

<div class="tabwid"><style>.cl-4e31d6b8{}.cl-4e2a7242{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-4e2d55a2{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4e2d55b6{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4e2d55b7{margin:0;text-align:center;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4e2d6682{width:0.2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4e2d668c{width:2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4e2d6696{width:0.2in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4e2d6697{width:2in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4e2d66a0{width:0.2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4e2d66aa{width:2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.tabwid {
  font-size: initial;
  padding-bottom: 1em;
}

.tabwid table{
  border-spacing:0px !important;
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  border-color: transparent;
  caption-side: top;
}
.tabwid-caption-bottom table{
  caption-side: bottom;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td, .tabwid th {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
.katex-display {
    margin: 0 0 !important;
}</style><table data-quarto-disable-processing='true' class='cl-4e31d6b8'><thead><tr style="overflow-wrap:break-word;"><th class="cl-4e2d6682"><p class="cl-4e2d55a2"><span class="cl-4e2a7242">Metric</span></p></th><th class="cl-4e2d668c"><p class="cl-4e2d55b6"><span class="cl-4e2a7242">Mean (900 samples)</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-4e2d6696"><p class="cl-4e2d55a2"><span class="cl-4e2a7242">DIC</span></p></td><td class="cl-4e2d6697"><p class="cl-4e2d55b7"><span class="cl-4e2a7242">339.9</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e2d6696"><p class="cl-4e2d55a2"><span class="cl-4e2a7242">pD</span></p></td><td class="cl-4e2d6697"><p class="cl-4e2d55b7"><span class="cl-4e2a7242">88.7</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e2d6696"><p class="cl-4e2d55a2"><span class="cl-4e2a7242">WAIC</span></p></td><td class="cl-4e2d6697"><p class="cl-4e2d55b7"><span class="cl-4e2a7242">274.9</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e2d6696"><p class="cl-4e2d55a2"><span class="cl-4e2a7242">deviance</span></p></td><td class="cl-4e2d6697"><p class="cl-4e2d55b7"><span class="cl-4e2a7242">251.2</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4e2d66a0"><p class="cl-4e2d55a2"><span class="cl-4e2a7242">pWAIC</span></p></td><td class="cl-4e2d66aa"><p class="cl-4e2d55b7"><span class="cl-4e2a7242">23.7</span></p></td></tr></tbody></table></div>

#### Posterior Samples & Densities {- .tabset}

##### Alpha & Beta {-}

Alpha ($\alpha$) is the estimate of the constant effect of larval density on recruitment, and beta ($\beta$) is the estimate of the linear effect of larval desnity on recruitment.

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

##### Delta {-}

Delta ($\delta$) is the estimate of the ice-width effect on survivability, i.e. total mortality, for each age.

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

##### Gamma {-}

Gamma ($\gamma$) is the estimate of the timing effect on survivabilty, i.e. total mortality, for each age.

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

##### Epsilon {-}

Epsilon ($\varepsilon$) is the estimate of the fish condition effect on survivability, i.e. total mortality, for each age. This affect is only present from 

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-2.png)![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-3.png)

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

##### Z0 {-}

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png)

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)

##### Variances {-}

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18-1.png)

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)


![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21-1.png)

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.png)

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23-1.png)

#### Correlations {- .tabset}

##### Correlation of Abundance-at-age {-}


```
## plot: [1, 1] [========>--------------------------------------------------------------------------------------------------------------------------------------] 6% est: 0s
## plot: [1, 2] [=================>-----------------------------------------------------------------------------------------------------------------------------] 12% est: 1s
## plot: [1, 3] [==========================>--------------------------------------------------------------------------------------------------------------------] 19% est: 1s
## plot: [1, 4] [===================================>-----------------------------------------------------------------------------------------------------------] 25% est: 1s
## plot: [2, 1] [============================================>--------------------------------------------------------------------------------------------------] 31% est: 1s
## plot: [2, 2] [=====================================================>-----------------------------------------------------------------------------------------] 38% est: 1s
## plot: [2, 3] [==============================================================>--------------------------------------------------------------------------------] 44% est: 0s
## plot: [2, 4] [=======================================================================>-----------------------------------------------------------------------] 50% est: 0s
## plot: [3, 1] [===============================================================================>---------------------------------------------------------------] 56% est: 0s
## plot: [3, 2] [========================================================================================>------------------------------------------------------] 62% est: 0s
## plot: [3, 3] [=================================================================================================>---------------------------------------------] 69% est: 0s
## plot: [3, 4] [==========================================================================================================>------------------------------------] 75% est: 0s
## plot: [4, 1] [===================================================================================================================>---------------------------] 81% est: 0s
## plot: [4, 2] [============================================================================================================================>------------------] 88% est: 0s
## plot: [4, 3] [=====================================================================================================================================>---------] 94% est: 0s
## plot: [4, 4] [===============================================================================================================================================]100% est: 0s
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24-1.png)

##### OSA Correlations {-}

![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25-1.png)

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26-1.png)

##### POSA Correlations {-}

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-27-1.png)

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-28-1.png)

#### Residuals {- .tabset}

##### Standard Residuals {-}

![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-29-1.png)

##### One-Step Ahead Residuals (OSA) by Cohort {-}


```
## Error in knitr::include_graphics(here("../output/osa_resid_cohort.png")): Cannot find the file(s): "../output/osa_resid_cohort.png"
```

##### One-Step Ahead Residuals (OSA) by Year {-}


```
## Error in knitr::include_graphics("../output/osa_resid_year.png"): Cannot find the file(s): "../output/osa_resid_year.png"
```

##### Pearson One-Step Ahead Residuals (POSA) by Cohort {-}

![plot of chunk unnamed-chunk-32](output/posa_resid_cohort.png)

##### Pearson One-Step Ahead Residuals (POSA) by Year {-}

![plot of chunk unnamed-chunk-33](output/posa_resid_year.png)


## Version 1.1 (2024-05-22) {-}

### Changes {-}

Three primary changes were made for this version.

First, I added a separate process error term to the recruitment to distinguish from the process error added to abundance, since abundance is fit to indices but recruitment is estimated directly and will likely have greater variance. Initially, it seemed that this is redundant (and likely nonidentifiable) because "processed" recruitment is a prior based on another prior, where

$$ \rho_y \sim dlnorm(rec_y, \sigma_{procrec}), $$
which is equivalent to

$$ \rho_y \sim dlnorm(dlnorm(rec_{y-1}, \sigma_{rec}), \sigma_{procrec}). $$

However, recruitment process error variance is markedly larger than all other variance estimates, and trends in recruitment changed significantly from the previous version (see below).

The second change was the inclusion of an estimated fish condition effect prior to and after data availability for condition. For the first 7 year of the assessment, fish condition, $CO_y$, follows a random walk forward in time,

$$ CO_y \sim dnorm(CO_{y-1}, \sigma_{CO}), $$

where 

$$ \sigma_{CO} \sim dlnorm(0, 100), $$
to keep the the variance near 1 (condition is normalized such that $\sigma \approx 1$, and a log-normal prior will ensure a positive value), and conditions start from the prior

$$ CO_{y=1} \sim dnorm(0, 1/1.44). $$
which is defined using the mean and inverse-square variance of the condition data from 1995-2021. Total mortality at age is now defined as

$$ 4 \gamma^a \frac{TI_y}{\delta^a} \left(1 - \frac{TI_y}{\delta^a}\right) + \varepsilon^a CO_y + \hat{Z}^a, $$
for all years, whre years prior to 1995 and afer 2021 use the estimated condition index. Note our definition for the additive mortality component changes.

The last and most drastic change was the change in total mortality estimates. Here, I changed the constant total mortality, $Z_0$, to a scaled, allometric mortality based on the *Lorenzen* equation for natural mortality given length [see @lorenzen2022a]. Here, I defined an allometric equivalent for *total* mortality as

$$ \hat{Z}^a = Z_{\infty}\left(\frac{L^a}{L_{\infty}}\right)^c,$$
where $Z_{\infty}$ is a parameter with an uninformed prior,

$$ Z_{\infty} \sim dunif(0.01, 5), $$
assuming fish mortality caps at 5 (which is reasonable giving previous estimates for the additive mortality). Also, $c = -1$ is a good approximation and tends to be a consistent scale across relationships between $lnM$ and $ln(L/L_{\infty})$, although without modelling random effects (which here I do not) $c$ may be less negative [see Table 2, @lorenzen2022b]. Here we assume mortality has no year covariate.

Our model for allometric mortality is based off the vonBertalanffy growth equation,

$$ L^a = L_{\infty}(1-exp(-K * a)).$$
Values for $K$ and $L_{\infty}$ are estimated by fitting length-age data to the above equation. Length-age data are stratified *by age* using **Proportionate Allocation** (see [Stratified Sampling](https://en.wikipedia.org/wiki/Stratified_sampling) and also [Proportional Allocation](https://www.sciencedirect.com/topics/mathematics/proportional-allocation)) using ages 1 to 4 (the max age available in the data) and a sample size of 1000. Length-age pairs were fit to estimate average length-at-ages for the data, such that

$$ length(a_i) \sim dnorm(L^a, \sigma_{L^a}) $$
where *i* is the index of the length-age sample. Several fits were made using ages 1 to 4 through 6 (always using 1000 samples), but RMSE was slightly smaller when ages were truncated at 4 likely due to the paucity of data for ages 5 and 6.

Additional changes were made to prior information. Priors on log abundance at age for the first year were changed to

\begin{equation}
   rec_{y=1} \sim dnorm(6.6, 1/3)
\end{equation}
\begin{equation}
   N^{a=2}_{y=1} \sim dnorm(5.9, 1/3)
\end{equation}
\begin{equation}
   N^{a=3}_{y=1} \sim dnorm(4.4, 1/3)
\end{equation}
\begin{equation}
   N^{a=4}_{y=1} \sim dnorm(1.3, 1/3).
\end{equation}

Previous values were set according to the abundance in billions \emdash these changes were made to reflect the abundance in millions.


### Diagnostics {- .tabset}



#### Stock Trends {- .tabset}

##### Index Estimates {-}

![plot of chunk unnamed-chunk-34](figure/unnamed-chunk-34-1.png)

##### Recruitment Estimates {-}

![plot of chunk unnamed-chunk-35](figure/unnamed-chunk-35-1.png)

##### Length-at-Age Estimates {-}

![plot of chunk unnamed-chunk-36](figure/unnamed-chunk-36-1.png)

##### Condition {-}

![plot of chunk unnamed-chunk-37](figure/unnamed-chunk-37-1.png)

#### Statistics {-}

<div class="tabwid"><style>.cl-55422b9c{}.cl-5537a85c{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-553b0e2a{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-553b0e34{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-553b0e3e{margin:0;text-align:center;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-553b1ef6{width:0.2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-553b1f00{width:2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-553b1f01{width:0.2in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-553b1f0a{width:2in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-553b1f0b{width:0.2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-553b1f14{width:2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.tabwid {
  font-size: initial;
  padding-bottom: 1em;
}

.tabwid table{
  border-spacing:0px !important;
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  border-color: transparent;
  caption-side: top;
}
.tabwid-caption-bottom table{
  caption-side: bottom;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td, .tabwid th {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
.katex-display {
    margin: 0 0 !important;
}</style><table data-quarto-disable-processing='true' class='cl-55422b9c'><thead><tr style="overflow-wrap:break-word;"><th class="cl-553b1ef6"><p class="cl-553b0e2a"><span class="cl-5537a85c">Metric</span></p></th><th class="cl-553b1f00"><p class="cl-553b0e34"><span class="cl-5537a85c">Mean (900 samples)</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-553b1f01"><p class="cl-553b0e2a"><span class="cl-5537a85c">deviance</span></p></td><td class="cl-553b1f0a"><p class="cl-553b0e3e"><span class="cl-5537a85c">8,251.9</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-553b1f01"><p class="cl-553b0e2a"><span class="cl-5537a85c">pD</span></p></td><td class="cl-553b1f0a"><p class="cl-553b0e3e"><span class="cl-5537a85c">158.0</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-553b1f01"><p class="cl-553b0e2a"><span class="cl-5537a85c">DIC</span></p></td><td class="cl-553b1f0a"><p class="cl-553b0e3e"><span class="cl-5537a85c">8,410.0</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-553b1f01"><p class="cl-553b0e2a"><span class="cl-5537a85c">pWAIC</span></p></td><td class="cl-553b1f0a"><p class="cl-553b0e3e"><span class="cl-5537a85c">30.1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-553b1f0b"><p class="cl-553b0e2a"><span class="cl-5537a85c">WAIC</span></p></td><td class="cl-553b1f14"><p class="cl-553b0e3e"><span class="cl-5537a85c">8,282.0</span></p></td></tr></tbody></table></div>

#### Posterior Samples & Densities {- .tabset}

##### Z {-}

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-39-1.png)

![plot of chunk unnamed-chunk-40](figure/unnamed-chunk-40-1.png)

##### Lengths {-}

![plot of chunk unnamed-chunk-41](figure/unnamed-chunk-41-1.png)

![plot of chunk unnamed-chunk-42](figure/unnamed-chunk-42-1.png)

##### Variances {-}

![plot of chunk unnamed-chunk-43](figure/unnamed-chunk-43-1.png)

![plot of chunk unnamed-chunk-44](figure/unnamed-chunk-44-1.png)

![plot of chunk unnamed-chunk-45](figure/unnamed-chunk-45-1.png)

![plot of chunk unnamed-chunk-46](figure/unnamed-chunk-46-1.png)

![plot of chunk unnamed-chunk-47](figure/unnamed-chunk-47-1.png)

## Version 1.2 (2024-05-22) {-}

### Changes {-}

I first attempted a random walk backward projection from condition data from 1995. This was an attempt to avoid priors on condition for 1985 and keep condition consistent between 1995 and 1994 (see condition trends for version 1.1). However, JAGS does not seem to allow random walks on data. Instead, I tried a simple normal prior on codition c=given the mean and variance of condition data\emdash which is standardized,

$$ eCO_y \sim dnorm(1, 1/1.44), $$

### Diagnostics {- .tabset}



#### Stock Trends {- .tabset}

##### Index Estimates {-}

![plot of chunk unnamed-chunk-48](figure/unnamed-chunk-48-1.png)

##### Recruitment Estimates {-}

![plot of chunk unnamed-chunk-49](figure/unnamed-chunk-49-1.png)

##### Condition {-}

![plot of chunk unnamed-chunk-50](figure/unnamed-chunk-50-1.png)

#### Statistics {-}

<div class="tabwid"><style>.cl-58389764{}.cl-5829dd00{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-582e937c{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-582e9390{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-582e9391{margin:0;text-align:center;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-582eabb4{width:0.2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-582eabbe{width:2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 1.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-582eabbf{width:0.2in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-582eabc8{width:2in;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-582eabc9{width:0.2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-582eabd2{width:2in;background-color:transparent;vertical-align: middle;border-bottom: 1.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.tabwid {
  font-size: initial;
  padding-bottom: 1em;
}

.tabwid table{
  border-spacing:0px !important;
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  border-color: transparent;
  caption-side: top;
}
.tabwid-caption-bottom table{
  caption-side: bottom;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td, .tabwid th {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
.katex-display {
    margin: 0 0 !important;
}</style><table data-quarto-disable-processing='true' class='cl-58389764'><thead><tr style="overflow-wrap:break-word;"><th class="cl-582eabb4"><p class="cl-582e937c"><span class="cl-5829dd00">Metric</span></p></th><th class="cl-582eabbe"><p class="cl-582e9390"><span class="cl-5829dd00">Mean (900 samples)</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-582eabbf"><p class="cl-582e937c"><span class="cl-5829dd00">deviance</span></p></td><td class="cl-582eabc8"><p class="cl-582e9391"><span class="cl-5829dd00">8,319.7</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-582eabbf"><p class="cl-582e937c"><span class="cl-5829dd00">pD</span></p></td><td class="cl-582eabc8"><p class="cl-582e9391"><span class="cl-5829dd00">161.7</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-582eabbf"><p class="cl-582e937c"><span class="cl-5829dd00">DIC</span></p></td><td class="cl-582eabc8"><p class="cl-582e9391"><span class="cl-5829dd00">8,482.1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-582eabbf"><p class="cl-582e937c"><span class="cl-5829dd00">pWAIC</span></p></td><td class="cl-582eabc8"><p class="cl-582e9391"><span class="cl-5829dd00">32.7</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-582eabc9"><p class="cl-582e937c"><span class="cl-5829dd00">WAIC</span></p></td><td class="cl-582eabd2"><p class="cl-582e9391"><span class="cl-5829dd00">8,352.5</span></p></td></tr></tbody></table></div>

# References
