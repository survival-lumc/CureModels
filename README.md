# Cure models in survival analysis: Issues with identifiability and competing risks

**Authors**: Hein Putter, Per Kragh Andersen

## Abstract

Cure models have become increasingly popular over the last couple of decades. An appealing element of cure models is the idea that part of the population is immune to the event of interest. This paper raises three objections to cure models. First, we argue that true cure cannot exist for humans if the event of interest includes death, as death is inevitable. If death is excluded, competing risks must be accounted for. In this case it is more natural to model the cumulative incidence of both the event of interest and competing events over time, avoiding unreliable estimates at infinity and using observed data. Second, cure models rely on identifying plateaus in survival curves at the tail, where data is sparse and unreliable. Third, the mixture cure model faces issues with practical identifiability of covariate effects in the incidence model (probability of cure), especially alongside proportional hazards assumptions in the latency model (survival model, conditionally on not being cured). We support our arguments with real and simulated data. Full data and code is available online.

## Usage 

Files for Section 2: "BMT-cure-short.R" for R code and "BMT-cure-short.pdf" for resulting pdf.

Files for Section 3: "Cure-fup.qmd" for Quarto markdown file, "cure.bib" for bibtex file and "Cure-fup.pdf" for resulting pdf.

File for Section 4: "Cure_model_simulation.R"
