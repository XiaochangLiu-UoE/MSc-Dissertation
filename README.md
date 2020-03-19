# MSc Dissertation at UoE

Main goal: Performance evaluation of two Virtual Screening (**VS**) tools (LIDAEUS and Autodock Vina). Based on the evaluation results, making attempts to improve the **VS** performance of LIDAEUS.

## A bit background

Virtual screening (**VS**) is a technique to select suitable molecules from a compound pool, which is widely used in the area of drug discovery and material science. According to the selection criteria, it can be divided into two groups, ligand-based VS (**LBVS**) and structure-based VS (**SBVS**). 

- **LBVS** uses features extracted from existing compounds, which are the chemical structures or pharmacophore models. 
- **SBVS** is to evaluate whether a compound is suitable in a constricted space (binding site). 

And this project was to evaluate the SBVS performance.

Strictly, the performance of an **SBVS** programme should be evaluated from three aspects, which are **docking power** (also called sampling power), **scoring power** and **ranking power** (also referred to screening power), respectively. 

**Docking power** and **scoring power** are the programme's capabilities to generate accurate ligand binding pose and corresponding affinity (normally, binding free energy), respectively. 

Generally, **ranking power** doesn't care about whether the programme can generate correct ligand poses or accurate affinities but only cares about the ranking orders of compounds resulting from the programmes (usually according to their affinities).

In this project, we only focused on the **ranking power** of the programmes.

## How to evaluate

### Benchmarking dataset and Metric

Certainly, the feasibility of evaluating whatever power is dependent on the availability of experiment-validated data. For **docking power** and **scoring power**, things are much straightforward because what we need are just the ligand-receptor co-crystal structures and the affinities of the corresponding ligands, of which the availability is quite high and the data are less biased. However, for **ranking power**, things become a bit more tricky. 

There are two ways to establish the benchmarking dataset for evaluating the **ranking power**. 

#### Ranking power or Scoring power

The first one is to build an ordered dataset consisting of compounds with known affinity. And the order should be ranked according to the affinities. So, the common metric for evaluation using such a benchmarking dataset is the **Pearson Correlation**, which means that we expect the programme could generate a compounds' list with similar order to the benchmarking dataset. But, in some cases, people used this method to evaluate the **scoring power** rather than the **rank power** (actually, it just depends on what perspective you take). 

The advantage of this method is apparent. That is the data are less biased, even though the affinity values could be in doubt. However, at least, they are all validated by experiments.

#### Ranking power or Screening power

The second method is to build a list of active compounds and inactive compounds (also referred to as decoys). We expect that the programme can generate a compounds' list, which ranks all the actives in front of the decoys regardless of the order of active compounds. Apparently, such an evaluation deems the **VS** programme as a binary classifier (with this regard, call it **screening power** evaluation is more suitable than calling it **ranking power** evaluation, but, of course, it is still evaluating the **ranking power**, partially). Therefore, people borrowed metrics from the area of classifier, which are the Receiver Operating Characteristic (**ROC**) curve and the Area Under the **ROC** curve (**AUROC**). And people tend to use **AUROC** and its variants because numbers are more intuitive than pictures. 

##### Artificial bias

The disadvantage of this method is even more obvious. First, how should we find and define the inactive compound? In some datasets (**DUD-E** and **DEKOIS** series), decoys are artificial, which might have similar physiochemical properties to actives but different 2D-topology, vice versa. Because of this, people are always asking whether these decoys are really decoys. Such an issue is common in data-driven drug discovery. That is the lack of true negative data. 

##### Good metrics, bad metircs

Second, does the **AUROC** really work properly or are these metrics good enough? The answer is, **sometimes**. Because, soon, people found that multiple shapes of the **ROC** curve can share one single **AUROC** value but one of these shapes represents a better **VS** performance than the others did. The underlying concern is the so-called **early recognition problem** in **VS**. To see whether a programme can address such a problem, people borrowed and developed more metrics, including partial AUC (**pAUC**), Robust Initial Enhancement (**RIE**), Boltzmann-Enhanced Discrimination of ROC (**BEDROC**), Enrichment Factor (**EF**), power metric, predictiveness curve, statistical analysis framework, etc (in fact, observing the slope of the **ROC** cure can do such work). Essentially and certainly, all these metrics have more or less issues, but people seemed to give up addressing this good/bad metric problem since 2015 (the year when pridictiveness curve was introduced, but it seems no research group used it).

So, **AUROC**, **EF** and **BEDROC** were used in this project, but ultimately, we abandoned the **BEDROC** results because we thought the hyperparameter alpha of **BEDROC** is too sensitive. And, honestly, the definition of the so-called **early recognition problem** is quite vague and, to some extent, unrealistic. How early is the "**early**"? **BEDROC** is good metric and many people use it, but, still, we want to evaluate the performance from two aspects, which are the overall classification performance (using **AUROC**) and the early recognition performance (using **EF**).

## What did I do









