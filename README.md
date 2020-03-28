# MSc Dissertation at UoE

Details are in the complete dissertation file. Examiners' feedback is also provided. Sample codes includes some scripts about a re-scoring function and a binding site classification method.
***
<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [A bit background](#a-bit-background)
- [How to evaluate](#how-to-evaluate)
	- [Benchmarking dataset and Metric](#benchmarking-dataset-and-metric)
		- [Ranking power or Scoring power](#ranking-power-or-scoring-power)
		- [Ranking power or Screening power](#ranking-power-or-screening-power)
			- [Artificial bias](#artificial-bias)
			- [Good metrics or bad metrics](#good-metrics-or-bad-metrics)
- [What did I do](#what-did-i-do)
	- [Manipulating the programme](#manipulating-the-programme)
	- [Weighting energy terms](#weighting-energy-terms)
	- [Re-scoring](#re-scoring)
	- [Binding site classification](#binding-site-classification)

<!-- /TOC -->
***
**Main goal**: Performance evaluation of two Virtual Screening (**VS**) tools (LIDAEUS and AutoDock Vina). Based on the evaluation results, making attempts to improve the **VS** performance of LIDAEUS.

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

The disadvantage of this method is even more obvious. First, how should we find and define the inactive compound? In some datasets (**DUD-E** and **DEKOIS** series), decoys are artificially made, which might have similar physiochemical properties to actives but different 2D-topology, vice versa. Because of this, people are always asking whether these decoys are really decoys. Such an issue is common in data-driven drug discovery. That is the lack of true negative data.

##### Good metrics or bad metrics

Second, does the **AUROC** really work properly or are these metrics good enough? The answer is, **sometimes**. Because, soon, people found that multiple shapes of the **ROC** curve can share one single **AUROC** value but one of these shapes represents a better **VS** performance than the others did. The underlying concern is the so-called **early recognition problem** in **VS**. To see whether a programme can address such a problem, people borrowed and developed more metrics, including partial AUC (**pAUC**), Robust Initial Enhancement (**RIE**), Boltzmann-Enhanced Discrimination of ROC (**BEDROC**), Enrichment Factor (**EF**), power metric, predictiveness curve, statistical analysis framework, etc (in fact, observing the slope of the **ROC** cure can do such work). Essentially and certainly, all these metrics have more or less issues, but people seemed to give up addressing this good/bad metric problem since 2015 (the year when pridictiveness curve was introduced, but it seems no research group used it).

So, **AUROC**, **EF** and **BEDROC** were used in this project, but ultimately, we abandoned the **BEDROC** results because we thought the hyperparameter $\alpha$ of **BEDROC** is too sensitive. And, honestly, the definition of the so-called **early recognition problem** is quite vague and, to some extent, unrealistic. How early is the "**early**"? **BEDROC** is good metric and many people use it, but, still, we want to evaluate the performance from two aspects, which are the overall classification performance (using **AUROC**) and the early recognition performance (using **EF**).

## What did I do
Except for conducting iterative virtual screenings against 81 targets using LIDAEUS and Vina respectively, I had made following attempts.

### Manipulating the programme

**VS** programmes are designed and implemented by human beings, but their behaviour is not deterministic or perfectly interpretable even if it is a deterministic VS programme like LIDAEUS. Let alone Vina, which uses a  stochastic optimization approach for pose searching. And therefore, they are, to some extents, like black boxes.

So, what I had done to theses two programmes was to assess the impacts of their parameter variation. And, for **SBVS**, the pre-defined constricted space (the binding site) is significant. LIDAEUS uses sitepoints, derived from the natural ligand, to guide docking, and Vina uses a 3D box, of which the centre and the volume (X * Y * Z) are set by the users. Guessing from the perspective of the developers, although the space is constricted, the programmes are still wished to search with quite a few degrees of freedom, which means they might just want the these parameters to guide the search instead of restricting. But it seems things are a bit different.

For LIDAEUS, in the case of screening against retinoid X receptor (**RXR**), different parameter settings can result in different **AUROC**, ranging from 0.30 to 0.57.

For Vina, when screening against B-cell lymphoma 2 (**BCL2**), the resulting **AUROC** can vary from 0.29 to 0.69 if we increase the box's volume from 1194 to 28,714 Å<sup>3</sup> (near the max volume warned by the official document).

Some groups also conducted this **AUROC** analysis in the screening against some targets when there are some known actives in the compound pool to show that their settings making the programmes function well.

But, when facing an unknown system (a new target with no active ligands), how should we adjust our programmes to their best states? Hopefully, if we know some key amino acid residues, we can filter some compounds according to whether they contact the key residues. Or, if the target has some homologies and we happen to know these homologies' active ligands, we might use these actives as pseudo actives. However, it is not always the case. In the case of the SARS-CoV-2 virus, it is over 96% homologous with a virus originating from Rhinolophus bat while only 76% with the original SARS-CoV. People surely know some active ligands of the SARS-CoV's main proteinase (**Mpro**), but for the **Mpro** of the Rhinolophus bat virus, they might not do enough research. Fortunately, there are only 12 mutations in the **Mpro** of the SARS-CoV-2 compared to the SARS-CoV, and there is no mutation in the binding pocket.

Still, preliminary assessment is important

### Weighting energy terms

The scoring function (**SF**) of Vina is tuned using the PDBbind. So, I attempted to parameterized the **SF** of LIDAEUS, which consists of three terms, the van der Waals energy (***vdW***), the hydrogen bond donors' energy (***HBD***) and the hydrogen bond acceptors' energy (***HBA***). The free energy formula in the output looks like this,

$$Free\ Energy = a × vdW + b × HBD + c × HBA \quad (eq.\ 1)$$

where, the weights, **𝑎**, **𝑏**, **𝑐**, are normally all set to 1. What I did is to systematically assign different weights to the formula and to subsequently re-rank the compound according to the new energy score and calculate the **AUROC**. Each weight ranges from 0 to 20 with an increment of 1, resulting in 9261 (21 * 21 * 21) outcomes for each target. But we know this could generate redundant combinations where ratio of a, b and c are the same. So, such combinations would only be calculated once. That is to say, we have 7514 weight combinations. Codes are followed.

		def create_jobs(step, reduced=True):
		    jobs = []
		    ratio_set = set()
		    if reduced:
		        for j in range(step):
		            for k in range(step):
		                for l in range(step):
		                    ratio = cal_ratio(j, k, l)
		                    if ratio not in ratio_set:
		                        ratio_set.add(ratio)
		                        jobs.append((j, k, l))
		        return tuple(jobs)
		    else:
		        for j in range(step):
		            for k in range(step):
		                for l in range(step)
		                    jobs.append((j, k, l))
		        return tuple(jobs)

		def cal_ratio(j, k, l):
			if j == 0 and k ==0 and l == 0:
				return (0, 0, 0)
			else:
				total = j + k + l
				j_ratio = j/total
				k_ratio = k/total
				l_raito = l/total
				return (j_ratio, k_ratio, l_raito)

Given that there is another **SF** of LIDAEUS, a knowledge-based **SF** meauring whether the docked ligand contacts the key amino acid residues and generating a score ranging from 0 to 1.0, I decided to combine this **SF** and the force fied-based **SF**. The weighted free energy formula looks like,

$$Free\ Energy = PIP_{x} × (a × vdW + b × HBD + c × HBA) \quad (eq.\ 2)$$

where the **PIP** stands for *pose interaction profile* and 𝑥 represents four types of the **PIP** scores according to different scoring criteria. After obtaining the entropy score, I did the same thing,

$$Free\ Energy = a × (vdW + HBD + HBA) − (b × T \Delta S_{water\ loss} + c × T \Delta S_{side\ chain}) \quad (eq.\ 3)$$

$$Free\ Energy = PIP_{x} × weighted(Enthalpy + Entropy) \quad (eq.\ 4)$$

Not surprisingly, for most targets, the **AUROC** was increased. But, in some cases, the best **AUROC** could only be 0.50 and this is because all the weights are set to 0 so that every compound shares the same score, which was deemed as a random selection by the **AUROC** calculation algorithm.

### Re-scoring

In the previous part, I mentioned that I had calculated the entropy and from the **eq. 3**, we could know that I calculate the entropy of the water loss and the flexible amio acid side chain. This entropy calculation is a re-scoring function, implemented by counting, which is so straightforward (codes in the sample code folder). The protocol is followed.

1. assign fix number of water molecules to the target and also the ligand, if they carry hydrogen bond donors/acceptors.
2. find the contacted residues within 3.5 Å around the ligand
   - count the water molecues between them
   - markdown the contacted residues
3. calcultate the entropy

$$Entroy = Number_{water\ loss} × Entropy_{water\ molecue} + Entropy_{side\ chain} \quad (eq.\ 5)$$

where the entropy values of both a single water molecule and the conformationally changed side chains, are set accoding to the work of [Dunitz (1994)](https://doi.org/10.1126/science.264.5159.670) and [Doig and Sternberg (1995)]( https://doi.org/10.1002/pro.5560041101), respectively.

Given the simplicity of this method, simply adding these entropic values to the original enthalpy scores did not improve LIDAEUS's ranking power, but, after assigning different weights (**eq. 3**), the **AUROC** increased and the outcome was better than only weighting the enthalpy (**eq. 1**).

### Binding site classification

In this project, I used LIDAEUS to screen 81 targets and no consistent performance of LIDAEUS was observed over any target families (like kinases). And therefore, to find a target family that LIDAEUS could function with similar performance, I tried to classify the binding site of different targets to find some patterns. Here, only one of the classification methods will be introduced. But, no matter which methods I used, the first step for me, was still to extract the binding sites. Unlike [Liu and Altman (2011)](https://doi.org/10.1371/journal.pcbi.1002326) who extracted the residues in a 7.5 Å-radius spherical space, to make it simple, I only extracted the contacted residues within 3.5 Å around the original natural ligands.

And then, instead of superposing the sites and calcualting the RMSD or measuring their volume to compare, I used a set of geometric descriptors ([Ballester et al., 2009](https://doi.org/10.1016/j.jmgm.2009.01.001)), which was used for ligand-based screening. Here is the protocol.

1. find and calculate four points, namely,
   - the centroid (ctd),
   - the closest atom to the ctd (cst),
   - the farthest atom from the ctd (fct)
   - the farthest atom from the fct (ftf)
2. calculate the distribution of euclidean distances from other atoms to the above four points (4 distributions in total)
3. calculate the related mean, variance and the skew of each atomic distance distribution (12 descriptors)
4. calculate the distance between each site by using the euclidean distance, namely,

$$\sqrt{\sum_{i=1}^{12}{(X_i-Y_i)^2}} \quad (eq.\ 6)$$

5. calculate the similarity using the following formula

$$\frac{1}{1 + \sqrt{\sum_{i=1}^{12}{(X_i-Y_i)^2}}} \quad (eq.\ 7)$$

In an improved published version ([Shave et al, 2015](https://doi.org/10.1371/journal.pone.0116570)), the number of descriptors was increased to 48, whereas, in my implementation, I used 36 descriptors derived by calculating the above 12 coefficients on sets of hydrophobic atoms (**vdW**), hydrogen bond acceptors (**HBA**) and hydrogen bond donors (**HBD**). So, for each site, the descriptor vector looks like,

$$
\begin{bmatrix}
vdW_{01}&vdW_{02}&{\cdots}&{vdW_{12}}&
HBA_{13}&HBA_{14}&{\cdots}&{HBA_{24}}&
HBD_{25}&HBD_{26}&{\cdots}&{HBD_{36}}\\
\end{bmatrix}
$$

And the similarity would be

$$\frac {1}{1 + \sqrt{\sum_{i=1}^{36}{(X_i-Y_i)^2}}} \quad (eq.\ 8)$$
