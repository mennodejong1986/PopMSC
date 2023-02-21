# PopMSC
Pipeline for multispecies coalescent (MSC) based analyses on populations, including haploblock detection, testing the equality of alternative quartet frequencies (using the D-statistic), and 'ingroup-score' calculation. This workflow may serve to detect ancestral population structure, as explained below.  

## Using the D-statistic to infer ancestral population structure

Since it’s introduction in 2011, the Four Taxon Test (also known as the D-statistic, ABBA-BABA test or f4-score) has become a standard genomics tool to detect gene flow. Here, I describe how to Four Taxon Test can alternatively also be used to reconstruct past population structure.      

The procedure described here has been developed for better understanding of the mechanisms underlying single locus and multi-locus discordances, in particular those observed in brown bears (Ursus arctos). The procedure allows to test between two potential explanations of mitonuclear discordances. I will refer to these competing explanations as the ‘omniscient autosome’-hypothesis and the ‘short-sighted autosome’-hypothesis.  

## The ‘omniscient autosome’-hypothesis
One potential explanation for mitonuclear discordances is that discordant mtDNA signals are confounded by single-locus stochastics. Being comprised of thousands of independently segregating loci, nuclear genomic datasets (i.e., autosomal and/or X-chromosomal data) are not susceptible to these confounding effects, and hence clarify the confusion introduced by single-locus analyses. 

To understand the underlying reasoning, consider an ancestral population termed R (for root) containing two mtDNA haplotypes (I and II), which splits into daughter populations A and B, which next subdivide into granddaughter populations A1 and A2, and B1 and B2 (Fig. 1). Theory predicts that owing to the stochastic fluctuations of allele frequencies from generation to generation, eventually and inevitable each granddaughter population will retain one mtDNA haplotype only. If the second population split event (where A splits into A1 and A2, and B splits into B1 and B2) occurs soon after the first population splits (ancestral into A and B), there is limited opportunity for allele frequency changes in populations A and B, meaning there will be little difference in the probabilities which with each granddaughter population will retain either allele. One likely outcome is that three populations (say A1, A2, and B1) retain haplotype I, while a fourth population (say B2) retains haplotype II (Fig. 1). If constructing a phylogenetic tree from this data, the tree topology will suggest that population B1 is closely related to populations A1 and A2, and only distantly related to population B2 (Fig. 1).

This example illustrates why inferring the phylogeographic history from mtDNA, or any other single locus (e.g., Y-chromosomal data), can be misleading. If, instead, working with multiple unlinked loci, such as the autosomal dataset, the single-locus chance effects will cancel each other out. For some loci, populations will retain the same allele; for other loci, they will not. On balance, the proportions of shared alleles will accurately reflect the evolutionary relationships between populations. The law of large numbers tells us to rely on the autosomal and X-chromosomal data, and to be sceptical about discordant signals obtained from single-locus datasets (such as mtDNA and Y-chromosomal data).   

![alt text](https://github.com/mennodejong1986/PopMSC/blob/main/ILS_scheme.png)
***Figure 1. Differential fixation of ancestral variation can cause single-locus discordances.*** *Schematic visualisation depicting how differential fixation of ancestral variation, also known as incomplete lineage sorting, can cause single-locus trees, such as mtDNA phylogenies, which have a different topology compared to the true species tree (as indicated in grey). This effect will be cancelled out when studying multiple-locus discordances. Blue and red lines indicate two alleles of a single locus.* 

![alt text](https://github.com/mennodejong1986/PopMSC/blob/main/Haploblock_analyses.png)

## The ‘short-sighted autosome’-hypothesis
If we were to believe the ‘omniscient autosome’-hypothesis, the discordances between the single-locus and multi-locus datasets can be dismissed as single-locus stochastics. It implies that inferences from mtDNA and Y-chromosomal data are superfluous at best, and misleading at worst. Is this really the case? While random differential fixation of ancestral variation undoubtedly affects the distribution of mtDNA and Y-chromosomal haplotypes, the hypothesis does not preclude alternative mechanisms causing discordances. One indication that another mechanism must indeed be at play, are the observed inconsistencies between the autosomal and X-chromosomal data sets. As both are m ulti-locus datasets, single-locus stochastics are clearly not to blame. Which other factor, then, could possibly explain the observed discordances?    

The hypothesis we came up with was: partial or complete non-overlapping temporal resolution. The four genomic regions, including the multi-locus datasets, zoomed in on different time intervals of brown bear phylogeographic history. To be more precise, we hypothesised that the autosomal dataset reflects contemporary population structure, while the non-recombining mtDNA and Y-chromosomal datasets, and to a lesser extent the X-chromosomal dataset, are better in retaining signals of a former population structure. 

Consider, as an example, the brown bear populations in Kamchatka and southwestern Alaska, on opposite sides of the Bering Strait. Until approximately 11000 years ago, these populations were connected by the Beringian land bridge. This shared ancestry is still evident from the sex chromosomal and mtDNA data sets, which all identify Kamchatka bears and southwestern Alaska bears as a monophyletic clade. Autosomal data, in contrast, suggest that Kamchatka bears cluster with Russian Far East populations, and all Alaska bears with North American populations – consistent with present-day population connectivity .  

Thus, the autosomal dataset is not omniscient. Instead, it is more like a short-sighted man who, facing an object, can clearly see and describe the front sections but has troubles to discern anything beyond. This information on brown bear population structure in the distant past would be missing if it were not for three far-sighted men which examine and describe these remaining body parts: the mitogenome, the X chromosome and the Y chromosome. 

## Why is the autosome short-sighted?
The short-sightedness of the autosome is caused by gene flow, which occurs when unrelated individuals meet and produce offspring following the crossing or removal of a migration barrier. This process overwrites, or updates, genetic information in the genomes of the donor population, therewith obscuring evolutionary relationships. 

Consider, again, the ancestral population which splits into daughter populations A and B and subsequently into granddaughter populations A1, A2, B1 and B2. This time we assume that the consecutive population splits are separated by a sufficient amount of time, such that genetic drift causes large allele frequency differences between populations A and B. As a result, a dataset of just a few unlinked loci will suffice to correctly and confidently infer the evolutionary relationship between the granddaughter populations as ((A1,A2),(B1,B2)) (Fig. 2a). 

The phylogenetic analyses become more challenging when individuals from population A2 start to mate with individuals from population B1 (Fig. 2b). This genetic exchange will gradually diminish the genetic differences between populations A2 and B1, such that at some point populations A2 and population B1 will have become more similar to each other than to their respective sister populations. From this point onwards, the best fitting binary tree model is the topology ((A2,B1),(A1,B2)) (Fig. 2c). 

As illustrated by this simple example, a binary tree does not necessary depict evolutionary relationships between individuals or populations, but instead present-day genetic distances only. The topology ((A2,B1),(A1,B2)) could, for instance, also have arisen in the absence of gene flow, namely if A2 and B1 shared a recent common ancestor. 

In an ideal world in which gene flow does not exist, phylogenetic reconstruction would be relatively straightforward. In practice, the possibility of gene flow can rarely be ruled out, even when working with a cross-species dataset, let alone when analysing an intraspecies dataset. This causes severe headaches to evolutionary biologists, because when allowing for gene flow any set of observations is usually consistent with a myriad of potential demographic scenarios, meaning that phylogenetic reconstruction alone does not suffice to infer the past. 

Depending on the magnitude and the duration, gene flow will gradually erase genomic evidence of the former population structure. Theory predicts that a new migration-drift equilibrium is reached after approximately 1/m generations, with m denoting the proportion of the recipient population consisting of new immigrants each generation. However, the pendulum will have swung even long before the new equilibrium sets in. Somewhere halfway along the transitional process so many genomic regions have been affected by gene flow, that binary trees outputted by phylogenetic analyses will no longer reflect the former population structure, but instead the new population structure.      

![alt text](https://github.com/mennodejong1986/PopMSC/blob/main/Dstats_scheme.png)
***Figure 2. Why does the D-statistic sometimes reveal gene flow, and in other cases ancestral population structure?*** *Shown are the expected outcomes of three evolutionary scenarios for an imaginary dataset of 1000 unlinked loci. a. Prior to the onset of gene flow (‘no gene flow’-scenario), the alternative quartet topology frequencies q2 and q3 occur in equal frequencies, being represented by 200 loci each. Hence the D-statistic is zero. b. Owing to genetic exchange between populations A2 and B1, the topology q3, which is A1,B2|A2,B1 increased in frequency. The D-score of -0.263 indicates that the alternative topologies q2 and q3 do no longer occur in equal frequencies, and suggests gene flow either between A2 and B1, or between A1 and B2. c. When migration rates between populations A2 and B1 are high (or alternatively: if ongoing for many generations – not depicted here), the topology frequencies may change so much, that topology q3 overtakes topology q1, and therewith become the most frequent topology. Now that this tipping point has been reached, the consensus tree no longer indicates the past population structure, but instead the contemporary population structure. Furthermore, the D-statistic does no longer indicate gene flow, but instead the former population structure.*      

## Using the D-statistic to correct the short-sightedness of the autosome
The autosome is short-sighted, but this is a symptom which can be cured. We can apply a method to correct the short-sightedness of the autosome, such that it can disclose more information, including about events in the distant past. In our analogy, we can think of this approach as designing glasses for the short-sighted man, so that he cannot only see the closer front part of the object he is examining, but also the structures further to the back. This ‘correction method’ is the Four Taxon Test, also known as D-statistic, ABBA-BABA test or f4-score. The test is typically applied to detect gene flow, but can also serve to detect the genetic signals of population structure prior to the onset of gene flow.

To understand the workings of the Four Taxon Test, consider again the demographic scenario in which an ancestral population splits into populations A and B at time t1, which subsequently, at time t2, subdivide into granddaughter populations A1, A2, B1 and B2 (Fig. 1). Consider, furthermore, that we genotyped one individual per granddaughter population, for a total of 1000 autosomal, unlinked loci. If the consecutive split events t1 and t2 were separated by a sufficient amount of time, the ancestral populations A and B will have retained only one ancestral allele per locus by the time of the second population split (t2). In that case, all 1000 locus trees will agree with the true topology: ((A1,A2),(B1,B2)), which we denote from here onwards as q1. (If populations A and B retained the same allele, these alleles will have diverged due to novel mutations, which allows to group A1 with A2, and B1 with B2).   

If, instead, the two consecutive split events t1 and t2 occur in quick succession, the populations A and B do not have time to retain only one ancestral allele per locus. A geneticist might say that lineage sorting is incomplete and that both populations are not reciprocally monophyletic. As a result, and as discussed already above, the daughter populations A1 and A2 – and similarly B1 and B2 – may inherit multiple alleles per loci (e.g., both allele I and allele II), and by chance eventually retain a different one. For example, populations A1 and B1 might retain allele I, while populations A2 and B2 happen to retain allele II. If so, the tree topology suggested by this locus is not the true topology, but instead: ((A1,B1),(A2,B2)), here denoted as q2. Other loci might suggest, for the same reason, the third potential topology, namely ((A1,B2),(A2,B1)), here denoted as q3 (Fig. 2).

Theory predicts that as long as the two split events t1 and t2 did not occur simultaneously, most loci will agree with the true topology (q1). With increasing time between t1 and t2, this proportion approaches 1. With decreasing time between t1 and t2, this proportion drops and converges to 1/3, but will never below 1/3rd. This majority rule allows to infer the true topology simply by determining which of the three topologies is most frequent. Theory also predicts that, because the retention of ancestral alleles by genetic drift is a random process, the discordant topologies q2 and q3 are equally likely to occur and will therefore appear in roughly equal frequencies. For instance, let us assume that for the above scenario, 600 out of 1000 loci indicate the true topology q1, while the alternative topologies q2 and q3 are represented by 200 loci each (Fig. 2a). By applying the majority rule, we can infer the species tree to be ((A1,A2),(B1,B2)). 

One special feature of these topology frequencies is that they are not affected by random allele frequency changes and hence do not change over time. That is, as long as there is no gene flow between populations. When gene flow does occur, the topology frequencies will change correspondingly. For instance, if individuals in population B1 receives immigrants from population A2, the frequency of loci with topology q3 – that is: ((A1,B2),(A2,B1)) – will increase, while the number of loci with topologies q1 and q2 will decrease. One potential intermediate outcome is that topologies q1, q2 and q3 are represented by 525, 175 and 300 loci respectively (Fig. 2b). Topology q1 is still the dominant topology, but notice that the two alternative topologies q2 and q3 do no longer occur in equal frequencies. This imbalance between q2 and q3 can be quantified, using the D-statistic, as follows: (175–300)/(175+300) = -0.263 (Fig 2b). Values highly different from zero, as this value, are strong evidence for gene flow. The negative value furthermore indicates an excess of loci with q3 topologies, suggesting that the gene flow occurred either between populations A1 and B2, or alternatively between populations A2 and B1. 

What happens when the genetic exchange between populations A2 and B1 continues? There will come a moment that the frequency of loci with topology q3 will have increased so much, at the expense of the other two topologies, that topology q3 overtakes topology q1. From this point onwards, the majority rule will make us believe that topology q3 is the true topology (Fig. 2c). Indeed, it will be around this moment in time that the tree topology outputted by tree-building methods, which collapse the genetic information present throughout the genome, will no longer reflect the ancestral population structure, but instead the new population structure. 

The rise of topology q3 also radically alters the meaning of the D-statistic. With topology q3 now being regarded the species tree, the D-statistic will indicate that topology q1 and q2 occur in unequal frequencies. For example, if the topologies q1, q2 and q3 are now represented by 420, 140 and 440 loci respectively, the D-statistic becomes (420–140)/(420+140) = 0.5 (Fig. 2c). In the traditional usage of the Four Taxon Test, this outcome would be regarded as strong evidence for gene flow, either between populations A1 and A2, or alternatively between populations B1 and B2. However, once the gene flow has caused a tipping-point where one of the two alternative quartet topologies has become the most frequent topology, this interpretation needs revision. The second most frequent topology is no longer the alternative topology favoured by gene flow, but instead the topology reflecting the former population structure – the situation prior to the onset of gene flow. 

Unfortunately, from the data itself it is impossible to determine whether the tipping-point has been reached or not. The generally valid interpretation is that a non-zero D-value is indicative of shared ancestry, without going into details on the exact underlying scenario.               


## Applying the test to the brown bear dataset
Back to the whole-genome dataset of brown bears. Because the Four Taxon Test can extract from genomic data information about events in the distant past, it allowed us to put our two competing hypotheses to the test. To recall, the ‘short-sighted autosome’-hypothesis posits that discordant geographical distributions of mtDNA and Y-chromosomal haplotypes are meaningful, lingering signals of past population structure. The ‘omniscient autosome hypothesis’, in contrast, posits that the observed discordances are nothing more than noise generated by single-locus stochastics. We applied the Four Taxon Test to determine for each separate case of brown bear single-locus discordance, which of the two hypotheses holds true. 

If mtDNA discordant signals arise from single locus-stochastics, the two alternative quartet topologies are expected to occur in equal frequencies across the genome. If, instead, the mtDNA discordant signals reflect a former population structure, an imbalance should exist between the two alternative quartet topology frequencies. These contrasting predictions provide the opportunity to formally test the validity of the two hypotheses by measuring, using the D-statistic, the equality of the two alternative quartet topology frequencies. 

To that end, we divided the brown bear genomes into many small, non-recombining regions, so called haploblocks. Because recombination occurs between these regions but not within these regions, these haploblocks are effectively independently evolving loci, just like the mitogenome and the Y-chromosome. We also phased the data, meaning that we identified two haplotypes for each diploid individual. In theory, these haploid, non-recombining haplotypes should behave very similar to the mitogenome and the Y-chromosome. But while the mitogenome and Y-chromosome are alone, the procedure of chopping the autosomal chromosomes into haploblocks yielded thousands of unlinked, non-recombining loci – a number suitable for statistical testing. To be more precise, we selected around 3000 loci, each of them containing sufficient genetic variation for robust phylogenetic inferences.               

Using this input dataset, we calculated the D-statistic for all possible population quartets among our brown bear individuals. As our individuals had been assigned to 31 different populations, we ran the analysis for 31465 quarters, which is the total number of possible combinations of four populations given a total dataset of 31 populations (i.e., 31 choose 4 = 31465). Consistent with the hypothesis of a shared ancestral Beringian population, the tests revealed inequality of alternative topology frequencies for quartets containing Kamchatka and western Alaskan bears. The tests also produced statistically significant evidence that Syrian bears share ancestry with European bears, even though at present they are genetically slightly more similar to Himalayan bears.

In contrast, no signals of shared ancestry (i.e., inequal alternative quartet topology frequencies) were observed for brown bears from southern Scandinavia and the Cantabrian Mountains in Spain, even though both populations carry mtDNA haplotype 1a. Traditionally, the presence of mtDNA haplotype 1a in these isolated populations has been considered evidence for postglacial recolonisation of southern Scandinavia from an Iberian refuge. The outcome of the Four Taxon Tests does not support this classic hypothesis. 

Likewise, no supporting evidence was found for the hypothesis of multiple migration waves onto the Japanese island Hokkaido. More specifically, Hokkaido bears which carry mtDNA haplotype 3a are not closer related to nearby mainland bears (which also carry haplotype 3a) than the other two populations on Hokkaido Island are (which carry mtDNA haplotypes 3b and 4).

## Other support for the ‘omniscient autosome’-hypothesis
The Four Taxon Test thus revealed multiple instances where brown bear mitonuclear discordances can be sufficiently explained with the ‘omniscient autosome’-hypothesis, according to which these discordances are simply single-locus stochastic noise produced by random differential fixation of ancestral variation. Other pieces of evidence also hint at that direction. One indication is that the divergence times of the different mtDNA haplotypes vastly predate the population split times suggested by the timing of geological events, such as the emergence of land bridge islands or the establishment of LGM glacial. Such deep split times are a typical signature of differential fixation of ancestral variation.  

Another indication is that the discontinuities between neighbouring mtDNA clades coincide with geographical boundaries which facilitate differential fixation. The mtDNA break in southern Scandinavia (between haplotypes 1a and 3a) occurs in an area of low bear population density, which separates two adjacent female concentration areas. This lowland extends for over 100 kilometres and is rarely crossed by females. Similarly, the occurrence of the three distinct mtDNA haplotypes (3a, 3b, 4) on Hokkaido Island coincides with the geographical ranges of three mountain areas, which are separated by unforested lowlands.

Although we assumed that the non-recombining autosomal haplotypes (which we used as input for the Four Taxon Tests) behave like mtDNA haplotypes, there is in fact a relevant difference. In diploid populations, both parents contribute two autosomal alleles to their offspring, adding up to a total of four alleles per parental cross. In contrast, in the case of mtDNA only one allele is transmitted by one of the parents only, adding up to one allele per parental cross. Thus, the effective population size (Ne) of mtDNA is only one-fourth that of autosomal DNA. This relatively low effective population size is yet another contributing factor to the random differential fixation of mtDNA haplotypes. The average time required for the loss of ancestral variation depends on the effective population size: the smaller the population, the less time needed. 


