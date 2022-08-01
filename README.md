# wf-SCA
Latch Bio SDK workflow of the Statistical Coupling Analysis (SCA) algorithm.


Statistical Coupling Analysis (SCA) is a statistical method for defining and extracting the connective pattern of coevolving amino acids in a protein family.
    
Latch Bio link: [SCA workflow](https://console.latch.bio/workflows/65969/info).  
    
# Statistical Coupling Analysis (SCA)

This package provides an implementation of the inference pipeline of Statistical Coupling Analysis (SCA) by the Ranganathan lab. This is a statistical method that defines and extracts the connective pattern of coevolving amino acids in a protein family. This algorithm returns the statistical structure of the multiple sequence alignment (MSA), SCA conserviation and evolution, and identify sectors. For more details please check out the paper by the Ranganathan group: ["Evolution-Based Functional Decomposition of Proteins" (Olivier R. et al., 2016)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004817). The four pre-define protein families are great starting points: G-coupled receptor family, S1A serine protease family, Dihydrofolate reductase (DHFR) family, and Beta-Lactamase enzyme family. 
    
## Citations
If you use the code or data in this package, please cite:
```
    @article{{rivoire2016evolution},
        title={Evolution-based functional decomposition of proteins},
        author={Rivoire, Olivier and Reynolds, Kimberly A and Ranganathan, Rama},
        journal={PLoS computational biology},
        volume={12},
        number={6},
        pages={e1004817},
        year={2016},
        publisher={Public Library of Science San Francisco, CA USA}
    }
```
    
## Liscense
The original SCA source code is free software distributed under the BSD 3-clause license, please see the file LICENSE for details. Copyright (C) 2019 Olivier Rivoire, Rama Ranganathan, and Kimberly Reynolds.
    
