# Distance Metric Learning 
This repository contains the python and MATLAB scripts for Distance Metric Learning for different datasets namely\
1. The **Braille Character** Dataset.
2. The **HAR-72** Haptic Acceleration Response Dataset.
## Brief Description
This is the code for the work done on **Distance **Metric** Learning** for Data **Clustering**. The report can be found [here](https://meghbhalerao.github.io/pdfs/Megh-Bhalerao-IITB-Internship-Report.pdf)
## Steps to run the code
1. Download **HAR-72/Braille** datasets along with their **confusion matrices**. 
2. To run the code to get the learned distance metric for the **braille** character dataset, run `python ./scripts/estimate.py`. This will print out and save the learned distance metric.
3. To run the code for the **HAR-72** dataset: `cd matlab_scripts/`, followed by running either of the `main.m` files, as each of them contains different variations of experiments done on the HAR-72 dataset. 
## References
1. [Distance metric learning, with application to clustering with side-information](https://ai.stanford.edu/~ang/papers/nips02-metric.pdf). _Eric P. Xing, Andrew Y. Ng, Michael I. Jordan and Stuart Russell. NIPS 2002._ 
2. [Tactile recognition of raised characters: A parametric study](https://link.springer.com/article/10.3758/BF03329767). _Jack M Loomis. Bulletin of the Psychonomic Society 1985, 23(1), 18-20._ 
3. [Distance Metric Learning for Large Margin Nearest Neighbor Classification](https://www.jmlr.org/papers/volume10/weinberger09a/weinberger09a.pdf). _Kilian Q. Weinberger and Lawrence K. Saul. JMLR Vol 10, 2009._ 
