# SNP-SELEX

This repo is designed to store key scripts used in the analysis of SNP-SELEX data. The codes were written by [André M. Ribeiro dos Santos](mailto:andremrsantos@gmail.com) and [Yunjiang Qiu](mailto:serein927@gmail.com). 

## PWM score calculation

Scripts used to calculate PWM score and ΔPWM score are in pwm folder. 

- **compute-ic.py** script to calculate information content of motifs.
- **pfm2pwm.py** script to transform pfm to pwm.
- **score.py** script to calculate ΔPWM score.
- **threshold.py** script to calculate thresholds for PWM score.

## deltaSVM model 

Snakemake pipelines used for deltaSVM model training and testing are included in deltasvm folder. To run it, you need to download fastq files from PRJEB9797 in ENA database. Scripts written by us are included and a few dependencies are needed to be obtained seperately.

- **gkmtrain and gkmpredict** are from [LS-GKM](https://github.com/Dongwon-Lee/lsgkm).
- **subsample** is part of [MEME](http://meme-suite.org/doc/fasta-subsample.html). 
- **deltasvm.pl** is included but it is copied from [deltaSVM](http://www.beerlab.org/deltasvm/).
