#HSM6AP

HSM6AP is a machine-learning framework to predict m6A based on Sequence and Gnomic coordinates.

# Dependency

-python3

#content

-our model: data of homo sapiens for model training
-gnomic coordinates feature: 63-dimension, besed on Gnomic coordinates. (using GeneDerivedFeature.R)
-sequence feature: based on sequence.
-result: including the result of each model and voting.

#usage

1.Please unzip pseInOne. Tar. Gz

predict m6A-containing sequences

The script Full_main.py is used to predict if  given sequence and gnomic coordinates contain m6A sites. 

2.gnomic coordinates

(Please change the path and name of the your test and the output file name is name + "split.txt" in GeneDerivedFeature.R)
```
Rscript GeneDerivedFeature.R 
```

The required arguments

-infa a fasta file for test samples

This script ouput the prediction scores for given sequences. 

3.```
python FullMain.py -infa input_fa
```




# HSm6AP

The HSM6AP consists of two models.
The first is a complete experimental model, which includes the experimental model trained by the horizontal combination of gene-derived features and sequence-derived features. The first model fully expresses the whole idea of the experiment.
The second is a model based on sequential feature training.The uploading of this model is mainly to take into account that some users need to predict methylation sites across species.

Since we include two training samples, they are full transcripts and tree grown mrnas (full repersents full transcript and mature represents mature.).

SingleNucleotideReselutionData.csv and independent.csv comes from WHISLE.


The download address of Pse-in-One is http://bioinformatics.hitsz.edu.cn/Pse-in-One2.0/citation/
