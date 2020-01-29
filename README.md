# Estimating viability and predicting yield in Dilution to Extinction Cultures.


Update to the Button et al. 1993 method for calculating viability of cultures.


## How to set up this codebase

**Step one:** Download this codebase to your machine

```
git clone https://github.com/thrash-lab/viability_test.git
```

**Step two:** Setup a conda environment to run it in

```
conda create -n viability python=3.6 pandas numpy tqdm matplotlib
```

**Step three:** Change to your conda environment
```
conda activate viability
```

**Step four:** Check the code runs
```
python viability_test/viability_test.py -h
```

## How to use this code for your Dilution to Extinction culturing

The pipeline contains two methods:

1. A method to estimate how many positive and pure wells you will get from a DTE experiment for a given taxon
2. A method to evaluate the viability of a culture based on the number of observed wells in a DTE experiment

### Estimating yield from a DTE experiment

Let's say you are going to run an experiment whereby you inoculate 500 wells with 2 cells per well from a mixed community. Your taxon of interest makes up 50% of your community (based on 16S rRNA amplicon data). You want to know how many wells of that particular taxon you are likely to observe, assuming the cells are 100% viable:

```
python viability_test/viability_test.py predict_wells -w 500 -i 2 -r 0.05 -v 1
```

Unless you specify `--threads` or `-p`, it will use all the processors on the available machine to peform its bootstrapping.

After a few seconds, you should see something like the following:
```
Predicting outcome of a DTE experiment using 4 processors and 9999 bootstraps

    You simulated a DTE where you inoculated 500 wells with an average of 2 cells per well.
    Your estimated relative abundance for your taxon of interest was:       50.000%
    Your estimated viability of your taxon of interest was:                 100.000%
    
    After 9999 simulations:
    The median number of positive (not taxon-specific) wells was:           432 (417-447, 95%CI)
    The median number of wells inoculated with one cell was:                135 (116-155, 95%CI)
    The median number of wells containing your taxon of interest was:       316 (295-337, 95%CI)
    The median number of pure wells for your taxon of interest was:         116 (98-135, 95%CI)
```

So, your experiment will yield between 98 and 135 pure wells for your taxon of interest.


### Evaluating viability from a DTE experiment
Let's say you have performed a DTE experiment and you have observed 5 positive wells for your taxon of interest, instead of the 98-135 wells you were expecting from above, assuming a viability of 100%. You can use the method `estimate_viability` to estimate what viability of your taxon would explain the much lower-than-expected yields:

(note, this function is doing a LOT of bootstrapping, so be patient, or throw it on a large box).

```
python viability_test/viability_test.py estimate_viability -w 500 -i 2 -o 5 -r 0.5
```

Your output should look something like this:
 ```
Estimating viability using 4 processors and 9999 bootstraps
Minimum and maximum values are 0.01 and 0.06, respectively - refining:
Testing values between 0.000 and 0.070 in increments of 0.001

        You simulated a DTE where you inoculated 500 wells with an average of 2 cells per well.
        A range of decreasing viability was tested using 9999 bootstraps per experiment.
        Your estimated relative abundance for your taxon of interest was:       50.000%
        You observed:                                                           5 wells of interest
        
        If viability of your taxon had been 100%, you would have expected:      116 (98-135, 95%CI)
        
        The viability estimates that explain your observed counts is between:   0.900% - 6.500%

```

## Managing Compute
If the number of wells is large AND the number of bootstraps is large, then the bootstrapping can take a long time to perform. If you find your machine is struggling, then you can simply reduce the number of bootstraps using the `-b` or `--n_bootstraps` parameters for either method. 999 bootstraps will get you pretty close to decent values.