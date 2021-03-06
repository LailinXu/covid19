# covid19

Simple fit to COVID-19 cases as a function of time.

## Prerequists 
[Python](https://www.python.org/).

[ROOT](https://root.cern.ch/) (an open-source data analysis package for HEP).
How to download and install: https://root.cern.ch/downloading-root


## How to run

```
python  fun_covid19.py
```

The input data is under `data` as `txt` files, which could be found online, like 
* https://www.worldometers.info/coronavirus/
* https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide


How it works:
* Read data and fill them into a histogram;
* Do a simple exponential fit;
* Draw histograms and fit functions.


