### MTHFD1 links folate metabolism to BRD4-mediated transcriptional regulation

Sara Sdelci, André F. Rendeiro, Philipp Rathert, Gerald Hofstätter, Anna Ringler, Herwig P. Moll, Wanhui You, Kristaps Klavins, Bettina Gürtl, Matthias Farlik, Sandra Schick, Freya Klepsch, Matthew Oldach, Pisanu Buphamalai, Fiorella Schischlik, Peter Májek, Katja Parapatics, Christian Schmidl, Michael Schuster, Thomas Penz, Dennis L. Buckley, Otto Hudecz, Richard Imre, Robert Kralovics, Keiryn L. Bennett, Andre Müller, Karl Mechtler, Jörg Menche, James E. Bradner, Georg E. Winter, Emilio Casanova, Christoph Bock, Johannes Zuber, Stefan Kubicek

<!-- **Paper**: [http://dx.doi.org/](http://dx.doi.org/) -->

This repository contains scripts used in the analysis of the data in the paper.

<br>

#### Analysis

Here are a few steps needed to reproduce it (more than I'd want to, I admit):

1. Clone the repository: `git clone git@github.com:epigen/mthfd1.git`
2. Install required software for the analysis:`make requirements` or `pip install -r requirements.txt`

If you wish to reproduce the processing of the raw data, run these steps:

2. Download the data from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105786).
3. Prepare [Looper](https://github.com/pepkit/looper) configuration files similar to [these](metadata/project_config.yaml) that fit your local system.
4. Run samples through the pipeline: `make preprocessing` or `looper run metadata/project_config_file.yaml`
6. Run the analysis: `make analysis`

Additionaly, processed (bigWig and narrowPeak files together with a gene expression matrix) are available from [GEO with accession number GSE105786](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105786).
