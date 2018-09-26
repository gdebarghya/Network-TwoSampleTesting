# Network-TwoSampleTesting
MATLAB implementations for two sample hypothesis tests for networks

The codes implements the methods in the paper:
- D. Ghoshdastidar and U. von Luxburg. Practical Methods for Graph Two-Sample Testing. In *Advances in Neural Information Processing Systems (NIPS)*, 2018

Please cite this paper if you use these codes/results in your work.

## Copyright notice
Copyright(c) 2018 [Debarghya Ghoshdastidar](https://gdebarghya.github.io)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Contents of the repository
The repository contains Matlab implementation of all methods and experiments, and data files used by the codes.
The use of the Matlab codes and files are given below:

#### Codes for two-sample tests (refer to paper for names):
* `GraphChi2Test.m`: code for Asymp-Chi2
* `NormalityTest.m`: code for Asymp-Normal
* `TracyWidomTest.m`: code for Asymp-TW
* `LowRankTests.m`: includes codes for Boot-ASE and Boot-EPA
* `ShufflingTests.m`: includes codes for Boot-Spectral and Boot-Frobenius

#### Codes for experiments:
* `experiment1a_main.m`: experiment for Fig 1
* `experiment1b_main.m`: experiment for Fig 4
* `experiment1c_main.m`: experiment for Fig 3
* `experiment2a_main.m`: experiment for Fig 2
* `experiment2b_main.m`: experiment for Fig 5
* `experiment_seizure.m`: experiment on EEG data (Tab 1-4)
* `experiment_oregon1.m`: experiment on Oregon networks using Asymp-Normal (Tab 5, Fig 6)
* `experiment_oregon2.m`: experiment on Oregon networks using Asymp-TW (Fig 7)

#### Others codes and data files:
* `genSparseGraph.m`: generates random graphs for experiments and for LowRankTests
* `plot_all_results.m`: used for generating figures
* `experiment_oregon_preprocess.m`: used for preprocessing the data using SNAP. 
* `oregon_networks.mat`: pre-processed data along with communities for SNAP Oregon dataset
* `seizure.csv`: contains data for the Epileptic seizures
* `TW_beta1_CDF.mat`: table of values for Tracy-Widom distribution function
* `results`: folder where all results will be saved

