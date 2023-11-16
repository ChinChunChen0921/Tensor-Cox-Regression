# Multi-omics integration analysis using tensor-based covariates in Cox regression modeling

## Installation and download files

1. Install [R](https://www.r-project.org/)
2. Please download the two R files (TensorCoxReg_Main.R, TensorCoxReg_Function.R) and the demo CRC data (crc_data.RData) into the same directory.

## Run

1. Change your directory in TensorCoxReg_Main.R to source TensorCoxReg_Function.R and load the crc_data.RData.

## Description about the CRC data.

1. The data contain information of 555 subjects.
2. The progression free survival time and status are the outcomes of the model. 
3. Each subject has three clinical covariates (age, gender, and stage) and 2-order tensor covariates of 9 genes across 3 omic platforms. 
4. Within each platform, the data are standardized to have mean 0 and standard deviation 1 across samples.

## Instruction about the main function TensorCox().

### Description

### Usage
  TensorCox()
### Arguments

### Value

* A

