# het-cl-episdet

<p>
  <a href="https://doi.org/10.1007/978-3-030-57675-2_38" alt="Publication">
    <img src="https://img.shields.io/badge/DOI-10.1007%2F978--3--030--57675--2--38-blue.svg"/></a>
    
</p>

This repository contains the implementations of exhaustive epistasis detection for second and third-order interaction searches, targeting Intel CPUs, GPUs and CPU+GPU systems. It supports single-objective and multi-objective evaluations with K2 score and Mutual Entropy scoring functions. The CPU implementations are parallelized by using OpenMP, while GPU kernels are deployed with the OpenCL programming model.

## What is Epistasis Detection?

Epistasis detection is a computationally complex bioinformatics application with significant societal impact. It is used in the search of new correlations between genetic markers, such as single-nucleotide polymorphisms (SNPs), and phenotype (e.g. a particular disease state).
Finding new associations between genotype and phenotype can contribute to improved preventive care, personalized treatments and to the development of better drugs for more conditions.

## Setup

### Requirements

* OpenCL (version xxxxx or more recent)
* OpenMP
* Intel Compiler (Not tested with other compilers)

### Compilation

Compiling binary (`<binary_name>`) for performing 2-way and 3-way searches using K2 Bayesian scoring:

* 2-way:
```bash
$ cd <some_folder> && make
```
* 3-way:
```bash
$ cd <some_folder> && make
```

Compiling binary (`<binary_name>`) for performing 2-way and 3-way searches using Mutual Entropy scoring:

* 2-way:
```bash
$ cd <some_folder> && make
```
* 3-way:
```bash
$ cd <some_folder> && make
```

Compiling binary (`<binary_name>`) for performing 2-way and 3-way searches using multi-objective evaluation (K2 score and Mutual Entropy):

* 2-way:
```bash
$ cd <some_folder> && make
```
* 3-way:
```bash
$ cd <some_folder> && make
```

## Usage example

Running a 3-way search with a synthetic dataset with 1000 SNPs (166,167,000 triplets of SNPs to evaluate) and 4000 samples:

```bash
$ ./<some_binary> 4000 1000 
```

Running a 3-way search with a dataset in .txt format with 1000 SNPs (166,167,000 triplets of SNPs to evaluate) and 4000 samples (2000 controls and 2000 cases):

```bash
$ ./<some_binary> <some_dataset>
```
The application receives the input file in a binarized format. 

<!--
This example is expected to take slightly less than 2 minutes to execute and to achieve a performance of above 25 tera sets (triplets) of SNPs processed per second (scaled to sample size), when executed on a system with a GeForce 2070S GPU.
Higher performance can be achieved when processing more challenging datasets with more SNPs.

The construction of contingency tables, a phase of epistasis detection searching that counts of occurrences of the possible genotypes in cases and controls resulting from combining pairs/triplets of SNPs, represents the most computationaly complex portion of the application.
Thus, running the same example with Mutual Information instead of K2 Bayesian scoring is expected to achieve comparable performance.


Important parameters such as the number of SNPs per block (`BLOCK_SIZE`) and the number of CUDA streams used to process different rounds (`NUM_STREAMS`) can be changed by modifying the Makefile.
Depending on the dataset characteristics, specializing these parameters (e.g. using a larger block size when processing datasets with more SNPs) can have a significant influence on the performance achieved.

The Makefile is expecting the CUTLASS library to be inside the project root directory in a folder named `cutlass`.
If you installed the library in a different directory, you must modify the Makefile accordingly.

Notice that the application expects that the input dataset is in a particular binarized format.
You can download an example dataset with 4096 SNPs and 262144 samples from <a href="https://drive.google.com/file/d/1htjD1QCr5_LEPo3udQEJ-5XUX4TK65JM/view?usp=sharing">here</a>.
Due to the way data is processed using matrix operations, the number of bits per {SNP, allele} in the dataset files (\*.bin) representing cases or controls (stored in different files) must be a multiple of 1024 bits. In situations where the number of cases or controls is not a multiple of 1024, the input binary dataset is expected to be padded with zeros (i.e. unset bits). 
-->



## In papers and reports, please refer to this tool as follows

Campos R., Marques D., Santander-Jiménez S., Sousa L., Ilic A. (2020) Heterogeneous CPU+iGPU Processing for Efficient Epistasis Detection. In: Malawski M., Rzadca K. (eds) Euro-Par 2020: Parallel Processing. Euro-Par 2020. Lecture Notes in Computer Science, vol 12247. Springer, Cham. https://doi.org/10.1007/978-3-030-57675-2_38.

BibTeX:

    @InProceedings{10.1007/978-3-030-57675-2_38,
    author="Campos, Rafael
    and Marques, Diogo
    and Santander-Jim{\'e}nez, Sergio
    and Sousa, Leonel
    and Ilic, Aleksandar",
    editor="Malawski, Maciej
    and Rzadca, Krzysztof",
    title="Heterogeneous CPU+iGPU Processing for Efficient Epistasis Detection",
    booktitle="Euro-Par 2020: Parallel Processing",
    year="2020",
    publisher="Springer International Publishing",
    address="Cham",
    pages="613--628",
    isbn="978-3-030-57675-2"
}

<!--For additional readings in high-throughput epistasis detection, you can take a look at our IPDPS 2020 and JSSPP 2020 papers.-->


