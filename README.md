PRactIP: Protein-RNA interaction prediction by integer programming
==================================================================

PRactIP can predict residue-base contacts between proteins and RNAs using sequence information and structural informations predicted from only sequences. The method can be applied to any protein-RNA pair, even when rich information such as 3D structure is not available.

Residue-base contact prediction is formalized as an integer programming problem. PRactIP predicts a residue-base contact map that maximizes a scoring function based on sequence-based features such as k-mer of sequences and predicted secondary structure. The scoring function is trained by a max-margin framework from known PRIs with 3D structures.


Requirements
------------
* C++11 compatible compiler
  (tested on Apple LLVM version 6.1.0 and GCC version 4.8.1)
* [GNU Linear Programming Kit](http://www.gnu.org/software/glpk/) (>=4.41),
  [Gurobi Optimizer](http://www.gurobi.com/) (>=2.0), or
  [ILOG CPLEX](http://http://www-01.ibm.com/software/integration/optimization/cplex/) (>=12.0)
* [CentroidFold](https://github.com/satoken/centroid-rna-package)
* [PSIPRED](http://bioinf.cs.ucl.ac.uk/psipred/)

Install
-------

For GLPK,

	./configure --with-glpk

For Gurobi,

	./configure --with-gurobi

For CPLEX,

	./configure --with-cplex

You may have to specify the include path and the library path by CPPFLAGS and LDFLAGS like

	env CPPFLAGS='-I/path/to/gurobi/include' LDFLAGS='-L/path/to/gurobi/lib' ../configure --with-gurobi

Then,

	make
	make install # optional

Usage
-----
### Prediction

Run the following three commands for an amino-acid sequence (pro.fa) and an RNA sequence (rna.fa) with FASTA format:

	runpsipredplus pro.fa
	centroid_fold rna.fa > rna.ss
	practip pro.fa rna.fa

The first two commands predict the secondary structures for the amino-acid sequence and the RNA sequence. Using the predicted secondary structures, PRactIP predicts residue-base contacts between the given two sequences as follows:

	>Score=8.74228
	134 29
	134 33
	134 82
	146 34
	146 58
	146 71
	171 19
	171 51
	171 82

The first line shows the score for the predicted contacts. The second line and later indicate residue-base contacts, where the fist column is the position of the residue in the amino-acid sequence (starting from 0) and the second is the position of the base in the RNA sequence.

### Training

You can predict residue-base contacts using the default parameters which we provided. However, you can also train parameters with your own dataset as the following:

	practip --train=output.param training.lst

`training.lst' contains known protein-RNA interactions with residue-base contacts for each line:

	pro0.fa rna0.fa pro0-rna0.dat
	pro1.fa rna1.fa pro1-rna1.dat
	...

The first column is the filename of FASTA formatted amino-acid sequences, the second column is the filename of FASTA formatted RNA sequences, and the last column is the filename of residue-base contacts between the corresponding amino-acid sequence and RNA sequence, whose format is the similar to the output of PRactIP prediction, that is, the fist column is the position of the residue in the amino-acid sequence (starting from 0) and the second is the position of the base in the RNA sequence.

`output.param` contains the trained parameters that can be used for the prediction as follows:

	practip --predict=output.param pro.fa rna.fa

References
----------
* Sato, K., Kashiwagi, S., Sakakibara, Y.: A max-margin model for predicting residue-base contacts in protein-RNA interactions, in press. [preprint](http://dx.doi.org/10.1101/022459)
