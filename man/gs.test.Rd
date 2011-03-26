\name{gs.test}
\alias{gs.test}
\title{Internal function for gene set testing.}
\usage{gs.test(A, X, Y, group, tests, permutation=FALSE, nrot=1000, design,
    allocation.matrix=FALSE, verbose=FALSE)}
\description{Internal function for gene set testing.}
\value{List of the following, for up- and for down-regulation: Matrix with testing results for every gene set in its rows and the applied gene set test in its columns.}
\references{Artmann, Stephan and Jung, Klaus and Bleckmann, Annalen and Beissbarth, Tim (submitted).
Detection of simultaneous group effects in microRNA expression and 
related functional gene sets.

Brunner, E. (2009) Repeated measures under non-sphericity.
Proceedings of the 6th St. Petersburg Workshop on Simulation,
605-609.

Jelle J. Goeman, Sara A. van de Geer, Floor de Kort, Hans C. van
Houwelingen (2004) A global test for groups of genes: testing
association with a clinical outcome. Bioinformatics 20, 93-99.

Jung, Klaus and Becker, Benjamin and Brunner, Edgar and Beissbarth, Tim (submitted).
Comparison of Global Tests for Functinoal Gene Sets in
Two-Group Designs and Selection of Potentially
Effect-causing Genes.

Majewski, IJ, Ritchie, ME, Phipson, B, Corbin, J, Pakusch, M,
Ebert, A, Busslinger, M, Koseki, H, Hu, Y, Smyth, GK, Alexander,
WS, Hilton, DJ, and Blewitt, ME (2010). Opposing roles of polycomb
repressive complexes in hematopoietic stem and progenitor cells.
_Blood_, published online 5 May 2010.

Mansmann, U. and Meister, R., 2005, Testing differential gene
expression in functional groups, _Methods Inf Med_ 44 (3).

Smyth, G. K. (2004). Linear models and empirical Bayes methods for
assessing differential expression in microarray experiments.
_Statistical Applications in Genetics and Molecular Biology_,
Volume *3*, Article 3.

Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader,
JE, and Smyth, GK (2010). ROAST: rotation gene set tests for
complex microarray experiments. _Bioinformatics_, published online
7 July 2010.}
\author{Stephan Artmann}
\arguments{\item{A}{Allocation matrix as in "miR.test" function.}
\item{X}{miRNA expression matrix as in `miR.test' function. Only necessary when allocation.matrix=TRUE.}
\item{Y}{mRNA expression matrix as in "miR.test" function.}
\item{group}{group as in `miR.test' function}
\item{tests}{Test applied, sie gene.set.tests}
\item{permutation}{Shall permutation procedure for global tests be applied? Put 'FALSE' to use approximate results or give a number for the number of permutations.}
\item{nrot}{Number of rotations of rotation tests. Defaults to 1000 to be able to show p-values as low as 10^-3.}
\item{design}{If specified, group will be ignored. Design matrix as used in `limma' package. Cannot be used with global tests.}
\item{allocation.matrix}{Logical, is A an allocation matrix with mRNAs in its columns and miRNAs in its rows, or is it an allocation data.frame?}
\item{verbose}{Defaults to FALSE. If TRUE, progress is printed.}}
