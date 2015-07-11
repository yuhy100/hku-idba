# What Is IDBA? #
A Practical Iterative De Bruijn Graph De Novo Assembler related to _**[Sequence assembly](http://en.wikipedia.org/wiki/Sequence_assembly)**_ problem in bioinfomatics.


---

# Current Release #

Please look at the downloads for current release.


---

# Project Description #

IDBA is a open source de novo assembler for next-generation short read sequences. It is fast, parallel and capable of assembling large scale genomic assembly such as human genome.


---

# Publications #
_**[Yu Peng, Henry Leung, S.M. Yiu, Francis Y.L. Chin. IDBA - A Practical Iterative de Bruijn Graph De Novo Assembler (accepted by RECOMB 2010)](http://hku-idba.googlecode.com/files/idba.pdf)**_

Abstract: The de Bruijn graph assembly approach breaks reads into k-mers before assembling them into contigs. The string graph approach forms contigs by connecting two reads with k or more overlapping nucleotides. Both approaches face the problem of false-positive vertices from erroneous reads, missing vertices due to non-uniform coverage and branching due to erroneous reads and repeat regions. A proper choice of k is crucial but for any single k there is always a trade-off: a small k favors the situation of erroneous reads and non-uniform coverage, and a large k favors short repeat regions. We propose an iterative de Bruijn graph approach iterating from small to large k capturing merits of all values in between. With real and simulated data, our IDBA algorithm is superior to all existing algorithms by constructing longer contigs with similar accuracy and using less memory. The running time of IDBA is comparable with existing algorithms.


---

# Contact #
If you find a bug, please report it _**[here](http://code.google.com/p/hku-idba/issues/list)**_.

E-mail: _**ypeng@cs.hku.hk**_

[more support infomation...](http://code.google.com/p/hku-idba/wiki/Help)

