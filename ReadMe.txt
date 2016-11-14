What is RBrel14.exe?

The program RBrel14.exe analyzes the spatial pattern of data that are in the form of spatially referenced counts. RBrel14.exe is a console-mode 32-Bit program designed to run under Windows 9x and newer. Input and output to RBrel14.exe is through simple ASCII data files.

RBrel14.exe is Copyright © 2008, Kelvin F. Conrad.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

You may contact Kelvin Conrad, the author of SADIEShell via e-mail at:
conradkf@hotmail.com


SADIE

SADIE is a acronym for Spatial Analysis by Distance IndicEs. The concepts underlying SADIE regard a set of data as represented by regions, within which the observed counts are either arranged effectively at random, or form local neighbourhoods of similarly-sized counts close to one another, termed clusters. In SADIE, spatial pattern is measured locally, at each sampled unit, through an index of clustering. Each unit with a count greater than the overall mean is assigned a patch cluster index, which by convention is positive. Each unit with count less than the overall mean is assigned a gap cluster index, which by convention is negative. Each index is computed to allow for the size of the count (abundance) at each sample unit.


RBrel14.exe

RBrel14.exe analyzes the spatial pattern of spatially referenced counts which are taken at specific locations such as numbers of moths in light-traps, numbers of plants in selected quadrats, where the two-dimensional location of the traps or quadrats is known. RBrel14.exe measures and detects the degree of clustering in the form of patches and gaps in the data.
The software uses the "red-blue techniques", described in: Perry, J.N., Winder, L., Holland, J.M. & Alston, R.D., (1999), Red-blue plots for detecting clusters in count data, Ecology Letters, 2, 106-113. The output of the software may be used in other graphics packages, such as Surfer® or Genstat®, to produce coloured graphical displays and maps of the clustering in the data.
In order to understand the output from the program it is essential to read Perry et al (1999).
Additionally, indices and randomization tests based on previous work are provided, specifically those based on the distance to regularity and the distance to crowding, as described in: Perry, J.N., (1998), Measures of spatial pattern for counts, Ecology, 79, 1008-1017.


Overview of the program

The program uses two files for input. If you aren't using a shell or the interactive mode, the first input file must have the name: RBni5.dat; this contains the spatial coordinates with the corresponding raw counts.
The second file must be called RBni8.dat; this contains two parameters that control program execution. Examples of files RBni5.dat and RBni8.dat accompany the software.

The program produces five files of output called RBno6.dat, RBno7.dat, RBno9.dat, RBno10.dat and cluster.dat.

The file RBno6.dat records a copy of the raw data, the parameter values selected, plus some basic summary statistics of the data. It is also used to log the path and name of the file read and to record any error messages from the program.

The file RBno7.dat contains the minimal output required for an analysis.

The file RBno9.dat contains output for further graphical analyses that must then be cut and pasted as input to some other package.

The file RBno10.dat contains the briefest summary of the most important of those indices and probabilities output to RBno7.dat.

The file Cluster.dat, contains x and y coordinates in columns 1 and 2, ordered by the corresponding cluster indices in column 3; you can ignore column 4; this file can be read straight into the SURFER® mapping program or similar software for mapping and interpolation of the clustering indices.

Note: (i) both files RBni5.dat and RBni8.dat must be present in the same folder as RBrel14.exe in order to run the program, unless you are using the interactive mode, AND none of the files RBno6.dat, RBno7.dat, RBno9.dat, RBno10.dat or Cluster.dat must be present when the program is run.
In the case that you are using the interactive mode, the files must not be in the project folder. If the files are present, you must rename or delete any existing versions of RBno7.dat, RBno9.dat, RBno10.dat and Cluster.dat. The presence of these files could cause RBrel14.exe to fail without any obvious error message, although an error message should be logged to Rbno6.dat.
Also note that, if present, RBno6.dat will be overwritten without warning.


Input File Structure

With the above provisos, running the program is easy! This is what you must do. First, put the n records in your data into file rbni5.dat, in the following form:

x co-ordinate 1 y co-ordinate 1 count 1
x co-ordinate 2 y co-ordinate 2 count 2
x co-ordinate 3 y co-ordinate 3 count 3
.                              .                      .
.                              .                      .
.                              .                      .
x co-ordinate n y co-ordinate n count n

Where the x & y co-ordinates for each pair, (xk, yk), k=1,...,n, should be read in as real numbers and the count, ck, k=1,...,n, should be read in, on the same line, as an integer, with no decimal point. No more than 2000 records can be analyzed in the current version.

Secondly, specify two parameters in the file RBni8.dat as follows: On line one, specify an integer seed, iseed, between 1 and 30,000, for the random number generator. Specifying the same seed in successive runs of the program will generate identical randomizations; specifying a different value will result in different randomizations. On line two, specify an integer value, k5psim, between 1 and 153 that will determine the number of randomizations done. The value of k5psim relates to how many blocks of 39 randomizations are performed. Within the program, the value of k5psim is multiplied by 39 to give the total number of randomizations performed; the result is denoted nsims in the output. If you can afford the time required to perform the randomizations, then you should use the largest value of k5psim possible, 153. Hence you should put the two values required into file RBni8.dat in the following form:
iseed
k5psim

Note: If you are using the interactive mode, you are prompted to enter iseed and k5psim, and RBni8.dat is not used.

More detailed instructions appear in Red_Blue_Handbook.pdf