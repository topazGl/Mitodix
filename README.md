# Mitodix - Mitosis detection
Supplementary output videos are found here: https://goo.gl/aJZKv7


## Cite:
T. Gilad, J. Reyes, J.Y. Chen, G. Lahav, and T. Riklin Raviv. "Fully Unsupervised Symmetry-Based Mitosis Detection in Time-Lapse Cell Microscopy". Bioinformatics. 2018 Dec 24.

Topaz Gilad, Mark-Anthony Bray, Anne E. Carpenter, and Tammy Riklin Raviv. "Symmetry-based mitosis detection in time-lapse microscopy." In Biomedical Imaging (ISBI), 2015 IEEE 12th International Symposium on, pp. 164-167. IEEE, 2015.

## Main functions:
- Our triangulation adaptation together with an example is found in: Code/Algorithms/StochasticDelaunay.m
- Our similarity measure is here: Code/Algorithms/wcorr2.m
- Integer prog. problem class is in: Code/Cells/IntProg.m
- Our symmetry-preserving registration is a function of the "Frame" object: Code/Cells/@Frame/GetCellsPatch.m

## Code Guidelines:
- Code files are inside the code folder. Place your input dataset in the "Data" folder.
- Inputs are both gray-scale sequence  and its binary / labeled segmentation. Note that the algorithm assumes touching cells are separated in the input segmentation (different connected components clusters).
- Make sure to edit Mitodix.config XML file in order to add you data set. Frame rate is only used to calculate frame rate for the output *.avi.
- main_detection.m has an example how to create and use a 'Mitodix' object to detect mitosis events in a time-lapse sequence.
