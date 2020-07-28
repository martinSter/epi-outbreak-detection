# Selecting nodes for outbreak detection with greedy optimization

The code in this repository is heavily based on beautiful (and fast!) C code written by Petter Holme. Thanks for sharing it. The code can be found here: https://github.com/pholme/tsir. By studying his code in detail, I became a better C programmer (still ways to go...).

The code can be compiled as follows:
1. Create a directory o for the object files.
2. Type make to compile with gcc.
3. Then run the code with `./main data/escort.csv 0.3 100 16704718213277798935`.

We first provide the network file, then the transmission probability, the average recovery rate, and finally the random seed for the random number generator. Note that you can generate the seed by running the Python script *rnd.py*.

