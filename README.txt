This program implements G-means clustering [1]

Set working directory to the folder with the data and source the R file to see plots of results from G-means. Plots are created every time there is a split, two views of x1 versus x2 and x2 versus x3.

To observe results on my simulated 2-d data, please comment out line 28 and source the R file. If you would like to use new test data, load it into a data frame called “data” after running first 28 lines, then continue running the code from line 29 onward.

Libraries required for this code are:
1. MASS
2. ellipse

One plot of x1 versus x2 after G-means completion is also attached.

[1] Hamerly, Greg, and Charles Elkan. "Learning the k in k-means." Advances in neural information processing systems. 2004.
