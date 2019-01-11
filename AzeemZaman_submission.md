# Harvard Data Science Animation Contest Sumbission




# LDA and data separation in the orthogonal space

**Author**: Azeem Zaman

**Affiliation**: PhD Student, Department of Statistics, Harvard University


**Artifact:**

![](https://raw.githubusercontent.com/azeemzaman/animation-contest/master/ArtifactExamples/AzeemZaman_artifact.gif)

**Code:** *Code to create this gif can be found [here](https://github.com/azeemzaman/animation-contest/blob/master/CodeExamples/AzeemZaman_code.R)*

### Explanation

Consider the problem of distinguishing two sets of points. If one wants to do this using a line, a very effective and widely used method to do this is support vector machines (SVM). SVM gives a method for finding a hyperplane that minimizes the number of points on the wrong sidie of the boundary. A special case of this problem was consider much earlier by Fisher. Suppose we have two normal distributions with the same covariance (but different means). We can propose a separating hyperplane (in two dimensions this is just a line) to distinguish the two distributions. Fisher's solution is known as Linear Discriminant Analysis (LDA), which uses the line that resulting from classifying a point based on a comparison of the likelihoods.

As mentioned above, LDA can easily been seen as an optimal separating hyperplane in the framework of SVM or as a likelihood based classifier. A nice geometric property of LDA is related to the separation of the projected distributions. Suppose we propose a separating hyperplane (the solid black line) and project the data onto the space orthogonal to this hyperplane (shown by the red line). The data points will still have a normal distribution. We can look at the separation of these normal distributions in the projected space and try to maximize it. It turns out this separation is maximized if we take the hyperplane to be the one given by LDA! 

The left half of the animation shows various potential separating hyperplanes and their orthogonal spaces. On the orthognal space we plot the densities of the projected data. These densities are shown in the right half of the figure for clarity. The separation is maximized when the black line is close to the LDA solution. When the black line is close to the direction between the two means, then there is very little separation.
