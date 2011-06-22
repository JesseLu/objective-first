Objective-first optimization package
====================================

This package is a research and development build with the goal of formulating a design algorithm for wave structures.


What is objective-first (Ob-1) optimization?
--------------------------------------------

Objective-first optimization is a unique strategy for finding a physical structure and field which meet a design objective. 

Understanding the idea behind Ob-1 optimization is most easily understood by comparing it with traditional optimization techniques.

Traditional: Satisfy physics while improving the value of the design objective.

Ob-1: Satisfy the design objective while decreasing the error in the physical equations (physics residual).


Validation steps
----------------

The following steps form a progressive validation procedure that allows one to determine whether or not a specific optimization strategy is viable.

1.  Cast the physics of the wave equation in a matrix-vector form and independently solve for both the field and unbounded structure variables using a generic linear algebra solver.
    This step validates the correctness of the matrices and, in doing so, the formulation of the problem.

2.  Produce an equivalent result using a gradient descent method on both field and unbounded structure variables independently.
    This step validates the accuracy of the gradient.

3.  Joint gradient-descent on field and unbounded structure variables.
    This is likely the first non-convex problem attempted and serves as a key milestone in evaluating whether or not Ob-1 will be effective.

4.  Joint gradient-descent on field and bounded structure variables.
    Use the level-set method to calculate the shape derivatives for the boundaries of the structure. Care will be needed in testing the shape derivatives. 


