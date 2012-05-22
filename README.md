TODO
====

Cloaking designs
----------------

*   Anti-reflection
*   Wrap-cloak
*   Front/back/side cloak
*   Channeler

Mimic designs
-------------

*   Uncloak (mimic the prev. cloaked object)
*   Mimic n = -1 material
*   Lens mimic
*   Lithography mimic


Purpose
=======

This package implements an "objective-first" approach to
the design of nanophotonic waveguide couplers in two dimensions.

Feel free to use and modify the code, for this reason
the implementation here is simple and heavily documented.

See the following paper for details on objective-first optimization:
    J. Lu, J. Vuckovic, 
    "Objective-first design of nanophotonic waveguide couplers," 
    (to be submitted)

Draft can be found at:
https://github.com/JesseLu/misc/blob/master/presentations/ob1_wg_paper/paper.pdf


Installation
============

This package requires the basic version of Matlab,
The CVX software (www.stanford.edu/~boyd/cvx)
is included verbatim in this repository as well.

To install, just unzip all files in a directory.
Then to use, just open matlab from the directory.

You can also run example.m to for a demo of the package.
To do this, simply type 'example' in the Matlab command line.


Problem Specification
=====================

This package attempts to solve the following nanophotonic design problem; given

1.  an arbitrary input waveguide mode on the left,
1.  an arbitrary output waveguide mode on the right, and
1.  a central "design box" between the two; 

find a dielectric structure within the "design box" which will convert from 
the input to the output waveguide mode as efficiently as possible.


Limitations
-----------

1.  There is no guarantee on the performance of the design, meaning that
    the user is not even guaranteed that the final design will improve upon
    the initial one.
    This is a natural outcome of employing an objective-first approach.
    In practice, this software routinely designs very high efficiency 
    (~95%) couplers.
1.  Although a discretized permittivity (or dielectric structure)
    is usually desired, the current version of the package only limits the
    values of the permittivity to a continuous range.
    The ability to produce completely binary structures is being actively
    developed for the next iteration of this package.


Usage
=====

Basic usage consists of the following three commands:

1.  setup(), determine the input and output waveguides modes;
1.  solve(), run the design algorithm;
1.  simulate(), determine the actual performance of the design.

See the example.m file for an example optimization.

Please consult the documentation of use help command, 
where command is setup, optimize, or simulate, 
for more information on how to use these commands.


License
=======

This code is public domain and comes with no guarantees. 
Feel free to use it however you like.

