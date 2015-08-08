Again, this is an **old** project previously under subversion, and pre-C++11.

## Connected Components

In the fall of 2010 I took a Computer Algebra course during which we slowly built up a Maple program to calculate the different partitions in the plane that the graphs of arbitrary bivariate polynomials made. The polynomials where over the rational numbers using infinite precision integers for numerators and denominators.

I wanted to write a program that solved the same problem, but written in an object-oriented way in C++, making good use of encapsulation, classes, namespaces, the STL, and template classes/functions. I used the Boost libraries and some of the new features in the Technical Report 1 (TR1) extensions to the C++ standard. I used the GCC compiler and GDB debugger.

The development work was done with the Code::Blocks IDE on a Windows desktop and an OSX laptop. I used subversion for revision control and backup, hosted on [assebla](https://www.assembla.com).

I wanted to document my code well and keep the documentation inline with the source, so used QT-style comments and doxygen to generate the documentation from the source. The comments are QT-style (as opposed to JavaDoc-style comments) because I was considering using QT for a graphical layer and thought there might be value in being able to use doxygen or the QT documentation generator.

I wanted to use infinite precision arithmetic so I needed extra libraries. After research I decided on using the GiNaC library for general unlimited-precision computation. GiNaC itself uses a lower level number library called the Class Library for Numbers (CLN). CLN itself can use a very low level library that is written mostly in assembly called the GNU Multi Precision (GMP) library. GMP is used by many computational software packages like Mathematica and Maple.

There is a GUI using the native GUI toolkit wxWidgets (formerly wxWindows).

## Requirements

The GiNaC, CLN, and GMP libraries.
