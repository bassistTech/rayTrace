# Python Ray Tracing Project

Francis Deck, 8-17-2024

## Subject matter

*If you were searching for the kind of "ray tracing" used in computer graphics, you came to the wrong place*.

This is ray tracing for modeling and analyzing optical designs such as lenses, mirrors, and gratings. Hopefully it's the same physics. The math should be similar, but computer graphics optimizes for speed, and optics design optimizes for accuracy. I've looked in both optics and graphics literature for useful formulas.

If you're a beginner at Python, using this package will be a challenge, but it also might be a good way to learn Python. If you're new at both optics and Python, prepare for a steeper learning curve, but have fun.

## Purpose

This project started as a learning exercise, to see if I could gather and document the formulas needed for optical ray tracing, and demonstrate some basic applications in a Jupyter notebook. In the future, I might use it to try out some unconventional design ideas for fun.

At some point it became obvious that keeping all of the code in notebooks was awkward. I moved the "guts" of the ray tracing code to a more traditional Python package.

I'm not trying to supplant commercial design software such as Zemax, Oslo, etc. By day, I use one of those programs, and am satisfied with it.

## What it does so far

I think I've gotten this to the point where I could use it for real designs, *provided that I extend the package with more Python code when I run into a roadblock*. 

* "Standard" conic section surface, which includes the spherical surface, and a cylindrical version
* "Paraxial" surface
* Refractive materials (glasses) and mirrors
* A crude glass table
* Plane diffraction grating
* Coordinate breaks
* Crude ray aiming
* Optimization is crude but works

## What it doesn't do

* I draw the line at creating any kind of GUI or general purpose "software." My reason is that doing so increases the size and complexity of the project, exponentially.
* Turning this into an optics textbook. That's too much, and there are already good textbooks.

## My evaluation so far

Since I use one of the "big" commercial packages at my day job, it's easy to identify pro's and con's versus my own code. The commercial packages still have some outstanding advantages that would be hard to live without, either for commercial design or education:

* Computation speed, when it's necessary to analyze lots of rays.

* My program doesn't handle configurations, tolerance analysis, thermal design, CAD output... The list goes on and on.

* Proliferation of features. The sheer number of these features embodies years of effort. Some of them incorporate domain knowledge that I don't possess, either because it's proprietary or I haven't learned it yet. If you want to add them to my package, you have to program them yourself.

* Worth asking: Are these features necessary, or are the commercial packages just bloated? My experience is that when you run into an interesting obstacle or problem, it's useful to look through the lists of operands for one that does exactly what you need. This is comparable to programmers having access to huge libraries.

* Using my package still requires proficiency in Python. In my opinion, nobody has come up with a "no code" or even "low code" free alternative to the big commercial packages.

## Overview of contents

**Entire repo** can be installed as a Python package on your machine using:

	git clone https://github.com/bassistTech/rayTrace.git
	cd rayTrace
	pip install -e .

The **-e** flag causes **pip** to install the package without moving it into your Python directory, so you can keep editing all of the files. I'd be surprised if you could use my package right now for anything interesting without modifying it, so it doesn't really belong in your *site-packages*.

**Docs** are notebooks that document the theory and evolution of the program, but are not meant for long term use, and are not kept up to date with changes to the library.

**Notebooks** are projects and tests using the Python package.

**Tutorials** are projects intended to walk you through using the package, step by step.

**Tests** are my own notebooks for verifying that I haven't broken anything.

## Coding style

I've used object-oriented programming (OOP) in a lot of projects, but this project doesn't use it. Objects can do a couple of things:

1. Encapsulate complicated data structures, while allowing the compiler (or code profiler) to detect incorrect uses. But we don't have a large number of structures: For now, just surfaces and rays. These are easily expressed as **dict**s and a big Numpy array.

2. Manage projects requiring many layers of abstraction, particularly GUIs. But this package is pretty flat. Most of the operations are only a couple of layers deep. I want to keep it that way.

**Surfaces** are specified in a traditional fashion, but using a **dict** structure instead of a tabular editor. A **dict** allows the list of surfaces to be *inhomogeneous*. What I mean is, each surface in the list can have properties specific to its surface type, such as conic, cylindrical, grating, etc.

**Geometry** is the same list of surfaces, but each surface gets more properties such as its location and orientation within a global coordinate system.

**Rays** are represented by a big multi-dimensional **Numpy** array. Two "columns" of the array are for the intersection point and direction vector (in global coordinates) of each ray. Additional "columns" track properties of the rays, including wavelength. This is kind of an abuse of of the array structure, but makes for quick computation. On the other hand, if it's not really speeding things up (which I intend to find out), then I'll create a cleaner structure.

The ray list is passed into functions, and modified in place. This is not functional programming, but seems appropriate for the typical analysis workflow. In any event, an equivalent C program would do the same thing, since a C function can't return an object.

An interesting bit of trivia is that Python is highly optimized for fast access of **dict** structures. I've got them in the inner loops of my ray tracing code, and my tests have shown that there's no speed penalty.

I think I've made effective use of **Numpy**, but it could always improve. I've tested a number of further improvements, most of which produce little or no gain at the expense of readability. Tiny gains are not worthwhile for me, because I don't want my code to be opaque to critique. I spent some time converting my code to work with **Numba**. Once again, little or no improvement, and the code became pretty barbarous.

## Bibliography of existing packages

No Python package is worth starting unless there's a perfectly good package that already exists. That's why there are half a million packages. Anyway, here are some other packages worth looking into:

https://github.com/Garchupiter/Kraken-Optical-Simulator is another Python ray tracing library. It strives for accuracy, and handles a reasonable range of design options. But it suffers from the same drawback as my package: Getting the ray tracing engine working is easy. Doing anything useful with the rays requires a lot of Python programming.

https://github.com/DCC-Lab/RayTracing looks powerful, but is not a 3-d ray tracer, and does not handle aspherics. I'm going to dig into this code and see if there's anything I can borrow.

https://github.com/mjhoptics/ray-optics is also two-dimensional tracing, but seems to have a lot of stuff.

https://github.com/quartiq/rayopt is three-dimensional tracing. I started digging through it, but it's deeply object oriented and was hard for me to understand. Also, it contains some warnings and TODOs, and doesn't appear to be maintained. I'd rather work with code that I can see to the bottom of, for this project.

https://github.com/mess42/pyrate seems worth looking into. Hard to figure out from the repo what it is and does.

https://github.com/rfrazier716/PyRayT looks promising, but is preliminary. I'm going to keep my eye on this.

https://github.com/chbergmann/OpticsWorkbench is a ray tracing add-on for FreeCAD. I haven't explored it at all, but I like the concept.

## Other references

https://tmurphy.physics.ucsd.edu/astr597/exercises/raytrace-3d.pdf has a closed form solution for the intersection of a ray and a conic section surface, with Python example code.

https://raytracing.github.io/books/RayTracingInOneWeekend.html is a book on ray tracing, for computer graphics. I'll be looking through it for insights on how I could speed up my program, since graphics tend to be quite performance intensive.

https://refractiveindex.info is where I got my glass catalog.
