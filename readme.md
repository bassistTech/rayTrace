# Python Ray Tracing Project

Francis Deck, 11-2-2023

## Subject matter

*If you were searching for the kind of "ray tracing" used in computer graphics, you came to the wrong place*.

This is ray tracing for modeling and analyzing optical designs such as lenses, mirrors, and gratings.

## Purpose

This project started as a learning exercise, to see if I could gather and document the formulas needed for optical ray tracing, and demonstrate some basic applications in a Jupyter notebook. And then, test my understanding by coding it and see if it works. In the future, I might use it to try out some unconventional design ideas for fun.

At some point it became obvious that keeping all of the code in notebooks was awkward. I moved the "guts" of the ray tracing code to a more traditional Python package.

I'm not trying to supplant commercial design software such as Zemax, Oslo, etc. By day, I use one of those programs, and am satisfied with it.

## What it does so far
* "Standard" conic section surface, which includes the spherical surface as a special case, and a cylindrical version
* "Paraxial" surface
* Refraction and reflection
* Plane diffraction grating
* Coordinate breaks
* Crude ray aiming
* Optimization

## Immediate plans
* Demonstrate more stuff
* Make breaking changes!
* As "features" emerge, move them out of the notebooks and into the Python library.

## Non plans
* Any kind of GUI or "software"

## My evaluation so far
Since I use one of the "big" commercial packages at my day job, it's easy to identify pro's and con's versus my own "package." The commercial packages still have some outstanding advantages that would be hard to live without:

* Computation speed. Analyses of any complexity require tracing *lots* of rays. This is important for optimization and image analysis.

* Proliferation of features. The sheer number of these features embodies years of effort. Some of them incorporate domain knowledge that I don't possess, either because it's proprietary or I haven't learned it yet. If you want to add them to my package, you have to program them yourself.

* Worth asking: Are these features necessary, or are the commercial packages just bloated? My experience is that when you run into an interesting obstacle or problem, it's extremely useful to look through the lists of operands for one that does exactly what you need. This is comparable to programmers having access to huge libraries.

## Overview of contents

**Entire repo** can be installed as a Python package on your machine using:

	git clone https://github.com/bassistTech/rayTrace.git
	cd rayTrace
	pip install -e .

The **-e** flag causes pip to install the package without moving it into your Python directory, so you can keep editing all of the files. I'd be surprised if you could use my package right now for anything interesting without forking it, so it doesn't really belong in your *site-packages*.

**Docs** are notebooks that document the theory and evolution of the program, but are not meant for long term use, and are not kept up to date with changes to the library. You can think of these notebooks as my lab notebooks, through which I developed the primitive code. **rayTraceSlow.ipynb** was the original notebook, but is written in mostly bare Python, so it's slow, but perhaps easier to read. The transition to array math using numpy can be a bit abrupt for a beginning programmer. **rayTrace.ipynb** is more or less the same program, but shows the transition to numpy.

**Notebooks** are projects and tests using the Python package.

## Bibliography of existing packages

No Python package is worth starting unless there's a perfectly good package that already exists. ;-) That's why there are half a million packages. Anyway, here are some other packages worth looking into:

https://github.com/DCC-Lab/RayTracing looks powerful, but is not a 3-d ray tracer, and does not handle aspherics. I'm going to dig into this code and see if there's anything I can borrow.

https://github.com/mjhoptics/ray-optics is two-dimensional tracing, but seems to have a lot of stuff.

https://github.com/quartiq/rayopt is three-dimensional tracing. I started digging through it, but it's deeply object oriented and was hard for me to understand. Also, it contains some warnings and TODOs, and doesn't appear to be maintained. I'd rather work with code that I can see to the bottom of, for this project.

https://github.com/mess42/pyrate seems worth looking into. Hard to figure out from the repo what it is and does.

https://github.com/rfrazier716/PyRayT looks promising, but is preliminary. I'm going to keep my eye on this.

https://github.com/chbergmann/OpticsWorkbench is a ray tracing add-on for FreeCAD. I haven't explored it at all, but I like the concept.

https://tmurphy.physics.ucsd.edu/astr597/exercises/raytrace-3d.pdf has a closed form solution for the intersection of a ray and a conic section surface, with Python example code.

## Other references

https://raytracing.github.io/books/RayTracingInOneWeekend.html is a book on ray tracing, for computer graphics. I'll be looking through it for insights on how I could speed up my program, since graphics tend to be quite performance intensive.