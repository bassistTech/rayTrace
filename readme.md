# Python Ray Tracing Project

Francis Deck, 11-2-2023

## Subject matter

*If you were searching for the kind of "ray tracing" used in computer graphics, you came to the wrong place*.

This is ray tracing for modeling and analyzing optical designs such as lenses, mirrors, and gratings.

## Purpose

This project started as a learning exercise, to see if I could gather and document the formulas needed for optical ray tracing, and demonstrate some basic applications in a Jupyter notebook. My goal is to find and document the formulas for optical ray tracing. And then, test my understanding by coding it and see if it works. In the future, I might use it to try out some unconventional design ideas for fun.

I'm not trying to supplant commercial design software such as Zemax, Oslo, etc. By day, I use one of those programs, and am perfectly happy with it.

## What it does so far
* "Standard" conic section surface, which includes the spherical surface as a special case, and a cylindrical version
* Refraction and reflection
* Plane diffraction grating
* Coordinate breaks

## Immediate plans
* Handling the basic properties of glass
* Demonstrate a reasonably complex example
* Demonstrate an analysis
* Demonstrate optimization

## My evaluation so far
Since I use one of the "big" commercial packages at my day job, it's easy to identify pro's and con's versus the my own "package." The commercial packages still have some outstanding advantages that would be hard to live without:

* Computation speed. Analyses of any complexity require tracing *lots* of rays. This is important for optimization, which is a super power when it comes to optics design.

* Proliferation of features. The sheer number of these features embodies years of effort. Some of them incorporate domain knowledge that I don't possess, either because it's proprietary or I haven't learned it yet. If you want to add them to my package, you have to program them yourself. Programming could easily end up taking more time than designing.

* GUI, as in, None.

As I mentioned, optimization is a super-power. So my package can make simple diagrams that illustrate principles, but is actually quite weak for serious design.

## Bibliography of existing packages

No Python package is worth starting unless there's a perfectly good package that already exists. ;-) That's why there are half a million packages. Anyway, here are some other packages worth looking into:

https://github.com/DCC-Lab/RayTracing looks powerful, but is not a 3-d ray tracer, and does not handle aspherics. I'm going to dig into this code and see if there's anything I can borrow.

https://github.com/mjhoptics/ray-optics is two-dimensional tracing, but seems to have a lot of stuff.

https://github.com/quartiq/rayopt is three-dimensional tracing. I started digging through it, but it's deeply object oriented and was hard for me to understand. Also, it contains some warnings and TODOs, and doesn't appear to be maintained. I'd rather work with code that I can see to the bottom of, for this project.

https://github.com/mess42/pyrate seems worth looking into. Hard to figure out from the repo what it is and does.

https://github.com/rfrazier716/PyRayT looks promising, but is preliminary. I'm going to keep my eye on this.