# MAE150A

This repository is designed for MAE 150A Intermediate Fluid Mechanics at UCLA.

To get started with jupyter, please do the following:

1. Download Julia to your computer from https://julialang.org

2. Start a Julia session.

3. Type `]` to enter the Julia package management system. At this prompt, type `add IJulia`. This downloads `jupyter`, which is a browser-based computing environment. The 'ju' part stands for Julia; the 'py' part stands for Python. You can return to the Julia prompt by pressing backspace.

If you already have jupyter set up on your computer, start here:

1. Enter the Julia package management system, if you are not already in it (see above). Type `add https://github.com/jdeldre/MAE150A`. This will add the MAE150A package to your own computer.

2. Now press backspace to return to the Julia prompt. At this prompt, type `using MAE150A`. This precompiles the package. Note: it will probably take several minutes.

3. After you have waited patiently for that to finish, type `MAE150A.open_notebooks()`. This will put you into the main index of notebooks. You can open each one and run it.

4. If you wish to use any of these notebooks as a starting point for your own modifications, please use the "Make a Copy..." option from the main jupyter menu.
