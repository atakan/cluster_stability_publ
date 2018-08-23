About this document
===================

I (MAG) decided to document the program and the steps taken to do so. I use
Sphinx possibly with various expansions.

I used sphinx-quicstart in docs/ directory. I needed to choose seperate
build and source directories since github did not like _build name,
probably because of suffix. So, I also changed the prefix for _static
etc. into ".", though this only changes the subdirectories of source
directory and keeps some other _static etc. names (this is problematic,
see the next paragraph). Finally, I chose
mathjax over png/svg math and Makefile over Windows command file and
opted to create a .nojekyll file.  I believe, I chose the defaults
for the rest of the options.

It turns out you *must* have a .nojekyll file in the *root directory*.
Otherwise jekyll interferes and ignores directories starting with an
underscore.

LaTeX to reStructuredText
-------------------------
I am much more used to writing in LaTeX for math, so I actually write in
LaTeX and convert the result to reStructuredText in vim, using the
following substitution commands::

    :%s/\$\(.\{-}\)\$/:math:`\1`/gc converts $...$ 
    :%s/\\\[\n\(\_.\{-}\)\\\]/^M.. math::^M   \1^M/gc  converts \[...\]

Note that ^M is obtained by Ctrl-V and Enter.
