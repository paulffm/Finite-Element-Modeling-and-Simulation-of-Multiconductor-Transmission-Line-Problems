The documentation is generated in this folder. The program [Sphinx](https://www.sphinx-doc.org/en/master/index.html) needs the present folder structure for the generation. 

If you habe Sphinx not installed, you can install it with anaconda with the command

``conda install -c anaconda sphinx``

After generation, the documentation files are saves in the `_build` folder. 

### Generation of the documentation

Open the Anaconda Prompt and activate the correct environment. Then exeute

``make html``

for a documentation in html format or 

`` make latexpdf``

for a documentation in LaTeX format and pdf. The results are stored in `_build/html` and `_build/latex`, respectively. 

### Settings

For both output formats, and for all the other possible output formats (see [here](https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-latex-output)), settings regarding the appearance are stored in the file `conf.py`. The main structure of the documentation is made in `index.rst`. This is a [reStructuredText](https://docutils.sourceforge.io/rst.html#user-documentation) file and anything can be added to it, als long it is in reStructuredText format. 

### Open the doc (html)

To open the html version of the documentation, start `_build/html/index.html` with your browser.
Alternatively, you can use the file `open_doc.bat` on Windows and the file `open_doc.sh` on Mac to open the documentation.

