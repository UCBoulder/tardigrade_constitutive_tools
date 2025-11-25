###############################
tardigrade\_constitutive\_tools
###############################

*******************
Project Description
*******************

A collection of tools useful for constitutive modeling. These tools are
intended to reduce the burden in creating a new constitutive model from
scratch enabling a faster turn-around for model development. These tools
should be as general as possible to avoid cluttering the database with
extraneous things.

Information
===========

Developers
==========

* Nathan Miller nathanm@lanl.gov
* Kyle Brindley kbrindley@lanl.gov

************
Dependencies
************

The developer dependencies are found in ``environment.txt``.

.. code-block:: bash

   $ conda create --name tardigrade_constitutive_tools-dev --file environment.txt

**************************
Building the documentation
**************************

.. warning::

   **API Health Note**: The Sphinx API docs are a work-in-progress. The doxygen
   API is much more useful

.. code-block:: bash

   $ pwd
   /path/to/tardigrade_constitutive_tools
   $ cmake -S . -B build
   $ cmake --build build --target Doxygen Sphinx

*****************
Build the library
*****************

1) Build just the library

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade_constitutive_tools
      $ cmake -S . -B build
      $ cmake --build build --target tardigrade_constitutive_tools

****************
Test the library
****************

.. code-block:: back

   $ pwd
   /path/to/tardigrade_constitutive_tools
   $ cmake -S . -B build
   $ cmake --build build --target tardigrade_constitutive_tools test_tardigrade_constitutive_tools
   $ ctest --test-dir build

*******************
Install the library
*******************

Build the entire project before performing the installation.

4) Build the entire project

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade_constitutive_tools
      $ cmake -S . -B build
      $ cmake --build build --target all

5) Install the library

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade_constitutive_tools
      $ cmake --install build --prefix path/to/root/install

      # Example local user (non-admin) Linux install
      $ cmake --install build --prefix /home/$USER/.local

      # Example install to an active conda environment
      $ cmake --install build --prefix $CONDA_PREFIX

***********************
Build the Conda package
***********************

.. code-block:: bash

   $ conda mambabuild recipe --no-anaconda-upload -c conda-forge --output-folder conda-bld

*****************
Local development
*****************

In some cases it is not convenient to pull down every repository required but it may be desired that local
versions of the repository are used. An example of when this may be needed is if development is across
multiple libraries and is proceeding faster than collaborators can check in results. In this case, and
outside of developers no-one should need to do this, a version of the code using local repositories can be
built.

To perform in-source builds of upstream libraries, the active Conda environment can NOT include installed versions of
the upstream libraries to be built in-source with the current project. It is possible to mix sources with some upstream
libraries coming from the active Conda environment and others built in-source from a Git repository. Developers may
build minimal working Conda environments from the Python Modules discussion.

1) Build and activate a minimal Conda development environment

   .. code-block:: bash

       $ conda env create --file configuration_files/environment.yaml
       $ conda activate environment

2) Define convenience environment variables

   .. code-block:: bash

      $ tardigrade_error_tools=/path/to/my/tardigrade_error_tools
      $ tardigrade_error_tools_version=origin/dev
      $ tardigrade_vector_tools=/path/to/my/tardigrade_vector_tools
      $ tardigrade_vector_tools_version=origin/dev

3) Perform the initial configuration. Note that the environment variables are mutually independent. Each variable can be
   used alone or in arbitrary combinations. The default values are found in the root ``CMakeLists.txt`` file. The ``PATH``
   variables can accept anything that the [``CMake``
   ``FetchContent``](https://cmake.org/cmake/help/latest/module/FetchContent.html) ``GIT_REPOSITORY`` option can accept.
   The ``GITTAG`` variables will accept anything that the [``CMake``
   ``FetchContent``](https://cmake.org/cmake/help/latest/module/FetchContent.html) ``GIT_TAG`` option can accept.

   .. code-block:: bash

      # View the defaults
      $ grep _TOOLS_ CMakeLists.txt
      set(TARDIGRADE_ERROR_TOOLS_PATH "" CACHE PATH "The path to the local version of tardigrade_error_tools")
      set(TARDIGRADE_ERROR_TOOLS_GITTAG "" CACHE PATH "The path to the local version of tardigrade_error_tools")
      set(TARDIGRADE_VECTOR_TOOLS_PATH "" CACHE PATH "The path to the local version of tardigrade_vector_tools")
      set(TARDIGRADE_VECTOR_TOOLS_GITTAG "" CACHE PATH "The path to the local version of tardigrade_vector_tools")

      $ Build against local directory paths and possible custom branch
      $ pwd
      /path/to/tardigrade_constitutive_tools
      $ mkdir build
      $ cd build
      $ cmake .. -DTARDIGRADE_ERROR_TOOLS_PATH=${my_tardigrade_error_tools} -DTARDIGRADE_ERROR_TOOLS_GITTAG=${tardigrade_error_tools_version} -DTARDIGRADE_VECTOR_TOOLS_PATH=${my_tardigrade_vector_tools} -DTARDIGRADE_VECTOR_TOOLS_GITTAG=${tardigrade_vector_tools_version}

4) Building the library

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade_constitutive_tools/build
      $ make

***********************
Contribution Guidelines
***********************

Git Commit Message
==================

Begin Git commit messages with one of the following headings:

* BUG: bug fix
* DOC: documentation
* FEAT: feature
* MAINT: maintenance
* TST: tests
* REL: release
* WIP: work-in-progress

For example:

.. code-block:: bash

   git commit -m "DOC: adds documentation for feature"

Git Branch Names
================

When creating branches use one of the following naming conventions. When in
doubt use ``feature/<description>``.

* ``bugfix/\<description>``
* ``feature/\<description>``
* ``release/\<description>``

reStructured Text
=================

[Sphinx](https://www.sphinx-doc.org/en/master/) reads in docstrings and other special portions of the code as
reStructured text. Developers should follow styles in this [Sphinx style
guide](https://documentation-style-guide-sphinx.readthedocs.io/en/latest/style-guide.html#).

Style Guide
===========

This project does not yet have a full style guide. Generally, wherever a style can't be
inferred from surrounding code this project falls back to
[PEP-8](https://www.python.org/dev/peps/pep-0008/)-like styles. There are two
notable exceptions to the notional PEP-8 fall back:

1. [Doxygen](https://www.doxygen.nl/manual/docblocks.html) style docstrings are
   required for automated, API from source documentation.
2. This project prefers expansive whitespace surrounding parentheses, braces, and
   brackets.
   * No leading space between a function and the argument list.
   * One space following an open paranthesis ``(``, brace ``{``, or bracket
     ``[``
   * One space leading a close paranthesis ``)``, brace ``}``, or bracket ``]``

An example of the whitespace style:

.. code-block:: bash

   my_function( arg1, { arg2, arg3 }, arg4 );

The following ``sed`` commands may be useful for updating white space, but must
be used with care. The developer is recommended to use a unique git commit
between each command with a corresponding review of the changes and a unit test
run.

* Trailing space for open paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([({[]\)\([^ ]\)/\1 \2/g' <list of files to update>

* Leading space for close paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([^ ]\)\([)}\]]\)/\1 \2/g' <list of files to update>

* White space between adjacent paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([)}\]]\)\([)}\]]\)/\1 \2/g' <list of files to update>
