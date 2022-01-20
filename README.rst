Calculating the Molecular Weight of a Chemical
--------------------------------------------------

Background
+++++++++++

The molecular weight (MW) can be algebraically calculated from any chemical formula that adheres to `chemical conventions <https://en.wikipedia.org/wiki/Chemical_formula>`_. 

The `ChemMW` object of ``ChemW`` parses a chemical formula string -- which consists of any combination of elements, however sensible or outrageous -- and precisely calculates the MW of the chemical formula, based upon the current physical accuracy that is embedded in the ``periodic_table`` dictionary of the ``chemicals`` module.

The `PHREEQdb` object of ``ChemW`` parses a `PHREEQ database <https://www.usgs.gov/software/phreeqc-version-3>`_ via the `ChemMW` object. The object exports a JSON file that consolidate the elements and minerals, and importantly mineral masses, of the database. This unique application of the `ChemMW` object has been applied as the pivotal means of predicting the mass of mineral scaling in the `ROSSpy module <https://pypi.org/project/ROSSpy/>`_ for reverse osmosis research.

The ``ChemW`` module is offered with the `MIT License <https://opensource.org/licenses/MIT>`_\.

Usage
++++++

+++++++++++++
installation
+++++++++++++

The following command are executed in a command prompt/terminal environment::
 
 pip install hillfit

+++++++++++
__init__
+++++++++++

The data environment, in a `Python IDE <https://www.simplilearn.com/tutorials/python-tutorial/python-ide>`_, is defined: 

.. code-block:: python

 import hillfit
 hf = hillfit.HillFit(x_data, y_data)

- *x_data* & *y_data* ``list`` or ``ndarray``: specifies the x-values & y-values, respectively, of the raw data that will be fitted with the Hill equation.

++++++++++++++++
fitting()
++++++++++++++++

The parameterized data is fitted to the Hill equation, with the following arguments and their default values:

.. code-block:: python

 hf.fitting(x_label = 'x', y_label = 'y', title = 'Fitted Hill equation', sigfigs = 6, view_figure = True)

- *x_label* & *y_label* ``str``: specifies the x-axis & y-axis labels, respectively, that will be applied to the regression plot for the raw data points and the fitted Hill equation.
- *title* ``str``: specifies the title of the regression plot for the raw data points and the fitted Hill equation.
- *sigfigs* ``int``: specifies the number of `significant figures <https://en.wikipedia.org/wiki/Significant_figures>`_ that will be used in printed instances of the fitted Hill equation.
- *view_figure* ``bool``: specifies whether the regression plot will be printed in the Python environment.

-----------------------------
Accessible content
-----------------------------
Many data sets and exported components of the fitted information are accessible through the ``hillfit`` model object. 

- *top*, *bottom*, *ec50*, & *nH* ``float``: The fitted parameters of the Hill equation are accessible via ``hf.top``, ``hf.bottom``, ``hf.ec50``, & ``hf.nH``, respectively.
- *fitted_xs* & *fitted_ys* ``list``: The x- and y-values of the fitted Hill equation are accessible via ``hf.x_fit`` & ``hf.y_fit``, respectively.
- *fitted_equation* ``str``: The fitted Hill equation, in an amenable string format for the `eval() built-in Python function <https://pythongeeks.org/python-eval-function/>`_ that allows the user to directly execute the string as a function for an "x" variable, is accessible via ``hf.equation``.
- *figure* ``matplotlib.figure``: The regression figure is available via ``hf.figure``.
- *x_data* & *y_data* ``ndarray``: The arrays of original data are available via ``hf.x_data`` & ``hf.y_data``, respectively.


++++++++++
export()
++++++++++

The fitted Hill equation, with its data points and parameters, and the regression information are exported to a designated folder through the following syntax and arguments:

.. code-block:: python

 hf.export(export_path = None, export_name = None)

- *export_path* ``str``: optionally specifies a path to where the content will be exported, where `None` selects the current working directory.
- *export_name* ``str``: optionally specifies a name for the folder of exported content, where `None` enables the code to design a unique folder name for the information.

Execution
+++++++++++

Hillfit is executed through the following sequence of the aforementioned functions, which is exemplified in the `example Notebook of our GitHub repository <https://github.com/freiburgermsu/hillfit/tree/master/examples>`_:

.. code-block:: python
 
 import hillfit
 hf = hillfit.HillFit(x_data, y_data)
 hf.fitting(x_label = 'test_x', y_label = 'test_y', title = 'Fitted Hill equation', sigfigs = 6, view_figure = True)
 hf.export(export_path = None, export_name = None)
