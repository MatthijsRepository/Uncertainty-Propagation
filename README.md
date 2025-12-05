## The calculation tool: what is it and what can you do with it?

The pythontool is a general uncertainty propagator that calculates uncertainty according to the Guide to the Uncertainty of Measurement, or GUM. 
The user specifies a set of variables in a text file, basic variables and derived variables.
- Basic variables are constants or measured variables
- Derived variables are specified by an equation and hence dependent on other variables. 
This way, the user defines an equation tree, or more strictly a directed acyclic graph. The variables are the nodes in this ‘tree’.

 ![Alt text](https://github.com/user-attachments/assets/89cbe44e-ae13-467f-9d0a-1dd895465d20)
 
In the textfile the user can specify the value of a variable, or express that the values are contained in a specific column of a CSV file. If no value is specified it means the value must be calculated from other variables.  
The user can also specify for each variable:
- the aggregation rule (summing or averaging)
- whether the variable is a rate of a quantity over time, or whether it is simply a quantity
- whether uncertainty is ‘maskable’. Uncertainty masking means that the uncertainty of a time series is not included where time-series is zero-valued. For instance to exclude uncertainty in G at night when calculating aggregated uncertainty. Note that setting this to true allows the uncertainty to be masked, not that it is masked. Use of masking must be specified in the jobscript.  

The user can specify uncertainty sources for each variable. Multiple uncertainty sources can act on the same variable. For each uncertainty source the user can specify:
- absolute or relative error
- distribution
- deviation (mean assumed 0)
- autocorrelation over time (currently 0 and 1 are supported, extension to exponential decay is planned)
- an optional multiplier to the magnitude, in the form of an equation. For instance, one can multiply the directional response error by a function of the solar zenith angle to more closely model its dependence on solar zenith angle.

## Points to keep in mind
This tool is meant to calculate the minimally achievable uncertainty in a quantity based on the specifications on the used measurement systems. In case the magnitude or characteristics of an uncertainty are unknown, it cannot be included in this calculation.  
Be aware that the results are always dependent on the used dataset. The uncertainty of a quantity $C = A * B$ with an uncertainty source in $B$, will always be dependent on the value(s) of $A$. To draw general conclusions, one needs to average or aggregate over a large amount of data.  
The tool can be used for ‘normal’ equations: regular arithmetic, mathematical operations, exponents, trigonometric functions and integrals. Be aware that the code takes as the sensitivity coefficient for a source simply the instantaneous value of the partial derivative of the measurement equation with respect to the variable the uncertainty source acts on – which is a first-order approximation. The validity of this approximation may be questionable in case of highly nonlinear functions and large relative uncertainties.


## Assumptions, calculation errors and other limitations
- Importantly, the tool works under the assumption that separate uncertainty sources are independent (not cross-correlated).
- Additionally, at present no smart-matching feature is implemented in case the same uncertainty source acts on a variable multiple times. For instance, when two variables in an equation are both dependent on the same temperature measurement, the uncertainty in the temperature measurement will appear twice in the total uncertainty equation. The code will treat these as separate sources with their own sensitivities, while in practice the sensitivities should be combined, leading to an error.  
In mathematical terms: $u_T^2 * (s_1^2 + s_2^2) \neq u_T^2 (s_1 + s_2)^2$.  
This means that in such cases there will be a slight under-estimation of the uncertainty.  
To fix this issue, one can scan the var.uncertainty.root_sources list for duplicates to detect such cases, and implement exception handling in case duplicate sources are detected. This should be done when calculating the total uncertainty and when retrieving the root source contribution split. 
- A final limitation of the code is the amount of data that can be processed at a time. Because of the nature of correlated uncertainties, the computational complexity of aggregating uncertainty over time scales quadratically with the included timeframe. It is therefore recommended to split calculations for large datasets into a set of smaller computations. 

It is recommended to expand the code to be able to include exponentially decaying autocorrelation, with an inclusion cutoff if the correlation is below a specified limit, resulting in a band correlation matrix. Correlated uncertainty aggregation over arbitrarily long timescales can then be performed iteratively and potentially parallelized for performance. Moreover, a short-circuit in case of zero correlation could be implemented to avoid redundant matrix multiplication with a diagonal matrix.  

Plotting functionality was designed with orientation to variables in mind. However, after implementation of looping through datasets it became apparent that plots are mostly based on data stored inside the `Results` object. No integrated plotting functionality is included for this object. TThus, plotting results from multiple job calls must be done manually.

## How the code works – general
The tool is developed in an object-oriented way. All variables are objects that store information on their own values, uncertainties, dependencies and properties such as calculation rules and state indicators.  
The operations on the variables, or on the equation tree as a whole, are performed by objects called engines. The calculation engine calculates variable values, the uncertainty engine calculates uncertainty, the time engine handles time matching between variables, et cetera.

#### User Interface
The JobHandler object is the main interface between the user and the internal functionality. The user can load equation trees and `CSVData` dataclass objects to this handler and specify which tasks should be executed. The job handler will then perform pre-execution checks, variable initialization, job execution and post-job result storing and data cleaning. The user needs to specify the preprocessing and main job routine to the jobhandler by defining a preprocessing and main function. The `JobHandler` object has wrapper functions for many of the main functionalities of the code that can be used inside these functions. Additionally, the JobHandler also has internal copies of all engines, so the user can also directly access the full engine functionality inside the preprocessing and main functions with the right syntax. The jobscript.py file contains an illustration on how this works.   
To load CSV data to variables, the code makes use of a `CSVData` dataclass. This is a small dataclass that contains the values, the start and end time and the timestep of the data. The `PandasCSVHandler` object can read a CSV to a Pandas dataframe, and can compile it to `CSVData` objects. It can also do this for a specific day. For convenience, data cleaning methods are owned by this dataclass.  
The usage of this dataclass is because of earlier architectural decisions. It is recommended to deprecate usage of this dataclass and move to a purely-pandas based method.

#### Dependencies
The code is purely python-based and makes use mostly of standard python libraries: `numpy`, `matplotlib`, `pandas`. The code also uses `pvlib` to perform solar zenith angle calculations. In case you don't want to use this functionality, simply do not use functionality related to solar zenith angles and comment out the lines in `solar_module.py` related to it.

## Overall workflow
- The user creates instances of the `JobHandler` and `PandasCSVHandler` classes.
- The user gives the equation tree text file to the Job Handler. The handler will compile an equation tree from this text file, check whether it is well-defined and non-circular. It will populate all variables with pointers to their dependencies, and it will prepare all variable equation executables.
- The user defines the data preprocessing function and the main job functions.
- The user reads out the CSV data, splits the data in smaller segments, such as a single day, and executes the code for each of these smaller segments. Executing once for a large block of data is inhibited by the quadratic scaling of the computations to account for correlations. As stated in the limitations section: short-circuiting in case of zero correlation or imposing a cutoff for low correlations to reduce computational complexity is currently not implemented.
- The user then simply loads the smaller segments of CSV data into the Job Handler and executes the predefined job.

#### Time series matching
When combining data from datasets with different temporal granularity, the code will try to perform a time harmonization. This means that the code will aggregate both timeseries to a timestep equal to the lowest common multiple of the involved timesteps in a calculation. Additionally, it will ensure that computations are performed with datapoints spanning the same time interval (i.e. a datapoint spanning 8:00-8:10 is not combined with data spanning 9:00-9:10).   
 
If a derived variable required a harmonization of its dependencies, information about this will be stored in this variable’s `harmonization_cache`, which contains a `TimeHarmonizationData` object for each of the variable’s dependencies. This object contains information on how to transform the dependency to the temporal granularity of the derived variable (the timestep increase factor, pruned edges, etc.).  

![Alt text](https://github.com/user-attachments/assets/6ad239fb-51ca-40a1-87e1-6db590520468)

#### Retroactive definition of basic variables
Some datasets only include the measured plane of array irradiance and not the pyranometer voltages and sensitivities. However, sometimes one wants to incorporate uncertainties in the voltages, sensitivities and POA irradiance separately. Users can leave basic variables empty as long as an equation is provide with which their values can be retroactively calculated from other variables, such as $V = G*S$ with $G$ and $S$ defined. Before execution, the code will populate the basic variable $V$ and then execute as if variable $V$ was defined regularly.
 
## Main classes and dataclasses – except engines
Most dataclasses are stored in the my_dataclasses.py file. The attributes of classes and dataclasses, including descriptions of what the attributes are, are all listed in the `__init__` and `__post_init__` functions of these classes.

#### The variable class
The variable class contains all information about a variable: values, equation, executables to calculate its values, dependencies, uncertainty, timedata and harmonization caches (we will return to this later) and metadata. Some of this data, such as uncertainty data, is stored in instances of other dataclasses that are owned by this variable.   
Variables have overloaded arithmetic and array_ufunc methods, such that one can easily perform regular and numpy calculations with variables. Be aware when using arithmetic on two timeseries-variables by a python-hardcoded statement such as `x = var_a + var_b`, that no automatic timeseries matching is performed in this case. The overloaded arithmetic only checks whether the arithmetic on the values arrays can be resolved, not whether the arrays span the same time interval.  
Variables further contain methods to report and plot their own values and uncertainties.

#### The VariableUncertainty dataclass
This dataclass contains all information related to the uncertainty of a variable. It is therefore always owned by a variable instance. Initially, the dataclass only contains a list of direct uncertainty sources acting on this variable. After uncertainty calculation this will contain information about the sources and magnitudes of all uncertainty at and downtree from this variable.

#### The UncertaintySource dataclass
This is a small dataclass that stores the characteristics of a single uncertainty source. It is also capable of constructing the uncertainty source’s temporal autocorrelation matrix.

#### The TimeHarmonizationData dataclass
This dataclass is stored to inform the code how the timeseries data of a specific dependency was changed when computing a variable’s values. This includes information on the timestep increase, how much data at the edges was discarded, etc.

#### The CSVData dataclass
A dataclass that is used to store CSV data. It is recommended to deprecate usage of this dataclass and move to a purely pandas-based framework.  
The dataclass contains a dictionary with data, and the common data timestep and time range. Additionally, this dataclass contains cleaning methods for this data, including checking for extreme, missing or negative values, interpolating replacements for these. In case the CSVData contains a column named `zenith`, one can also use functionality to compare locations of nan or nonzero values to the solar zenith angle, and detect missing data during the day or erroneous nonzero data during the night.


## The engines
The code makes use of 4 main engines. Engines act on a registry of variables and are designed to perform specific functions for the user.

#### The Equation engine
The equation engine takes an uninitialized set of variables and is converts it to a working equation tree. After input parsing, the variables only contain their equation in the form of a string. 
- The equation engine can read the equation string, extract dependencies from it, match these with the variables in the given (or an internal) variable registry and populate the variables with pointers to the variables they depend on, creating a recursive tree.
- The equation engine can check whether the equation tree is well-defined, in the sense that there are no missing dependencies and no circular definitions.
- The equation engine can read the equation string, converts it to a format sympy can work with and subsequently convert it to a sympy equation format. Sympy is a python library for symbolic mathematics, and can be used to create executables from equations and take partial derivatives.
- The equation engine handles the creation of executables and partial derivative executables from the equation for the user.

#### The Calculation engine
The calculation engine handles execution of equations, partial derivative equations and time aggregations. In short: it calculates values, but not uncertainties.   
The calculation engine can be told to calculate the values of a variable, and it will recursively ensure all dependencies are calculated. It will also make use of an internal time engine to ensure that dependencies are always time-harmonized before calculation. It is also used to calculate the values of the variable’s partial derivatives.

#### The Time engine
The time engine’s main functionality is to ensure dependencies are time-harmonized before executables are called. The variable’s executables simply take numeric or array-valued input to perform arithmetic, the time engine’s responsibility is to ensure that values that are given in these equation span the same time intervals.  
This means it has functionality to:
- check whether dependencies are time harmonious,
- rebin timeseries and prune ends if this is not the case,
- build TimeHarmonizationData objects for future reference, in case the harmonization must be repeated.  
Additionally, the time engine can be used to perform a hard temporal resolution decrease for a variable to irreversibly bring it to a coarser temporal granularity. Be aware that this procedure is irreversible and destructive: information will be lost. For relative uncertainty: calculating uncertainty and aggregating to a new granularity is different from aggregating a variable and calculating the aggregate’s uncertainty! The former procedure is always more precise.

#### The Uncertainty engine
The uncertainty engine is responsible for uncertainty calculation, propagation and aggregation. Recall that a variable’s uncertainty is defined by sources acting upon it directly (‘direct sources’) and sources acting on its dependencies (‘down-tree sources’). Down-tree sources are always multiplied by a sensitivity, which in first order equals the partial derivative of the variable with respect to the dependency this down-tree source acts on.  
When calculating uncertainty of a variable, we must obtain the uncertainty timeseries for each source acting on, or down-tree from, the desired variable. For absolute uncertainties these timeseries assume a single value, for relative these timeseries will be varying.
As stated previously: aggregating a timeseries and calculating the uncertainty of the aggregate, yields different results from calculating the uncertainty in the original temporal resolution and aggregating the uncertainty! Moreover, due to how correlation works, aggregating the values of an uncertainty to a final temporal granularity in two steps yields different results from aggregating it to this granularity in a single step.  
Therefore, the magnitude of the uncertainty due to a source is always kept in the source’s original temporal resolution. Aggregation to a new temporal granularity is only done when requested, in one step.   
When calculating the uncertainty of a variable, the engine will first retrieve the values of all uncertainty sources times their sensitivities, in their original temporal resolution. The sources are then individually temporally aggregated to the correct resolution, and then combined to arrive at a total uncertainty.  
The retrieval of the uncertainties per source is a recursive procedure that is depth-first, and works from the leaves upwards. So, when at a given node, the code will first retrieve the down-tree uncertainties times their sensitivities before adding the direct uncertainties of the present node and passing them all upwards.   

Suppose a variable is of hourly resolution and the uncertainty source is of 1-minute resolution. This means the sensitivity to this source will have hourly resolution. Then each group of 60 uncertainty values corresponding to each hour will be multiplied by the same sensitivity.  
The uncertainty engine keeps track of the relative resolution difference factors through an “upsample factors” registry. At each node all this information is stored inside the variable’s own VariableUncertainty dataclass. So, at each node the user has information on the sources acting on the variable, the values of non-aggregated uncertainties times the sensitivities, the upsample factors and the propagation path for each uncertainty source.



