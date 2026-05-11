## The calculation tool: what is it and what can you do with it?

The pythontool is a general uncertainty propagator that calculates uncertainty according to the Guide to the Uncertainty of Measurement, or GUM. 
The user specifies a set of variables in a text file, basic variables and derived variables.
- Basic variables are constants or measured variables
- Derived variables are specified by an equation and hence dependent on other variables. 
This way, the user defines an equation tree, or more strictly a directed acyclic graph. The variables are the nodes in this ‘tree’.

 ![Alt text](https://github.com/user-attachments/assets/89cbe44e-ae13-467f-9d0a-1dd895465d20)
 
In the textfile the user can specify the value of a variable, or express that the values are contained in a specific column of a pandas dataframe. If no value is specified it means the value must be calculated from other variables.  
The user can also specify for each variable:
- the aggregation rule (summing or averaging)
- whether the variable is a rate of a quantity over time, or whether it is simply a quantity
- whether uncertainty is ‘maskable’. Uncertainty masking means that the uncertainty of a time series is not included where time-series is zero-valued. For instance to exclude uncertainty in G at night when calculating aggregated uncertainty. Note that setting this to true allows the uncertainty to be masked, not that it is masked. Use of masking must be specified in the jobscript.  

The user can specify uncertainty sources for each variable. Multiple uncertainty sources can act on the same variable. For each uncertainty source the user can specify:
- absolute or relative error
- distribution
- deviation (mean assumed 0)
- autocorrelation over time (currently 0 and 1 are supported, extension to linear/exponential decay is planned)
- an optional multiplier to the magnitude, in the form of an equation. For instance, one can multiply the directional response error by a function of the solar zenith angle to more closely model its dependence on solar zenith angle.

## Points to keep in mind
This tool is meant to calculate the minimally achievable uncertainty in a quantity based on the specifications on the used measurement systems. In case the magnitude or characteristics of an uncertainty are unknown, it cannot be included in this calculation.  
Be aware that the results are always dependent on the used dataset. The uncertainty of a quantity $C = A * B$ with an uncertainty source in $B$, will always be dependent on the value(s) of $A$. To draw general conclusions, one needs to average or aggregate over a large amount of data.  
The tool can be used for ‘normal’ equations: regular arithmetic, mathematical operations, exponents, trigonometric functions and discrete integrals. Be aware that the code takes as the sensitivity coefficient for a source simply the instantaneous value of the partial derivative of the measurement equation with respect to the variable the uncertainty source acts on – which is a first-order approximation. The validity of this approximation may be questionable in case of highly nonlinear functions and large relative uncertainties.  
Also please keep in mind that the calculated uncertainty is only as good as the least accurate approximation used in the code. While the uncertainty propagation works with first-order approximations on the sensitivities, other approximations may be of even more significant influence. For instance, when estimating module temperature from ambient temperature, be aware that you are introducing another approximation to the model. When using such an approximation, it does not make sense to rigorously propagate the uncertainty in the measured wind speed to the final PR uncertainty, since the approximation itself may be a far greater source of (unquantified) uncertainty. Using such an approximation should be done with the intent of 'synthesizing' your own back-of-module temperature to extend the available dataset with. The temperature uncertainty should then be added to the synthesized back-of-module temperature measurement, not on the ambient temperature measurement.


## Assumptions, calculation errors and other limitations
- Importantly, the tool works under the assumption that separate uncertainty sources are independent (not cross-correlated).
- Additionally, at present no smart-matching feature is implemented in case the same uncertainty source acts on a variable multiple times. For instance, when two variables in an equation are both dependent on the same temperature measurement, the uncertainty in the temperature measurement will appear twice in the total uncertainty equation. The code will treat these as separate sources with their own sensitivities, while in practice the sensitivities should be combined, leading to an error.  
In mathematical terms: $u_T^2 * (s_1^2 + s_2^2) \neq u_T^2 (s_1 + s_2)^2$.  
This means that in such cases there will be a slight under-estimation of the uncertainty.  
To fix this issue, one can scan the var.uncertainty.root_sources list for duplicates to detect such cases, and implement exception handling in case duplicate sources are detected. This should be done when calculating the total uncertainty and when retrieving the root source contribution split. 
- A final limitation of the code is the amount of data that can be processed at a time. Because of the nature of correlated uncertainties, the computational complexity of aggregating uncertainty over time scales quadratically with the included timeframe. For trivial correlations (0 or 1), the code short-circuits to the analytical result to avoid matrix calculations. However, keep computational complexity in mind when working with slowly decaying temporal correlations.

It is recommended to further expand the code with decaying autocorrelation, with an inclusion-cutoff if the correlation is below a specified limit, leading to the correlation matrix being a band matrix. Correlated uncertainty aggregation over arbitrarily long timescales can then be performed iteratively and potentially parallelized for performance. 

## How the code works – general
The tool is developed in an object-oriented way. All variables are objects that store information on their own values, uncertainties, dependencies and properties such as calculation rules and state indicators.  
The operations on the variables, or on the equation tree as a whole, are performed by objects called engines. The calculation engine calculates variable values, the uncertainty engine calculates uncertainty, the time engine handles time matching between variables, et cetera.

#### User Interface
The `JobHandler` object is the main interface between the user and the internal functionality. The user can load equation trees and pandas dataframes objects to this handler and specify which tasks should be executed. The job handler will then perform pre-execution checks, variable initialization, job execution and post-job result storing. The user needs to specify the main job routine to the jobhandler by defining a `main` function executed when the `execute` function is called. The `JobHandler` object has wrapper functions for many of the main functionalities of the code that can be used inside the `main` function. Additionally, the JobHandler also has internal copies of all engines, so the user can also directly access the full engine functionality inside the `main` function with the right syntax. The jobscript.py file contains an illustration on how this works.   
To load timeseries data to variables, the JobHandler makes use of a custom data backend in which pandas dataframes are stored and that maps variables to the dataframes and columns containing their data. This datahandler also contains some rudimentary data cleaning and checking methods. 
Days can be blacklisted in the datahandler, such that the data is masked or execution fails if calculations involve a specific day.
Execution of the job is done through the `JobHandler.execute` function. Upon calling this function, the user can specify whether execution should be performed over the entire dataset or a subset of it, and how blacklisted days should be handled. Calling this function will cause the state of the equation tree to be reset, calculations to be performed and desired results, with some metadata, to be written to `JobHandler.Results`. The state of the equation tree will remain in place once execution is finished, for custom access and control.

#### Dependencies
The code is purely python-based and makes use mostly of standard python libraries: `numpy`, `matplotlib`, `pandas`. The code makes use of `SymPy` for the creation of executables of a variable's equation and to take symbolic partial derivatives of the variable's equation, which can subsequently be turned into executables, to calculate sensitivity coefficients. The code also uses `pvlib` to perform solar zenith angle calculations. In case you don't want to use this functionality, simply do not use functionality related to solar zenith angles and comment out the lines in `solar_module.py` related to it.

## Overall workflow
- The user defines an equation tree following the syntax of the example equation tree.
- The user creates a `JobHandler` instance.
- The user loads their data as pandas dataframes and passes these into the Job Handler.
- The user gives the equation tree text file to the Job Handler. The handler will compile an equation tree from this text file, check whether it is well-defined and non-circular. It will populate all variables with pointers to their dependencies, and it will prepare all variable equation executables.
- The user defines the `main` job function and potentially a preprocessing function if the pandas dataframes are not already preprocessed.
- The user executes the script by calling the `JobHandler.execute` function. In the function arguments, the user can specify the desired timerange of the execution and other execution parameters, such as the results group data is written to.
- After execution, written results can be accessed via dedicated Job Handler retrieval functions. The equation tree also retains its state and can be manually accessed through `JobHandler.variables`. 
 
The `main` function specifies which variables and uncertainties must be calculated, which data should be stored in the results storage and any further logic determining whether a result is considered successful.

#### Time series matching
When combining data from datasets with different temporal granularity, the code will try to perform a time harmonization. This means that the code will aggregate both timeseries to a timestep equal to the lowest common multiple of the involved timesteps in a calculation. Additionally, it will ensure that computations are performed with datapoints spanning the same time interval (i.e. a datapoint spanning 8:00-8:10 is not combined with data spanning 9:00-9:10).   
 
If a derived variable required a harmonization of its dependencies, information about this will be stored in this variable’s `harmonization_cache`, which contains a `TimeHarmonizationData` object for each of the variable’s dependencies. This object contains information on how to transform the dependency to the temporal granularity of the derived variable (the timestep increase factor, pruned edges, etc.).  

![Alt text](https://github.com/user-attachments/assets/6ad239fb-51ca-40a1-87e1-6db590520468)

#### Retroactive definition of basic variables
Some datasets only include the measured plane of array irradiance and not the pyranometer voltages and sensitivities. However, sometimes one wants to incorporate uncertainties in the voltages, sensitivities and POA irradiance separately. Users can leave basic variables empty as long as an equation is provided with which their values can be retroactively calculated from other variables, such as $V = G*S$ with $G$ and $S$ defined. Before execution, the code will populate the basic variable $V$ and then execute as if variable $V$ was defined regularly.

#### Handling of timesums
Timesums are effectively the time-integration of timeseries data. To start, we must note an important difference between timesums and aggregations. An aggregation is the rebinning of timeseries data to a coarser temporal resolution, and is usually done for automatic timeseries harmonization. When a variable is aggregated for a calculation, the uncertainty of the associated timeseries remains in their original temporal resolution – no details are lost in that regard.  
On the other hand, the timesum operation completely integrates the timeseries and associated uncertainties. The result is not a timeseries of values with timeseries of uncertainties, but a single value with single values for all uncertainties.  
Moreover, **the timesum operation converts rates to quantities**, while aggregation does not, even if the aggregation results in the timeseries being a single bin.  
 
The internal handling of a timeseries in an equation is special. If a timeseries is defined in an equation, the code will detect it and create a new variable for this timeseries. Thus, an equation for `PR` like `TS(‘Pout’ / ‘P0’) / TS(‘G’ / ‘G_STC’)` will have two dependencies, the internally created variables `TS(‘Pout’ / ‘P0’)` and `TS(‘G’ / ‘G_STC’)`. Similarly, even if a variable `A` is defined with an equation of just `TS(‘Pout’ / ‘P0’)`, the code will create a separate variable for this timesum, and the equation for `A` trivially refers to this variable.  
The reason behind this is for the code to be able to treat timesum calculations separately. In the above example for `PR`, the creation of intermediate variables helps during the evaluation of the equation for `PR`: instead of needing to internally resolve two timesum statements during calculation, instead the codes calculates the timesums first as dependencies, and then calculates the PR through a simple division.  
The created timesum variables have as equation the equation specified inside the timesum, and an `is_timesum` flag that is set to `True`. The calculation and uncertainty engines will first treat the equation inside the timesums regularly, check for this flag and then handle the final aggregation of the values if this `is_timesum` flag is set to `True`.  
 
This has further advantages. First of all, this allows the engines to make use of the regular framework for matching timeseries data of different temporal resolution. For instance, a timeseries `B` =  `TS(‘G’*’C_25’)` will first treat the equation `’G’*’C_25’` with time-matching, and then aggregate the final timeseries.  
A second advantage is that it allows to separately define aggregation rules. In general: the logic is that the intermediate timesum variables created by the equation engine will contain information on how to calculate the values of the timesum (aggregation rules, whether the summed quantity is a rate) while the variable above it, declared by the user, contains information on how to treat the result of the timesum in further computations.   

## Short API overview: main classes and dataclasses – except engines
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


## Short API overview: engines
The code makes use of 4 main engines. Engines act on a registry of variables and are designed to perform specific functions for the user.

#### The Equation engine
The equation engine takes an uninitialized set of variables and converts it to a working equation tree. After input parsing, the variables only contain their equation in the form of a string. 
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
Additionally, the time engine can be used to perform a hard temporal resolution decrease for a variable to irreversibly bring it to a coarser temporal granularity. Be aware that this procedure is irreversible and destructive: information will be lost.

#### The Uncertainty engine
The uncertainty engine is responsible for uncertainty calculation, propagation and aggregation. Recall that a variable’s uncertainty is defined by sources acting upon it directly (‘direct sources’) and sources acting on its dependencies (‘down-tree sources’). Down-tree sources are always multiplied by a sensitivity, which in first order equals the partial derivative of the variable with respect to the dependency this down-tree source acts on.  
When calculating uncertainty of a variable, we must obtain the uncertainty timeseries for each source acting on, or down-tree from, the desired variable. For absolute uncertainties these timeseries assume a single value, for relative these timeseries will be varying.
As stated previously: aggregating a timeseries and calculating the uncertainty of the aggregate, yields different results from calculating the uncertainty in the original temporal resolution and aggregating the uncertainty! Moreover, due to how correlation works, aggregating the values of an uncertainty to a final temporal granularity in two steps yields different results from aggregating it to this granularity in a single step.  
Therefore, the magnitude of the uncertainty due to a source is always kept in the source’s original temporal resolution. Aggregation to a new temporal granularity is only done when requested, in one step.   
When calculating the uncertainty of a variable, the engine will first retrieve the values of all uncertainty sources times their sensitivities, in their original temporal resolution. The sources are then individually temporally aggregated to the correct resolution, and then combined to arrive at a total uncertainty.  
The retrieval of the uncertainties per source is a recursive procedure that is depth-first, and works from the leaves upwards. So, when at a given node, the code will first retrieve the down-tree uncertainties times their sensitivities before adding the direct uncertainties of the present node and passing them all upwards.   

Suppose a variable is of hourly resolution and the uncertainty source is of 1-minute resolution. This means the sensitivity to this source will have hourly resolution. Then each group of 60 uncertainty values corresponding to each hour will be multiplied by the same sensitivity.  
The uncertainty engine keeps track of the relative resolution difference factors through an “upsample factors” registry. At each node all this information is stored inside the variable’s own VariableUncertainty dataclass. So, at each node the user has information on the sources acting on the variable, the values of non-aggregated uncertainties times the sensitivities, the upsample factors and the propagation path for each uncertainty source.


## Suggested changes and expansions
Handling of non-trivial autocorrelations. Most workflow is already in place.
- Adding builders for desired types of correlation matrix to `UncertaintySource`.
- Implementing correlation limits to keep computational complexity manageable.
- For calculations with large datasets: make efficeint use of sparseness of correlation matrix to keep calculations fast.
 
Implement an `UncertaintyPackage` dataclass for use in the `UncertaintyEngine`. Currently the uncertainty sources, weighted sensitivities, upsample factors and propagation paths are passed in tuple.
- Create a dedicated dataclass to replace this tuple with.
- Create dedicated routines to add a new propagation layer, and retrieve uncertainty data from this dataclass, instead of handling this in-engine. Cleaner separation and dedicated update pipelines reduce risk of potential mistakes.
- Implement usage of dataclass in the engine.
- Make the data class time-aware for each uncertainty dataseries it contains.
- Potentially create routine for calculating sensitivities from the weighed uncertainties by dividing these by their root direct uncertainty timeseries.
 
Implement dedicated routines for partial recalculations, specifically for recalculating and repropagating a single uncertainty source. This would support detailed investigations into a single uncertainty source, without other uncertainties needing to be recalculated.
- Recommended to implement `UncertaintyPackage` first.
 
Improve time-awareness of the equation tree. Currently, gaps in the timeseries are not allowed. This can lead to small gaps in data invalidating large timespans (e.g. a few minutes invalidate an entire day).
- Keeping calculations numpy-based is recommended, since automatic pandas time-index matching can complicate calculations with data of differing time resultions.
- Instead, perhaps allow `nan` data to be passed to the variables, but account for it during aggregations. Missing data could also be passed as 0, where appropriate.
 
Implement routine to write data to backend pandas dataframes, for complete input/output using pandas.
- By writing from variable state to data backend, the tool can be used to populate new columns in existing loaded dataframes, or create new dataframes with calculated data. 



