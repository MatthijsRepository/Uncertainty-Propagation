from engines.my_dataclasses import Variable #, UncertaintySource
import numpy as np
from copy import deepcopy


class UncertaintyEngine:
    """ 
    This engine handles the calculation of uncertainty.
    
    The engine is typically invoked by calling the ``calculateTotalUncertainty`` function on a variable.
    The engine will then calculate and propagate all uncertainty downtree while properly handling time-aggregation.
    Engine also populates ``variable.uncertainty`` objects with (sub)results. Also makes use of these previously calculated results, if available.
    The uncertainty engine will invoke functionality of the equation, calculation and time engines wherever required.
    
    Attributes
    ----------
    equation_engine: EquationEngine
        Equation engine used to take partial derivatives and build executables on the fly.
    calculation_engine: CalculationEngine
        Calculation engine used to calculate the partial derivatives during uncertainty evaluation.
    
    Notes
    -----
    - Before ``calculateTotalUncertainty`` can be invoked on a derived (non-basic) variable, its values must have been calculated by a calculation engine.
    - Uncertainty timeseries are always kept in their root temporal resolution, even if they are combined with data of different temporal resolution.
      This is because partial aggregation of uncertainty timeseries is a destructive operation from an information perspective,
      aggregating a time-aggregation is mathematically incorrect. Time aggregations are therefore always calculated from root-resolution.
    - Uncertainties are allowed to be masked. This means that their uncertainties are set to zero whenever their parent variable is 0.
      This functionality is used to exclude uncertainty contributions from measurements at night.
      This setting only affects variables that are specified to be 'maskable' in the input script.
    """
    def __init__(self, equation_engine=None, calculation_engine=None):
        self.equation_engine = equation_engine
        self.calculation_engine = calculation_engine
    
    def _calculateUncertaintySourceValues(self, var, source, equation_engine=None):
        """ 
        Helper function that calculates and populates the uncertainty for a given uncertainty source.
        
        For relative uncertainty sources, uses the variables of the given variable to calculate the uncertainty.
        If the uncertainty source has an internal parent variable registered, checks whether this matches the given variable.
        In case the uncertainty source is defined through a more complex equation, through ``multiplier``, 
        the internal equation engine is used to build the required executable.
        
        Parameters
        ----------
        var: Variable
            Variable that owns the uncertainty source.
        source: UncertaintySource
            Uncertainty source for which to calculate the uncertainty values.
        
        Raises
        ------
        ValueError
            If the source has a registered parent variable that does not match the passed variable in ``var``.
        RuntimeError
            If the source is defined through a more complex ``multiplier`` equation that must be evaluated, 
            but the uncertainty engine has no equation engine to build the executable with.
        """
        #If the source has a parent variable given, check whether is matches the passed variable
        if source.parent_variable is not None and source.parent_variable is not var:
            raise ValueError(f"Error: passed variable {var.name} and registered parent variable {source.parent_variable.name} of uncertainty source {source.name} do not match. Cannot calculate uncertainties.")
        
        if source.is_relative:
            source.values = source.sigma * var.values
        else:
            source.values = source.sigma
        
        #This block handles complex uncertainty sources with a multiplier.
        if source.multiplier is not None:
            #Check if there is already an executable in place. If not: we create it now
            if source.executable is None:
                
                if equation_engine is None:
                    equation_engine = self.equation_engine
                if equation_engine is None:
                    raise RuntimeError(f"Cannot calculate the uncertainty values of uncertainty source {source.name} of variable {var.name}. \
                                     Please provide the uncertainty engine with an equation engine to interpret the source equation.")
                #Preparing source dependencies
                source.equation = source.multiplier
                source.dependencies = {}
                equation_engine.populateVariableDependencyNames(source)
                equation_engine.populateVariableDependencies(source)
                
                #Preparing and executing equation
                equation_engine.buildVariableExecutable(source)
            #Execute executable
            args = source.dependencies.values()
            rescale_values = source.executable(*args)
            
            #Rescaling the values by the multiplier values
            source.values *= rescale_values
        source.values = abs(source.values)
        
    def _prepareVariableDirectUncertainties(self, var):
        """ 
        Helper function that ensures all direct uncertainty sources acting on ``var`` have their values calculated.
        
        Parameters
        ----------
        var: Variable
            Variable for which the direct uncertainties must be calculated.
            
        Raises
        ------
        ValueError
            If no values are defined or calculated for the passed variable.
        """
        if var.values is None:
            raise ValueError(f"Cannot prepare uncertainty for variable {var.name}, please evaluate the variable itself first!")
        if var.uncertainty.direct_uncertainties_calculated is True:
            return
        #Set is_calculated flag to true
        var.uncertainty.direct_uncertainties_calculated = True
        #Calculate magnitudes of all direct uncertainties if there are any
        if len(var.uncertainty.direct_uncertainty_sources)==0:
            return
        else:
            for source in var.uncertainty.direct_uncertainty_sources:
                source.parent_variable = var
                self._calculateUncertaintySourceValues(var, source)
            
    def prepareAllDirectUncertainties(self, variables):
        """ 
        Prepares the direct uncertainties for all variables in the given variable set.
        
        Parameters
        ----------
        variables: list
            Variables for which the direct uncertainties must be calculated.
        
        """
        for var in variables:
            self._prepareVariableDirectUncertainties(var)
    
    def prepareDownTreeDirectUncertainties(self, var):
        """ 
        Recursively ensures all direct uncertainties acting on ``var``, and all variables downtree, are calculated.
        Recursion occurs depth-first and stops traversing a node when a variable has no dependencies.
        
        Parameters
        ----------
        var: Variable
            Starting node of the recursive step.
        
        Notes
        -----
        - Recursion has no protection for cyclically defined equation trees. This will not cause problems for well-defined trees. 
          Upon initialization, the equation engine checks if the tree is well-defined.
        """
        self._prepareVariableDirectUncertainties(var)
        if var.is_basic:
            return
        else:
            for dep in var.dependencies.values():
                self.prepareDownTreeDirectUncertainties(dep)

    def _getDependencyPartialsValues(self, var, equation_engine=None, calculation_engine=None):
        """ 
        Get the absolute values of a variable's partial derivatives with respect to each dependency.
        
        Uses an equation engine to ensure the executables of all partial derivatives are in place for the passed variable.
        Uses a calculation engine to execute all partial derivatives.
        
        Parameters
        ----------
        var: Variable
            Variable for which the partial derivatives are to be calculated.
        equation_engine: EquationEngine, optional
            Equation engine used to build partial derivative executables. Required if uncertainty engine has no internal equation engine.
        calculation_engine: CalculationEngine, optional
            Calculation engine used to evaluate all partial derivatives. Required if uncertainty engine has no internal calculation engine.
        
        Raises
        ------
        RuntimeError
            If no equation or calculation engines are owned by the uncertainty engine, or given in the function call.
        
        Returns
        -------
        dict[str, float or array_like]
            Dictionary containing as keys the variables to which the partial derivatives are taken,
            and as values the float or array_like calculated values of the partial derivatives.
        """
        if equation_engine is None:
            equation_engine = self.equation_engine
            if equation_engine is None:
                raise RuntimeError("Please provide an equation engine when calling this function.")
        if calculation_engine is None:
            calculation_engine = self.calculation_engine
            if calculation_engine is None:
                raise RuntimeError("Please provide a calculation engine when calling this function.")
        #Populate variable partial executables
        equation_engine.buildPartialDerivativeExecutables(var)
        partials_dict = calculation_engine.executeAllPartials(var, absolute_values=True, store_results=False, force_recalculation=False)
        return partials_dict
    
    def _convertNestedListTo2DArray(self, lst, forced_length=None):
        """ 
        Deprecated: Converts a nested list of various sizes to a 2D array block - extends scalars to the length of the rest of the array.
        
        Parameters
        ----------
        lst: list
            2D nested list to be converted to a square array.
        forced_length: int or None, default = None
            Size of the output matrix. If ``None``, the output matrix will have the size of the longest list in ``lst``.
        
        Returns
        -------
        np.ndarray[dtype=float]
            2D array conversion of the nested input list, with padded values where necessary.
        """
        if forced_length is None:
            forced_length = max((np.size(v) if np.ndim(v)>0 else 1) for v in lst)
        lst = [np.full(forced_length, v) if (np.isscalar(v) or (isinstance(v, np.ndarray) and v.size == 1))  else v for v in lst]
        return np.vstack(lst).astype(float)

    def _initializeUncertaintyPropagation(self, var, mask):
        """ 
        Short-circuits uncertainty calculation if the variable has been flagged as certain, or if its uncertainties have been previously calculated.
        Helper for ``UncertaintyEngine.getWeightedRootUncertainties`` function.
        
        Parameters
        ----------
        var: Variable
            Variable for which uncertainty data retrieval has been called.
        mask: bool
            Whether uncertainty is allowed to be masked. 
            Previously calculated uncertainties only match the requested data if they were calculated with the same masking setting.
        
        Returns
        tuple
            If the helper determines that previously computed results can be used.
            For more info on this tuple, see ``UncertaintyEngine.getWeightedRootUncertainties``
            Tuple has structure
            (
                list[UncertaintySource],
                list[ np.ndarray ],
                list[int],
                list[ list[int] ],
                list[ list[Variable] ]
            )
        None
            If the helper determines that the uncertainties must be (re)calculated by the engine.
            
        Notes
        -----
        - Recommended to replace the uncertainty data tuple with a dedicated dataclass, with integrated
          routines for merging, adding recursion layers and retrieving human-readable recursion stacks.
        """
        if var.uncertainty.is_certain:
            return [], [], [], [], []
        if var.uncertainty.total_uncertainty_calculated and var.uncertainty.root_weighted_uncertainties is not None:
            #We can only return previously calculated results if they both have the same mask settings, otherwise we have to recalculate
            if var.uncertainty.is_masked is mask:
                return (deepcopy(var.uncertainty.root_sources),
                        deepcopy(var.uncertainty.root_weighted_uncertainties),
                        deepcopy(var.uncertainty.root_total_upsample_factors),
                        deepcopy(var.uncertainty.root_local_upsample_factors),
                        deepcopy(var.uncertainty.root_propagation_paths))
        return None

    def _rootWeightedUncertaintyCalculator(self, var, dep_name, new_sensitivities, dep_weighted_uncertainties, 
                                           total_upsample_factors, local_upsample_factors):
        """ 
        Applies sensitivity to uncertainty while accounting for temporal resolution differences, and updates the temporal upsample factors.
        
        For a given variable and dependency, this function prunes the weighted uncertainties to the right length such that start and end times match.
        Multiplies sections of length 'total_upsample_factor' of the old sensitivities by its corresponding entry in the ``new_sensitivities``.
        In other words: if ``var`` has timesep of a factor ``total_upsample_factor=x`` greater than the timestep of the original uncertainty 
        then each set of x consecutive entries of this source's weighted uncertainty should be reweighted by the same factor in ``new_sensitivities``.
        Note that ``dep_weighted_uncertainties`` has the original temporal resolution of each uncertainty source, not necessarily that of the variable or dependency. 
        
        Parameters
        ----------
        var: Variable
            Variable that the uncertainty is being calculated for.
        dep_name: str
            Name of the dependency that the uncertainties are propagated upward from.
        new_sensitivities: dict[str, float or np.ndarray]
            Dictionary of dependency names with their respective partial derivatives of ``var``.
        dep_weighted_uncertainties: list[ np.ndarray ]
            The weighted uncertainties of the dependency.
        total_upsample_factors: list[int]
            The total upsample factors for each uncertainty compared to their original temporal resolution.
        local_upsample_factors: list[ list[int] ]
            The stack of upsample factors for each uncertainty at each recursion step.
            
        Returns
        -------
        tuple
            Updated weighted uncertainties, total and local upsample factors.
            Tuple has structure
            (
                list[ np.ndarray ],
                list[int],
                list[ list[int] ]
            )
        
        Notes
        -----
        - This function updates the weighted_uncertainties, total_upsample_factors and local_upsample_factors.
        """
        
        local_upsample_factor = 1
        #If a temporal resolution decrease was performed at the calculation of var, then we should prune the weighted uncertainties and update the total upsample factor accordingly
        if var.harmonization_cache is not None:
            harmonization_data = var.harmonization_cache[dep_name]
            for i in range(len(dep_weighted_uncertainties)):
                #Get low and high offset indices
                low_index, high_index = harmonization_data.getTotalOffsetSteps()
                #First we prune the weighted uncertainties using the original total upsample factor
                dep_weighted_uncertainties[i] = dep_weighted_uncertainties[i][int(low_index * total_upsample_factors[i]) : \
                                                                              int(high_index * total_upsample_factors[i]) ]
                #now we update the total upsample factor
                local_upsample_factors[i] += [harmonization_data.upsample_factor]
                total_upsample_factors[i] *= harmonization_data.upsample_factor
                #Here we extend the new sensitivities by copying each entry 'total_upsample_factors[i]' times
                shaped_new_sensitivities = np.kron(new_sensitivities[dep_name], np.ones(int(total_upsample_factors[i])))
                #now the two arrays are both of the same shape and we can multiply the two arrays directly
                dep_weighted_uncertainties[i] = dep_weighted_uncertainties[i] * shaped_new_sensitivities
        else:
            for i in range(len(dep_weighted_uncertainties)):
                #If no time harmonization was performed we can append a local upsample factor of 1 to the stacks
                local_upsample_factors[i] += [1]
                #Here we extend the new sensitivities by copying each entry 'total_upsample_factors[i]' times
                shaped_new_sensitivities = np.kron(new_sensitivities[dep_name], np.ones(total_upsample_factors[i]))
                #now the two arrays are both of the same shape and we can multiply the two arrays directly
                dep_weighted_uncertainties[i] = dep_weighted_uncertainties[i] * shaped_new_sensitivities
        return dep_weighted_uncertainties, total_upsample_factors, local_upsample_factors

    
    def _handleDirectUncertaintyData(self, var, mask):
        """ 
        Initializes the uncertainty stack for the direct uncertainty sources acting on ``var``. 
        Helper to the ``getWeightedRootUncertainties`` function.
        
        Parameters
        ----------
        var: Variable
            Variable for which to initialize an uncertainty stack.
        mask: bool
            Boolean indicating whether uncertainties should be masked if their parent variable is zero.
        
        Returns
        -------
        tuple
            Initialized uncertainty stack for the direct uncertainties acting on ``var``.
            For more info on this tuple, see ``UncertaintyEngine.getWeightedRootUncertainties``
            Tuple has structure
            (
                list[UncertaintySource],
                list[ np.ndarray ],
                list[int],
                list[ list[int] ],
                list[ list[Variable] ]
            )
        
        Notes
        -----
        - Recommended to replace the uncertainty data tuple with a dedicated dataclass, with integrated
          routines for merging, adding recursion layers and retrieving human-readable recursion stacks.
        """
        #Initialize the length of our timeseries
        n_values = 1 if isinstance(var.values, (float,int)) else len(var.values)
        
        uncertainty_mask = np.ones(n_values)
        if mask and var.is_maskable and not np.isscalar(var.values):
            uncertainty_mask = np.zeros(n_values)
            uncertainty_mask[np.nonzero(var.values)] = 1

        #Initialize the relevant objects using the direct uncertainty sources of this variable
        all_sources, all_weighted_uncertainties, all_total_upsample_factors, all_local_upsample_factors, all_propagation_paths = [], [], [], [], []
        for source in var.uncertainty.direct_uncertainty_sources:
            all_sources                 += [source]
            all_weighted_uncertainties  += [source.values * uncertainty_mask]
            all_total_upsample_factors  += [1]
            all_local_upsample_factors  += [[1]]
            all_propagation_paths       += [[var]]
        return all_sources, all_weighted_uncertainties, all_total_upsample_factors, all_local_upsample_factors, all_propagation_paths

    def getWeightedRootUncertainties(self, var, mask):
        """ 
        Recursively retrieves or calculates all uncertainties acting on or downtree from var, and passes them on as an uncertainty stack.
        
        Calculates direct uncertainties acting on ``var``.
        Retrieves all uncertainties downtree from ``var``, and ensures they are properly multiplied with their sensitivities.
        Applies time-aggregation if the present variable is a timesum, and resets the total upsample factors in this case.
        Recursion occurs depth-first and stops traversing a node when:
            - it encounters a basic variable.
            - it encounters a variable flagged as certain.
            - it encounters a variable with uncertainty already calculated with the same masking setting.
        
        Parameters
        ----------
        var: Variable
            Variable for which the uncertainty stack is requested.
        mask: bool
            Whether uncertainties to be retrieved should be masked if their parent variable is zero and maskable.
        
        Returns
        -------
        tuple
            Calculated uncertainty stack.
            Tuple has structure
            (
                list[UncertaintySource],
                list[ np.ndarray ],
                list[int],
                list[ list[int] ],
                list[ list[Variable] ]
            )
            - list of ``UncertaintySource`` objects acting on or downtree from ``var``.
            - list of arrays containing the weighted uncertainties (uncertainties multiplied by sensitivities) of each source.
            - list of integer total upsample factors: difference factor between timestep of ``var`` and uncertainty's root timestep.
            - list of lists containing the integer upsample factors along the propagation path, for each source.
            - list of lists containing the propagation path of each uncertainty source to their root variable.
            
        Notes
        -----
        - Recommended to replace the uncertainty data tuple with a dedicated dataclass, with integrated
          routines for merging, adding recursion layers and retrieving human-readable recursion stacks.
        """
        #Prepare the direct uncertainties acting on this variable
        if not var.uncertainty.direct_uncertainties_calculated:
            self._prepareVariableDirectUncertainties(var)
        
        #Initialization helper short-circuits the function in trivial cases (uncertainty already calculated, no uncertainty)
        args = self._initializeUncertaintyPropagation(var, mask)
        if args is not None:
            return args
        
        #Initialize containers
        all_sources, all_weighted_uncertainties, all_total_upsample_factors, all_local_upsample_factors, all_propagation_paths = [], [], [], [], []
        
        #First we retrieve any downtree uncertainties
        if not var.is_basic:
            partial_derivatives = self._getDependencyPartialsValues(var)
            for dep_name in var.dependency_names:
                #Recursively retrieve uncertainty data from the dependencies
                dep_sources, dep_weighted_uncertainties, dep_total_upsample_factors, dep_local_upsample_factors, dep_propagation_paths = \
                    self.getWeightedRootUncertainties(var.dependencies[dep_name], mask=mask)
                #If there are no uncertainties for this dependency we skip it immediately
                if len(dep_sources)==0:
                    continue
                
                #Update the sensitivities block-wise by blockwise-multiplying the previous weighted uncertainties with new sensitivities
                dep_weighted_uncertainties, dep_total_upsample_factors, dep_local_upsample_factors = \
                    self._rootWeightedUncertaintyCalculator(var, dep_name, partial_derivatives,
                                                            dep_weighted_uncertainties, dep_total_upsample_factors,
                                                            dep_local_upsample_factors)
                
                #Update propagation paths
                for path in dep_propagation_paths:
                    path += [var]
                #Append to the relevant containers
                all_sources                += dep_sources
                all_weighted_uncertainties += dep_weighted_uncertainties
                all_total_upsample_factors += dep_total_upsample_factors
                all_local_upsample_factors += dep_local_upsample_factors
                all_propagation_paths      += dep_propagation_paths
        
        #Now we append the direct uncertainty sources to our containers
        dir_sources, dir_weighted_uncertainties, dir_total_upsample_factors, dir_local_upsample_factors, dir_propagation_paths = \
            self._handleDirectUncertaintyData(var, mask=mask)
        all_sources                += dir_sources
        all_weighted_uncertainties += dir_weighted_uncertainties
        all_total_upsample_factors += dir_total_upsample_factors
        all_local_upsample_factors += dir_local_upsample_factors
        all_propagation_paths      += dir_propagation_paths
        
        #In case the variable is a timesum we are at a destructive node in our equation tree.
        #Therefore we must pass the timesummed root uncertainties here and reset the upsample factors
        if var.is_timesum:
            all_weighted_uncertainties = self.timeSumWeightedRootUncertainties(var, all_sources, all_weighted_uncertainties,
                                                                               all_local_upsample_factors, all_propagation_paths)
            all_total_upsample_factors = [1 for _ in all_total_upsample_factors]
        #Pass on the package
        return all_sources, all_weighted_uncertainties, all_total_upsample_factors, all_local_upsample_factors, all_propagation_paths
    
    def _calculateUncertaintyAggregation(self, weighted_uncertainties, source):
        """ 
        Function performing the aggregation of uncertainty timeseries data using source autocorrelation.
        
        If source correlation is trivial (0 or 1), short-circuits calculation to avoid matrix multiplications.
        Aggregation occurs along last axis of the ``weighted_uncertainties`` array. 
        For full aggregation, simply pass a 1D array. For partial aggregation, ensure 2D array of desired shape beforehand.
        
        Parameters
        ----------
        weighted_uncertainties: np.ndarray
            Array containing uncertainty timeseries data.
        source: UncertaintySource
            Uncertainty source corresponding to the uncertainty timeseries.
        
        Returns
        -------
        np.ndarray or float
            Array of 1 dimension less than ``weighted_uncertainties``. Aggregated weighted uncertainty array.
            If ``weighted_uncertainties`` is 1D, returns a scalar.
        """
        #Check for trivial instances, 0 or 1 autocorrelation
        if isinstance(source.correlation, (float, int)):
            if source.correlation == 0:
                return np.sqrt(np.sum(weighted_uncertainties**2, axis=-1))
            if source.correlation == 1:
                return np.sum(weighted_uncertainties, axis=-1)
        #If case is not trivial, build correlation matrix of correct size and calculate the vector-matrix-vector product
        size = np.shape(weighted_uncertainties)[-1]
        corr_matrix = source.getCorrelationMatrix(size)
        return np.sqrt(np.vecdot(weighted_uncertainties, np.matvec(corr_matrix, weighted_uncertainties)))
    
    def timeSumWeightedRootUncertainties(self, calling_var, sources, weighted_uncertainties, local_upsample_factors, propagation_paths):
        """ 
        Calculates aggregation factor and performs complete timesum of the uncertainty data.
        
        Timesums aggregate uncertainty timeseries according to source correlation and variable aggregation rules.
        Timesums convert rates to quantities by multiplying with the timestep (e.g. power will become energy).
        For each source, function parses backward through propagation path and applies a correction factor if upsampling occured for an intensive variable,
        see Notes for further explanation.
        Propagation path parsing stops either at the root, or if a previous timesum is encountered in the path.
        
        Parameters
        ----------
        calling_var: Variable
            Variable calling the timesum.
        sources: list[UncertaintySource]
            List containing all uncertainty sources acting on ``calling_var``.
        weighted_uncertainties: list[ np.ndarray ]
            Array containing uncertainty timeseries data.
        local_upsample_factors: list[ list[int] ]
            List of lists containing the integer upsample factors along the propagation path, for each source.
        propagation_paths: list[ list[Variable] ]
            List of lists containing the propagation path of each uncertainty source to their root variable.
        
        Returns
        -------
        list[ float ]
            The time-aggregated uncertainty due to each uncertainty source.
        
        Notes
        -----
        - Axiom: combining an intensive variable with an extensive variable results in an extensive variable, unless otherwise specified.
        - Time aggregation of intensive variables corresponds to taking the average.
        - Strictly speaking, one should do a dimensional analysis of the sensitivity it is multiplied with. However, as a rule of thumb,
          it is assumed that if the total is an extensive quantity and the variable is intensive, then the partial derivative with respect
          to the variable must still be extensive. Hence the product will be extensive. This logic also works the other way around.
        - If a partial aggregation occurs somewhere in the propagation path, we must account for this with a separate factor.
          Let 'A'  intensive and 'B' extensive, then 'A*B' is extensive and must be summed upon aggregation. 
          If 'B' has twice the timestep of 'A', then the uncertainty timeseries of 'A*B' is in the resolution of 'A'.
          We must therefore divide the reuslt by a factor 2, to account for the partial aggregation of 'A' to the timestep of 'B'.
        """
        new_weighted_uncertainties = []
        
        for i, source in enumerate(sources):
            #Calculate the correction factor for the aggregation of intensive variables and rates
            #Each upsampling by a factor f at the node of an intensive variable or a rate will add a factor 1/n to the total
            #Thus, we loop back through the stack (i.e. move downwards):
            aggregation_correction_factor = 1
            #We start our traversal one level below the calling variable
            for j, var in enumerate(reversed(propagation_paths[i][:-1])):
                #Timesums are destructive nodes, stop our propagation here
                if var.is_timesum:
                    break
                #If no aggregation rule is defined we skip this node
                if var.aggregation_rule is None:
                    continue
                if var.aggregation_rule.startswith("ave") or var.is_rate:
                    #Note: we take index j+1 here because we skipped the first variable in the propagation path
                    aggregation_correction_factor *= 1/local_upsample_factors[i][-(j+1)] 
            
            #Perform the time aggregation
            result = self._calculateUncertaintyAggregation(weighted_uncertainties[i], source)
            
            #Handle rules for the calling timesum
            if calling_var.aggregation_rule.startswith("ave"):
                aggregation_correction_factor *= 1/len(calling_var.non_aggregated_values)
            if calling_var.is_rate:
                aggregation_correction_factor *= calling_var.aggregation_step
            
            #Apply correction factor                
            result *= aggregation_correction_factor
            new_weighted_uncertainties.append(result)
        return new_weighted_uncertainties
                 
    def aggregateWeightedRootUncertainties(self, sources, weighted_uncertainties, total_upsample_factors, 
                                           local_upsample_factors, propagation_paths):
        """ 
        Performs time-aggregation on an uncertainty stack (set of root sources, weighted uncertainties and upsample factors),
        by applying time aggregation for all upsample factors in the ``local_upsample_factors`` list.
        
        Time aggregation is a destructive procedure: aggregating a time aggregation leads to incorrect results.
        
        Parameters
        ----------
        sources: list[UncertaintySource]
            List of ``UncertaintySource`` objects acting on or downtree from ``var``.
        weighted_uncertainties: list[np.ndarray]
            List of arrays containing the weighted uncertainties (uncertainties multiplied by sensitivities) of each source.
        total_upsample_factors: list[int]
            List of integer total upsample factors: difference factor between timestep of ``var`` and uncertainty's root timestep.
        local_upsample_factors: list[ list[int] ]
            List of lists containing the integer upsample factors along the propagation path, for each source.
        propagation_paths: list[ list[Variable] ]
            List of lists containing the propagation path of each uncertainty source to their root variable.
        
        Returns
        -------
        list[np.ndarray]
            List of arrays containing the uncertainty data in their new temporal resolution.
        
        Notes
        -----
        The notes on uncertainty aggregation of ``UncertaintyEngine.timeSumWeightedRootUncertainties`` apply here.
        """
        #Note, we cannot make this function fully numpy in general, because the arrays inside weighted_uncertainties can be of different lengths
        new_weighted_uncertainties = []
        
        for i, source in enumerate(sources):
            factor = total_upsample_factors[i]
            #We take a shortcut if the total upsample factor is 1, then no rebinning has to take place for this source
            if factor == 1:
                new_weighted_uncertainties.append(weighted_uncertainties[i])
                continue
            #Else: we rebin using the correlation matrix
            
            #Calculate the correction factor for the aggregation of intensive variables
            #Each upsampling by a factor f at the node of an intensive variable will add a factor 1/n to the total
            #Thus, we loop back through the stack:
            aggregation_correction_factor = 1
            for j, var in enumerate(reversed(propagation_paths[i])):
                #Timesums are destructive nodes, stop our propagation here
                if var.is_timesum:
                    break
                if var.aggregation_rule.startswith("ave") or var.is_rate:
                    #Note the difference with the timesum variant: there we skip the last variable in the backpropagation
                    #Here we do not skip the last step, hence we index with j
                    aggregation_correction_factor *= 1/local_upsample_factors[i][-j]
            
            #Calculate new uncertainties, append to list
            wu = weighted_uncertainties[i].reshape((-1, factor))
            result = self._calculateUncertaintyAggregation(wu, source)
            
            #Apply correction factor, append to lists
            result *= aggregation_correction_factor
            new_weighted_uncertainties.append(result)
            
        return new_weighted_uncertainties

    def calculateTotalUncertainty(self, var, mask=False):
        """ 
        Calculates the uncertainty stack (root weighted uncertainties, upsample factors, propagation paths) and total uncertainty for the given variable.
        Function requests uncertainty stack by calling the recursive ``UncertaintyEngine.getWeightedRootUncertainties`` function.
        Populates the uncertainty information of the called variable.
        
        Parameters
        ----------
        var: Variable
            Variable for which the uncertainty must be calculated.
        mask: bool, default=False
            Whether masking of uncertainties should be applied if the source's parent variable is zero, and is maskable.
            Intended to exclude uncertainty of e.g. zero readings during the night to be included in calculations.            
        
        Returns
        -------
        np.ndarray or float
            Calculated uncertainty (timeseries) of the variable values.
        """
        if var.uncertainty.total_uncertainty_calculated is True:
            return
        if not var.uncertainty.direct_uncertainties_calculated:
            self._prepareVariableDirectUncertainties(var)
        
        #Retrieve root sources, weighted uncertainties and upsample factors
        var.uncertainty.root_sources, var.uncertainty.root_weighted_uncertainties, \
            var.uncertainty.root_total_upsample_factors, var.uncertainty.root_local_upsample_factors, \
                var.uncertainty.root_propagation_paths = self.getWeightedRootUncertainties(var, mask=mask)
        
        #Populate masking setting
        var.uncertainty.is_masked = mask
        
        #Here we aggregate all uncertainties to the temporal resolution of the called variable
        aggregated_weighted_uncertainties = self.aggregateWeightedRootUncertainties(var.uncertainty.root_sources, var.uncertainty.root_weighted_uncertainties, 
                                                                                    var.uncertainty.root_total_upsample_factors, var.uncertainty.root_local_upsample_factors,
                                                                                    var.uncertainty.root_propagation_paths)
        aggregated_weighted_uncertainties = np.array(aggregated_weighted_uncertainties, dtype=float)
        
        #If there are no root sources we can break our calculation here
        if len(var.uncertainty.root_sources)==0:
            var.uncertainty.total_uncertainty_calculated = True
            var.uncertainty.is_certain = True
            return
        
        var.uncertainty.aggregated_weighted_uncertainties   = aggregated_weighted_uncertainties
        var.uncertainty.total_uncertainty                   = np.sqrt(np.sum(aggregated_weighted_uncertainties**2, axis=0))
        var.uncertainty.total_uncertainty_calculated        = True
        
        #If the uncertainty is just a single value, we replace the length-1 array by the numeric value
        if isinstance(var.uncertainty.total_uncertainty, np.ndarray) and len(var.uncertainty.total_uncertainty)==1:
            var.uncertainty.total_uncertainty = var.uncertainty.total_uncertainty[0]
        
        return var.uncertainty.total_uncertainty
    
    def calculateRootContributions(self, var):
        """ 
        Calculates the fractional contribution of each root source to the total variance in the variable.
        
        Parameters
        ----------
        var: Variable
            Variable for which the uncertainty contribution split is desired.
        
        Raises
        ------
        RuntimeError
            If the uncertainty of the variable has not yet been calculated.
        
        Returns
        -------
        np.ndarray or None
            An array or array of timeseries arrays of the fractional contributions to the total uncertainty of each source.
            Returns None if the variable has no uncertainty.
        """
        if not var.uncertainty.total_uncertainty_calculated:
            raise RuntimeError(f"Cannot split uncertainty of variable {var.name} to contributions. Please calculate the total uncertainty first.")
        
        if var.uncertainty.is_certain:
            return None
        
        result = np.divide(var.uncertainty.aggregated_weighted_uncertainties**2,
                           var.uncertainty.total_uncertainty**2,
                           out=np.zeros_like(var.uncertainty.aggregated_weighted_uncertainties),
                           where=(var.uncertainty.total_uncertainty != 0))
        return result
    
    def calculateRootContributions_EXCEL_METHOD(self,var):
        """ Split the total uncertainty between root contributions according to the method used in the ASTM-G213-17 excel spreadsheet.
            This method split the contribution of a source first between all other uncertainty sources acting on the source's parent variable
            this is subsequently multiplied by the contribution of the uncertainty in this parent variable to the total uncertainty in 'var' """
        root_variables = np.array([path[0].name for path in var.uncertainty.root_propagation_paths])
        unique_root_variables = list(set(root_variables))
        
        new_split = np.zeros(np.shape(var.uncertainty.aggregated_weighted_uncertainties), dtype=float)
        total_variable_uncertainties = {}
        
        
        total_sum = 0
        for unique_root in unique_root_variables:
            indices = np.where(root_variables == unique_root)
            #Calculate the total uncertainty times the sensitivity of all sources on this variable
            total_variable_uncertainties[unique_root] = np.sqrt(np.sum(var.uncertainty.aggregated_weighted_uncertainties[indices,:]**2, axis=1))[0]
            total_sum += total_variable_uncertainties[unique_root]
        
        for unique_root in unique_root_variables:
            #Get all indices in aggregated_weighted_uncertainty corresponding to this variable
            indices = np.where(root_variables == unique_root)
            #Calculate the total uncertainty times the sensitivity of all sources on this variable
            #total_w_variable_uncertainty = np.sqrt(np.sum(var.uncertainty.aggregated_weighted_uncertainties[indices,:]**2, axis=1))[0]
            #Calculate direct sum of weighted uncertainties from this variable
            summed_w_variable_uncertainty = np.sum(var.uncertainty.aggregated_weighted_uncertainties[indices,:], axis=1)[0]
            for i in indices[0]:
                numerator = var.uncertainty.aggregated_weighted_uncertainties[i] * total_variable_uncertainties[unique_root]
                denominator = summed_w_variable_uncertainty * total_sum
                new_split[i] = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=(denominator != 0))
        return new_split
    
    
    def plotRootContributions(self, var):
        if var.timestep is not None:
            self._plotTimeSeriesRootContributions(var)
            return
        #In case the variable is not a timeseries
        root_split = self.calculateRootContributions(var)
        
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        fig = plt.figure(figsize=(15,6), dpi=100)
        ax = plt.subplot(111)
        
        labels = [source.name for source in var.uncertainty.root_sources]
        bottom = 0
        for i, value in enumerate(root_split):
            ax.bar(0, value, bottom=bottom, label=labels[i])
            bottom += value
        
        ax.set_xlim(-1.5, 1.5)
        ax.get_xaxis().set_visible(False)
        ax.grid(axis='y')
        
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5))
        #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        plt.title(f"Contribution split between uncertainty sources of variable {var.name}")
        plt.ylabel("Percentage contribution split")
        plt.show()
        
    def _plotTimeSeriesRootContributions(self, var):
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        
        root_split = self.calculateRootContributions(var)
        
        
        time_axis = var.getTimeAxis()
        
        fig = plt.figure(figsize=(15,6), dpi=100)
        ax = plt.subplot(111)
        ax.stackplot(time_axis, *root_split, labels=[source.name for source in var.uncertainty.root_sources])
        
        ax.grid()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            
        # Put a legend to the right of the current axis
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5))
        
        plt.title(f"Contribution split between uncertainty sources of variable {var.name}")
        plt.xlabel("Time")
        plt.ylabel("Percentage contribution split")
        plt.show()
        
    def plotAbsoluteRootContributions(self, var, k=2, ylims=None, return_ax=False):
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        
        root_split = self.calculateRootContributions(var) * var.uncertainty.total_uncertainty * k
        
        time_axis = var.getTimeAxis()
        
        fig = plt.figure(figsize=(15,6), dpi=100)
        ax = plt.subplot(111)
        ax.stackplot(time_axis, *root_split, labels=[source.name for source in var.uncertainty.root_sources])
        
        if ylims is not None:
            ax.set_ylim(ylims[0], ylims[1])
        ax.grid()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            
        # Put a legend to the right of the current axis
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(0.848, 0.5))
            
        ax.set_title(f"Total uncertainty of variable {var.name}, k={k}")
        ax.set_xlabel("Time")
        ax.set_ylabel("Total uncertainty")
        if return_ax:
            return ax
        else:
            plt.show()
        
    def plotRelativeRootContributions(self, var, k=2, ylims=None, return_ax=False):
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        
        absolute_split = self.calculateRootContributions(var) * var.uncertainty.total_uncertainty * k
        root_split = np.divide(absolute_split,
                               var.values,
                               out=np.zeros_like(absolute_split),
                               where=(var.values != 0))
        root_split *= 100
        root_split[np.where(root_split>20)] = 20
        
        #root_split = self.calculateRootContributions(var) * var.uncertainty.total_uncertainty * k / var.values
        
        
        time_axis = var.getTimeAxis()
        
        fig = plt.figure(figsize=(15,6), dpi=100)
        ax = plt.subplot(111)
        ax.stackplot(time_axis, *root_split, labels=[source.name for source in var.uncertainty.root_sources])
        
        if ylims is not None:
            ax.set_ylim(ylims[0], ylims[1])
        ax.grid()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            
        # Put a legend to the right of the current axis
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5))
            
        ax.set_title( f"Uncertainty relative to total signal, k={k}")
        ax.set_xlabel("Time")
        ax.set_ylabel("Relative contribution split [%]")
        if return_ax:
            return ax
        else:
            plt.show()    
    
   
    
    def getWeightedRootUncertainties_OLD(self, var, store=False):
        """ Get all weighted root uncertainties of a variable. Specifically, for all uncertainty sources downtree,
            it returns an array of the weighted uncertainties with respect to the present variable, in the original temporal resolution of the uncertainty.
            Thus, if the temporal timestep of this variable is 10 times greater than that of a root uncertainty,
            each 10 consecutive entries in the weighted root uncertainty will be multiplied by the same top-level sensitivity.
            Returns a list of uncertainty source objects, a list of their corresponding weighted uncertainties, 
            and the total temporal resolution difference factor between the source and the present variable """
        if not var.uncertainty.direct_uncertainties_calculated:
            self._prepareVariableDirectUncertainties(var)
        
        if var.uncertainty.is_certain:
            return [], [], [], [], []
        
        if var.uncertainty.total_uncertainty_calculated and var.uncertainty.root_weighted_uncertainties is not None:
            return var.uncertainty.root_sources, var.uncertainty.root_weighted_uncertainties, var.uncertainty.root_total_upsample_factors, var.uncertainty.root_local_upsample_factors, var.uncertainty.root_propagation_paths
        
        #Initialize the length of our timeseries
        n_values = 1 if isinstance(var.values, (float,int)) else len(var.values)
        
        #Initialize the relevant objects using the direct uncertainty sources of this variable
        var.uncertainty.all_uncertainty_sources = var.uncertainty.direct_uncertainty_sources
        if len(var.uncertainty.direct_uncertainty_sources)>0:
            all_weighted_uncertainties = [source.values * np.ones(n_values) for source in var.uncertainty.direct_uncertainty_sources]
        else:
            all_weighted_uncertainties = []
        all_total_upsample_factors, all_local_upsample_factors, all_propagation_paths = [], [], []
        for _ in var.uncertainty.direct_uncertainty_sources:
            all_total_upsample_factors += [1]
            all_local_upsample_factors += [[1]]
            all_propagation_paths += [[var]]
        
        #If the variable is not a basic variable, we recursively retrieve all required data from the dependencies
        if not var.is_basic:
            new_sensitivities = self._getDependencyPartialsValues(var)
                        
            for dep_name in var.dependency_names:
                #Recursively retrieve uncertainty data from the dependencies
                dep_sources, dep_weighted_uncertainties, total_upsample_factors, local_upsample_factors, propagation_paths = self.getWeightedRootUncertainties(var.dependencies[dep_name])
                #If there are no uncertainties for this dependency we skip it immediately
                if len(dep_sources)==0:
                    continue
                
                #Update the sensitivities block-wise by blockwise-multiplying the previous weighted uncertainties with new sensitivities
                dep_weighted_uncertainties, total_upsample_factors, local_upsample_factors = self._rootWeightedUncertaintyCalculator(var, dep_name, new_sensitivities,
                                                                                                                                     dep_weighted_uncertainties, total_upsample_factors,
                                                                                                                                     local_upsample_factors)
                #Update propagation paths
                for path in propagation_paths:
                    path += [var]
                
                #Append to the relevant containers
                var.uncertainty.all_uncertainty_sources += dep_sources
                all_weighted_uncertainties += dep_weighted_uncertainties
                all_total_upsample_factors += total_upsample_factors
                all_local_upsample_factors += local_upsample_factors
                all_propagation_paths      += propagation_paths
        
            #In case the variable is a timesum we are at a destructive node in our equation tree.
            #and we must pass the timesummed root uncertainties here.
            if var.is_timesum:
                all_weighted_uncertainties = self.timeSumWeightedRootUncertainties(var.uncertainty.all_uncertainty_sources, all_weighted_uncertainties,
                                                                                   aggregation_rule=var.aggregation_rule)
                all_total_upsample_factors = [1 for _ in all_total_upsample_factors]

        return var.uncertainty.all_uncertainty_sources, all_weighted_uncertainties, all_total_upsample_factors, all_local_upsample_factors, all_propagation_paths
                
                
      
    

    
    
    
        


