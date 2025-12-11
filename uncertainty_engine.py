from my_dataclasses import Variable #, UncertaintySource
import numpy as np
from copy import deepcopy


class UncertaintyEngine:
    def __init__(self, variables, equation_engine=None, calculation_engine=None, time_engine=None):
        self.variables = variables
        self.equation_engine = equation_engine
        self.calculation_engine = calculation_engine
        self.time_engine = time_engine
    
    def _calculateUncertaintySourceValues(self, var, source):
        """ Calculates the uncertainty values for a given uncertainty source """
        if source.is_relative:
            source.values = source.sigma * var.values
        else:
            source.values = source.sigma
            
        if source.multiplier is not None:
            #Check if there is already an executable in place. If not: we create it now
            if source.executable is None:
                if self.equation_engine is None:
                    raise ValueError(f"Cannot calculate the uncertainty values of uncertainty source {source.name} of variable {var.name}. \
                                     Please provide the calculation engine with an uncertainty engine to interpret the source equation.")
                #Preparing source dependencies
                source.equation = source.multiplier
                source.dependencies = {}
                self.equation_engine.populateVariableDependencyNames(source)
                self.equation_engine.populateVariableDependencies(source)
                
                #Preparing and executing equation
                self.equation_engine.buildVariableExecutable(source)
            #Execute executable
            args = source.dependencies.values()
            rescale_values = source.executable(*args)
            
            #Rescaling the values by the multiplier values
            source.values *= rescale_values
        source.values = abs(source.values)
        
    def _prepareVariableDirectUncertainties(self, var):
        """ Prepares all direct uncertainties acting on variable var. """
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
        """ Prepares the direct uncertainties for all variables in the given variable set """
        for var in variables:
            self._prepareVariableDirectUncertainties(var)
    
    def prepareDownTreeDirectUncertainties(self, var):
        """ Prepares the direct uncertainties acting on var, and all variables downtree from var """
        self._prepareVariableDirectUncertainties(var)
        if var.is_basic:
            return
        else:
            for dep in var.dependencies.values():
                self.prepareDownTreeDirectUncertainties(dep)

    def _getDependencyPartialsValues(self, var, equation_engine=None, calculation_engine=None):
        """ Executes partial derivative executable builder form the equation engine and executes calculation of partial values, 
            returns dictionary of partial derivative values per dependency """
        if equation_engine is None:
            equation_engine = self.equation_engine
        if calculation_engine is None:
            calculation_engine = self.calculation_engine
        #Populate variable partial executables
        equation_engine.buildPartialDerivativeExecutables(var)
        partials_dict = calculation_engine.executeAllPartials(var, absolute_values=True, store_results=False, force_recalculation=False)
        return partials_dict
    
    def _convertNestedListTo2DArray(self, lst, forced_length=None):
        """ Converts a nested list of various sizes to a 2D array block - extends scalars to the length of the rest of the array """
        if forced_length is None:
            forced_length = max((np.size(v) if np.ndim(v)>0 else 1) for v in lst)
        lst = [np.full(forced_length, v) if (np.isscalar(v) or (isinstance(v, np.ndarray) and v.size == 1))  else v for v in lst]
        return np.vstack(lst).astype(float)

    def _initializeUncertaintyPropagation(self, var, mask):
        """ initializer for the getWeightedRootUncertainties function, prepares direct uncertainties and helps to short circuit the main retriever in trivial cases """
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
        """ For a given variable and dependency, this function will prune the weighted uncertainties to the right length, and multiply sections
            of length 'total_upsample_factor' of the old sensitivities by its corresponding entry in the new_sensitivities
            in other words: if variable var has timesep of a factor 'total_upsample_factor=x' greater than the timestep of the original uncertainty 
            then each set of x consecutive entries of the source weighted uncertainty should be reweighted by the same factor in new_sensitivities 
            important to note is that dep_weighted_uncertainties has the original temporal resolution of the uncertainty source, not necessarily that of the variable or dependency """
        
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
        """ Handles the direct uncertainty sources for the getWeightedRootUncertainties function. """
        #Initialize the length of our timeseries
        n_values = 1 if isinstance(var.values, (float,int)) else len(var.values)
        
        uncertainty_mask = np.ones(n_values)
        if mask and var.is_maskable and not np.isscalar(var.values):
            nz = np.nonzero(var.values)[0]
            if len(nz)>0:
                uncertainty_mask[:nz[0]] = 0  
                uncertainty_mask[nz[-1]:] = 0
                
        
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
        """ Get all weighted root uncertainties of a variable. Specifically, for all uncertainty sources downtree,
            it returns an array of the weighted uncertainties with respect to the present variable, in the original temporal resolution of the uncertainty.
            Thus, if the temporal timestep of this variable is 10 times greater than that of a root uncertainty,
            each 10 consecutive entries in the weighted root uncertainty will be multiplied by the same top-level sensitivity.
            Returns a list of uncertainty source objects, a list of their corresponding weighted uncertainties, 
            and the total temporal resolution difference factor between the source and the present variable """
        
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
        
    def timeSumWeightedRootUncertainties(self, calling_var, sources, weighted_uncertainties, local_upsample_factors, propagation_paths):
        """ This function performs a full timesum of the uncertainty of all root sources while keeping it split by source
            Temporal autocorrelation is included. Cross-correlation between sources is not included - sources are assumed independent """
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
            corr_matrix = source.getCorrelationMatrix(len(weighted_uncertainties[i]))
            result = np.sqrt(np.vecdot(weighted_uncertainties[i], np.matvec(corr_matrix, weighted_uncertainties[i])))
            
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
        """ Performs a time-aggregation on a retrieved set of root sources, weighted root uncertainties and upsample factors
            In other words, this function brings the weighted root uncertainties for a variable retrieved by the getWeightedRootUncertainties function
            and aggregates all uncertainty arrays to the timestep of the variable.
            Note that the time aggregation is a destructive procedure in general; time aggregation of time aggregates only makes sense for fully (un)correlated error sources """
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
            
            #Building correlation matrix and uncertainty vector
            corr_matrix = source.getCorrelationMatrix(factor)
            wu = weighted_uncertainties[i].reshape((-1, factor))
            #Calculate new uncertainties, append to list
            result = np.sqrt(np.vecdot(wu, np.matvec(corr_matrix, wu)))
            result *= aggregation_correction_factor
            new_weighted_uncertainties.append(result)
            
        return new_weighted_uncertainties

    def calculateTotalUncertainty(self, var, recurse=True, mask=False):
        """ Calculates the total uncertainty, and uncertainty split per source, for the given variable
            The calculation populates the root uncertainties, propagation paths and upsample factors of all root ucnertainties downtree """
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
        """ Returns the fractional contribution of the variance of each root source to the total variance in the variable """
        if not var.uncertainty.total_uncertainty_calculated:
            raise ValueError(f"Cannot split uncertainty of variable {var.name} to contributions. Please calculate the total uncertainty first.")
        
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
        
    def plotAbsoluteRootContributions(self, var, k=1):
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        
        root_split = self.calculateRootContributions(var) * var.uncertainty.total_uncertainty * k
        
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
        ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(0.848, 0.5))
        
        plt.title(f"Total uncertainty of variable {var.name}, k={k}")
        plt.xlabel("Time")
        plt.ylabel("Total uncertainty")
        plt.show()
        
    def plotRelativeRootContributions(self, var, k=2):
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
        
        ax.set_ylim(0,10)
        ax.grid()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            
        # Put a legend to the right of the current axis
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5))
        
        plt.title(f"Uncertainty relative to total signal, k={k}")
        plt.xlabel("Time")
        plt.ylabel("Relative contribution split [%]")
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
                
                
      
    

    
    
    
        


