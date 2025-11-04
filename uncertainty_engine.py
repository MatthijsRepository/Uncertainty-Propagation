from my_dataclasses import Variable #, UncertaintySource
import numpy as np


class UncertaintyEngine:
    def __init__(self, variables, equation_engine=None, calculation_engine=None, time_engine=None):
        self.variables = variables
        self.equation_engine = equation_engine
        self.calculation_engine = calculation_engine
        self.time_engine = time_engine
        
        
    def _calculateDirectUncertainty(self, var):
        """ Calculates the magnitudes of the direct uncertainties for each source and returns these in an array; or 0 if there is no uncertainty """
        if len(var.uncertainty.direct_uncertainty_sources)==0:
            return 0
        else:
            #Prepare uncertainty array
            if isinstance(var.values, np.ndarray):
                u_array = np.zeros((len(var.uncertainty.direct_uncertainty_sources), len(var.values)))
            else: #values is float or int
                u_array = np.zeros((len(var.uncertainty.direct_uncertainty_sources), 1))
            
            #Fill the uncertainty array
            for index, u in enumerate(var.uncertainty.direct_uncertainty_sources):
                #names_list.append(u.name) ###!!!
                if u.is_relative is False:
                    u_array[index] = u.sigma #for absolute uncertainty: cast value directly into array
                else:
                    u_array[index] = u.sigma * var.values
            return u_array
                    
    def _retrieveDependencyUncertainties(self, var, recurse):
        """ Retrieve total uncertainties for all dependencies, returns dictionary """
        if var.is_basic:
            return 0          

        #Initialize dictionary of dependencies and associated uncertainties
        total_dep_uncertainties = {}
    
        #Loop through dependencies and populate uncertainties dictionary
        for dep in var.dependencies.values():
            #Check if dependency has uncertainty sources and whether these are calculated already
            if not dep.uncertainty.is_calculated:
                if recurse:
                    self.calculateUncertainty(dep, recurse)
                else:
                    raise ValueError(f"Uncertainty of dependency {dep.name} of variable {var.name} not yet calculated, please ensure this or enable recursive calculation")
            total_dep_uncertainties[dep.name] = dep.uncertainty.total_uncertainty
        return total_dep_uncertainties
    
    def _retrieveDependencyCorrelations(self, var, auto_calculate, recurse, force_recalculation=False):
        """ Retrieve correlations of the dependencies; optionally auto-calculates necessary elements, recurses through tree or recalculates all values """
        if var.is_basic:
            return 0
        
        correlations = []
        #Loop through dependencies and retrieve correlations
        for dep in var.dependencies.values():
            #Check if uncertainty calculation is performed, optionally auto calculate this
            if not dep.uncertainty.is_calculated:
                if auto_calculate or force_recalculation:
                    self.calculateUncertainty(var, recurse=True)
                else:
                    raise ValueError(f"Cannot retrieve correlation of dependency {dep.name} of variable {var.name} - uncertainty not calculated yet. Please enable auto-calculation or call calculation manually.")
            
            #If no uncertainty: simply return 0 as correlation. This will not interfere with the actual uncertainty correlation: contribution of this dependency will be 0.
            if dep.uncertainty.is_certain: 
            #if dep.uncertainty.total_uncertainty == 0: ###!!!
                correlations.append(np.array([0]))
                continue
            #Now we retrieve or recursively calculate the correlation of this dependency
            if (dep.uncertainty.correlation is None and recurse) or force_recalculation:
                cor = self.calculateCorrelation(dep, auto_calculate=auto_calculate, recurse=recurse, force_recalculation=force_recalculation)
            else:
                cor = dep.uncertainty.correlation
            correlations.append(cor)
        return self._convertNestedListTo2DArray(correlations)
                
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
    
    def _convertNestedListTo2DArray(self, lst):
        """ Converts a nested list of various sizes to a 2D array block - extends scalars to the length of the rest of the array """
        max_length = max((np.size(v) if np.ndim(v)>0 else 1) for v in lst)
        lst = [np.full(max_length, v) if (np.isscalar(v) or (isinstance(v, np.ndarray) and v.size == 1))  else v for v in lst]
        return np.vstack(lst).astype(float)
      
    def calculateUncertainty(self, var, recurse=False):
        """ Calculates (the time series of) the total uncertainty of a variable 
            calculates total direct uncertainty from sources, calculates weighted contributions of dependencies from the equation
            combines uncertainties through squared sums """
        print(f"Calculating uncertainty of variable {var.name}")
        #Populate direct uncertainty array in variable uncertainty data
        var.uncertainty.direct_uncertainties = self._calculateDirectUncertainty(var)
        #Calculate total direct uncertainty squared (per timestep)
        total_direct_uncertainty_sq = np.sum(var.uncertainty.direct_uncertainties**2, axis=0)
        #Store total direct uncertainty in the uncertainty object for later reference
        var.uncertainty.total_direct_uncertainty = np.sqrt(total_direct_uncertainty_sq)
        
        #Calculate uncertainty due to dependencies
        if not var.is_basic:
            #dep_names, total_dep_uncertainties = self._retrieveDependencyUncertainties(var, recurse)
            total_dep_uncertainties = self._retrieveDependencyUncertainties(var, recurse)
            dep_partial_values = self._getDependencyPartialsValues(var)
            
            #Get a list of ordered dependency names, and an array of weighted contributions
            dep_names = dep_partial_values.keys()
            weighted_dep_uncertainties = [dep_partial_values[dep_name]*total_dep_uncertainties[dep_name] for dep_name in dep_names]
            #Convert list to array: we have to extend constants to fill the whole range of values in case we have time series
            weighted_dep_uncertainties = self._convertNestedListTo2DArray(weighted_dep_uncertainties)
            
            #Populate dependency uncertainty block of the variable uncertainty
            var.uncertainty.dependency_uncertainty_names = dep_names
            var.uncertainty.weighted_dependency_uncertainties = weighted_dep_uncertainties
            
            #Calculate total squared uncertainty due to dependencies
            total_dep_uncertainty_sq = np.sum(weighted_dep_uncertainties**2, axis=0)
            var.uncertainty.total_uncertainty = np.sqrt(total_direct_uncertainty_sq + total_dep_uncertainty_sq)
        else:
            var.uncertainty.total_uncertainty = var.uncertainty.total_direct_uncertainty
        
        #Flag the uncertainty to be calculated. We do this here to avoid errors in the correlation calculator in case we are dealing with a timesum
        var.uncertainty.is_calculated = True
        #Check whether the variable is certain, if so we flag this and short circuit to return
        if np.sum(var.uncertainty.total_uncertainty)==0:
            var.uncertainty.is_certain = True
            return var.uncertainty.total_uncertainty
        
        #If variable is a timesum: we need to aggregate uncertainty
        if var.is_timesum:
            #self.calculateCorrelation(var, auto_calculate=True, recurse=recurse)
            var.uncertainty.total_uncertainty = self.timeAggregateTotalUncertainty(var, auto_calculate=True)
        return var.uncertainty.total_uncertainty
    
    def splitDirectUncertaintyContributions(self, var):
        """ Splits the fractional contribution to the direct uncertainty between the direct uncertainty sources """
        if not var.uncertainty.is_calculated:
            raise ValueError(f"Cannot split uncertainty contributions for variable {var.name}, please calculate uncertainties first")
        
        #If the variable has no direct uncertainty sources we return immediately to avoid populating by NaN's
        if len(var.uncertainty.direct_uncertainty_sources)==0: ###!!! is this useful?
            return 
        
        #We divide each direct uncertainty by the total direct uncertainty (for each timestep) or set 0 if total direct uncertainty is 0 at that timestep
        total_direct_sum = np.sum(var.uncertainty.direct_uncertainties, axis=0)
        var.uncertainty.direct_uncertainties_contributions = np.divide(var.uncertainty.direct_uncertainties, 
                                                                       total_direct_sum, 
                                                                       out=np.zeros_like(var.uncertainty.direct_uncertainties), 
                                                                       where=(total_direct_sum != 0))
        return
    
    def splitTotalUncertaintyContributions(self, var):
        """ Splits the fractional contributions to the total uncertainty between direct sources and all dependencies """
        if not var.uncertainty.is_calculated:
            raise ValueError(f"Cannot split uncertainty contributions for variable {var.name}, please calculate uncertainties first")
        #First we calculate the sum of direct uncertainty and uncertainty from dependencies
        if var.uncertainty.weighted_dependency_uncertainties is not None:
            total_sum = np.sum(var.uncertainty.weighted_dependency_uncertainties, axis=0) + var.uncertainty.total_direct_uncertainty
        else:
            total_sum = var.uncertainty.total_direct_uncertainty
        
        #If no direct uncertainties are present we do not have to calculate this part
        if len(var.uncertainty.direct_uncertainty_sources)>0:
            #We calculate the contribution of the direct uncertainty sources, while protecting for divide by 0
            var.uncertainty.total_direct_uncertainty_contribution = np.divide(var.uncertainty.total_direct_uncertainty,
                                                                              total_sum,
                                                                              out=np.zeros_like(var.uncertainty.total_direct_uncertainty),
                                                                              where=(total_sum!=0))
        #We similarly calculate the contribution per dependency
        if var.uncertainty.weighted_dependency_uncertainties is not None:
            var.uncertainty.dependency_uncertainties_contributions = np.divide(var.uncertainty.weighted_dependency_uncertainties,
                                                                               total_sum,
                                                                               out=np.zeros_like(var.uncertainty.weighted_dependency_uncertainties),
                                                                               where=(total_sum!=0))
        return
    
    def splitToSourceContributions(self, var, force_recalculation=False):
        """ Splits uncertainty contributions down to all root uncertainty sources """
        if not var.uncertainty.is_calculated:
            raise ValueError(f"Cannot split uncertainty contributions for variable {var.name}, please calculate uncertainties first")
        
        #If the root split is already calculated and forced recalculation is not chosen, we simply return previously calculated values
        if var.uncertainty.root_uncertainty_contribution_split is not None and not force_recalculation:
            return var.uncertainty.root_uncertainty_contribution_split
        
        #If var is basic the root uncertainties are simply the direct uncertainties
        if var.is_basic:
            if var.uncertainty.direct_uncertainties_contributions is None:
                self.splitDirectUncertaintyContributions(var)
            return var.uncertainty.getSourceNames(), var.uncertainty.direct_uncertainties_contributions
        
        #If we are dealing with a derived variable the propagation is more complicated
        
        #Check if variable uncertainty split has already been performed, and if not: perform uncertainty splits
        if var.uncertainty.direct_uncertainties_contributions is None:
            self.splitDirectUncertaintyContributions(var)
        if var.uncertainty.total_direct_uncertainty_contribution is None or var.uncertainty.dependency_uncertainties_contributions is None:
            self.splitTotalUncertaintyContributions(var)
            
        sources = []
        root_contributions_scaled = []
        #Start with the direct contributions - we only do this if direct contributions are present
        if len(var.uncertainty.direct_uncertainty_sources)>0:
            sources.extend(var.uncertainty.getSourceNames())
            root_contributions_scaled.append(var.uncertainty.direct_uncertainties_contributions * var.uncertainty.total_direct_uncertainty_contribution)
        
        #Recursively retrieve root contributions of all dependencies and scale these to the relative dependency contribution
        for index, dep_name in enumerate(var.uncertainty.dependency_uncertainty_names):
            #Split each dependency recursively to source contributions!
            dep_sources, dep_root_split = self.splitToSourceContributions(var.dependencies[dep_name], force_recalculation=force_recalculation)
            #If there is no root uncertainty split: there is no uncertainty in this entire branch - continue
            if dep_root_split is None:
                continue
            #rescale the root split of the dependency to the relative weight of this dependency in the total uncertainty
            dep_root_split = np.multiply(dep_root_split, var.uncertainty.dependency_uncertainties_contributions[index])
            sources.extend(dep_sources)
            root_contributions_scaled.append(dep_root_split)
        
        #If the variable has no root uncertainty sources, we simply return None
        if len(sources)==0:
            return
        else:
            #Populate the variable uncertainty dataclass
            var.uncertainty.root_uncertainty_sources = sources
            var.uncertainty.root_uncertainty_contribution_split = self._convertNestedListTo2DArray(root_contributions_scaled)
            return var.uncertainty.root_uncertainty_sources, var.uncertainty.root_uncertainty_contribution_split
        
    
    
    def timeAggregateRootSplit(self, var, harmonization_data, force_recalculation=True):
        if not var.uncertainty.is_calculated:
            raise ValueError(f"Cannot split uncertainty contributions for variable {var.name}, please calculate uncertainties first")
            
        #if var.uncertainty.root_uncertainty_contribution_split is not None and not force_recalculation:
        #    return var.uncertainty.root_uncertainty_contribution_split
        
        new_length = int((harmonization_data.high_index - harmonization_data.low_index) / harmonization_data.upsample_factor)
        
        new_direct_uncertainties = np.zeros((len(var.uncertainty.direct_uncertainty_sources), new_length))
        
        #Get time aggregation of direct uncertainty sources
        for i, source in enumerate(var.uncertainty.direct_uncertainty_sources):
            u_source = var.uncertainty.direct_uncertainties[i, harmonization_data.low_index:harmonization_data.high_index].reshape((new_length, harmonization_data.upsample_factor))
            corr_matrix = np.ones((harmonization_data.upsample_factor, harmonization_data.upsample_factor)) * source.correlation
            np.fill_diagonal(corr_matrix, 1)
            
            u_new = np.tensordot(u_source, corr_matrix, axes=(1,0))
            u_new = np.vecdot(u_new, u_source)
            new_direct_uncertainties[i] = np.sqrt(u_new)
            
            
        #Dependencies
        #for i, dep_name in enumerate(var.uncertainty.dependency_uncertainty_names):
            #We check if the calculation of the present variable required time harmonization
            #If this is the case, a cache of TimeHarmonizationData objects should be present
            #If it is not present, we can return the regular uncertainty root split
            #if var.harmonization_cache is None:
                
            
            
            #Find correct time harmonization data object
            
            #recurse
            #multiply outcome by the dependency sensitivity
            #aggregate further, by source
        
        
        #if var.is_basic:
        #    if var.uncertainty.direct_uncertainties_contributions is None:
        #        self.splitDirectUncertaintyContributions(var)
        #    return var.uncertainty.getSourceNames(), var.uncertainty.direct_uncertainties_contributions
        
        
        
        #Get absolute values of the uncertainty split

        return
    
            
    def _correlationCalculationRequirementsHelper(self, var, auto_calculate, recurse, force_recalculation):
        """ check if necessary elements are already calculated, and optionally auto-calculate these """
        #If force_recalculation is True we simply recalculate all required elements here and return, otherwise we pass the regular checks
        if force_recalculation:
            self.calculateUncertainty(var, recurse=True)
            self.splitDirectUncertaintyContributions(var)
            self.splitTotalUncertaintyContributions(var)
            return
        
        #Pass regular checks, auto-calculate dependencies if wanted
        #Check if uncertainty already calculated or not
        if not var.uncertainty.is_calculated:
            if auto_calculate:
                self.calculateUncertainty(var, recurse=True)
            else:
                raise ValueError(f"Please calculate uncertainty for variable {var.name} before calculating correlation.")
        
        #Check if direct uncertainty source split is in place
        if len(var.uncertainty.direct_uncertainty_sources)>0:
            if var.uncertainty.direct_uncertainties_contributions is None:
                if auto_calculate:
                    self.splitDirectUncertaintyContributions(var)
                else:
                    raise ValueError(f"Please calculate depedency uncertainty contribution split for variable {var.name} before calculating correlation.")
            if var.uncertainty.total_direct_uncertainty_contribution is None:
                if auto_calculate:
                    self.splitTotalUncertaintyContributions(var)
                else:
                    raise ValueError(f"Please calculate total contribution split for variable {var.name} before calculating correlation.")
        
        #Check if dependency uncertainty split is in place
        if var.uncertainty.dependency_uncertainty_names is not None and var.uncertainty.dependency_uncertainties_contributions is None:
            if auto_calculate:
                self.splitTotalUncertaintyContributions(var)
            else:
                raise ValueError(f"Please calculate direct uncertainty contribution split for variable {var.name} before calculating correlation.")
        
    def calculateCorrelation(self, var, auto_calculate=False, recurse=True, force_recalculation=False):
        """ Calculates the correlation of the uncertainty for each timestep """
        #Check if correlation already calculated
        if var.uncertainty.correlation is not None and force_recalculation is False:
            return var.uncertainty.correlation
        
        #Check if all requirements for the calculation are in place, optionally auto-calculate these
        self._correlationCalculationRequirementsHelper(var, auto_calculate, recurse, force_recalculation)

        #Calculate correlation per timestep
        #Direct contributions
        if len(var.uncertainty.direct_uncertainty_sources)>0:
            direct_correlations = np.array([u.correlation for u in var.uncertainty.direct_uncertainty_sources])
            direct_correlation_contributions = np.multiply(direct_correlations[:,np.newaxis], var.uncertainty.direct_uncertainties_contributions)
            direct_correlation_contributions = np.multiply(direct_correlation_contributions, var.uncertainty.total_direct_uncertainty_contribution)
            direct_correlation_contributions = np.sum(direct_correlation_contributions, axis=0)
        else:
            direct_correlation_contributions = 0
        
        #Dependency contributions
        if var.uncertainty.dependency_uncertainties_contributions is not None:
            dependency_correlations = self._retrieveDependencyCorrelations(var, auto_calculate=auto_calculate, recurse=recurse, force_recalculation=force_recalculation)
            dependency_correlation_contributions = np.multiply(dependency_correlations, var.uncertainty.dependency_uncertainties_contributions)
            dependency_correlation_contributions = np.sum(dependency_correlation_contributions, axis=0)
        else:
            dependency_correlation_contributions = 0
        
        var.uncertainty.correlation = direct_correlation_contributions + dependency_correlation_contributions
        return var.uncertainty.correlation
    
    def timeAggregateRootSourceSplit(self, var, harmonization_data):
        #Splitting uncertainty time aggregate by root source only makes sense if we correctly time aggregate each root source with its own correlation
        
        #NOTE: think about the following as well:
            #Each root contribution is rescaled (perhaps even multiple times). Can we temporally integrate per root source and aggregate then? That is what is required.
            #I'd assumse so. Consider a single root source with a given correlation on its own. We don't see all the rescale factor right? We simply see that time series. We should be able to simply aggregate that.
        
        
        #Time-aggregating dependency uncertainty as a total is a destructive process
        
        
        return
    
    def partialAggregation(self, var, harmonization_data):
        """ Performs time-aggregation by temporal subspaces, defined by the harmonization_data TimeHarmonizationData dataclass instance """
        
        #Check if we have to include fractional parts
        #If so, we immediately fail
        if harmonization_data.low_fraction != 1:
            raise ValueError("Partial uncertainty aggregation of variable {var.name} failed: cannot handle fractional rebinning for uncertainty!")
            
        new_length = int((harmonization_data.high_index - harmonization_data.low_index) / harmonization_data.upsample_factor)
        print(new_length)
        
        uncertainty_segments = var.uncertainty.total_uncertainty[harmonization_data.low_index:harmonization_data.high_index].reshape((new_length, harmonization_data.upsample_factor))
        
        correlation_matrices = var.uncertainty.correlation[harmonization_data.low_index:harmonization_data.high_index].reshape((new_length, 1, harmonization_data.upsample_factor))
        correlation_matrices = np.repeat(correlation_matrices, repeats=harmonization_data.upsample_factor, axis=1)
        correlation_matrices[:,np.eye(harmonization_data.upsample_factor,dtype=bool)] = 1
        
        new_uncertainties = np.matvec(correlation_matrices, uncertainty_segments)   #Mv
        new_uncertainties = np.vecdot(uncertainty_segments, new_uncertainties)      #v^T (Mv)
        new_uncertainties = np.sqrt(new_uncertainties)
       
        new_correlations = np.average(var.uncertainty.correlation[harmonization_data.low_index:harmonization_data.high_index].reshape((new_length, harmonization_data.upsample_factor)), axis=1)
        ###!!! should maybe be a weighted average
        print('WARNING: at the moment we do incorrect aggregation of correlation')
        
        if var.aggregation_rule == "average":
            new_uncertainties /= harmonization_data.upsample_factor
        
        return new_uncertainties, new_correlations
    
    
    
    
    def timeAggregateTotalUncertainty_NEW(self, var):
        self.calculateCorrelation()
        return
    
    def timeAggregateTotalUncertainty(self, var, auto_calculate=False):
        if not var.uncertainty.is_calculated:
            if auto_calculate:
                self.calculateUncertainty(var, recurse=True) ##Recursion error !! ###!!!
            else:
                raise ValueError(f"Cannot calculate aggregated uncertainty for variable {var.name}: please calculate total uncertainty first!")
        if var.uncertainty.correlation is None:
            if auto_calculate:
                self.calculateCorrelation(var, auto_calculate=auto_calculate, recurse=True)
            else:
                raise ValueError(f"Cannot calculate aggregated uncertainty for variable {var.name}: please calculate correlation first!")
        if isinstance(var.uncertainty.total_uncertainty, float):
            raise ValueError(f"Cannot time aggregate uncertainty of variable {var.name} since its total uncertainty is not array-like but float.")
            
        #Construct correlation matrix
        corr_matrix = np.repeat(var.uncertainty.correlation[np.newaxis,:], repeats=len(var.uncertainty.correlation), axis=0)
        corr_matrix = (corr_matrix + corr_matrix.T) / 2
        np.fill_diagonal(corr_matrix, 1)
        
        #Perform calculation
        aggregate_uncertainty = np.sqrt(var.uncertainty.total_uncertainty @ corr_matrix @ var.uncertainty.total_uncertainty)
        
        #Aggregation rule implementation: ###!!!
        if var.aggregation_rule == "integrate": ###!!! OLD
            aggregate_uncertainty *= var.aggregation_step
        if var.aggregation_rule == "average":
            aggregate_uncertainty /= len(var.uncertainty.total_uncertainty)
        
        return aggregate_uncertainty
    
    
    
        


