from my_dataclasses import Variable #, UncertaintySource
import numpy as np


class UncertaintyEngine:
    def __init__(self, variables, equation_engine=None, calculation_engine=None, time_engine=None):
        self.variables = variables
        self.equation_engine = equation_engine
        self.calculation_engine = calculation_engine
        self.time_engine = time_engine
    
    
    def _prepareVariableDirectUncertainties(self, var):
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
                if source.is_relative:
                    source.values = source.sigma * var.values
                else:
                    source.values = source.sigma
        
            
    def prepareAllDirectUncertainties(self, variables):
        for var in variables:
            self._prepareVariableDirectUncertainties(var)
    
    def prepareDownTreeDirectUncertainties(self, var):
        self._prepareVariableDirectUncertainties(var)
        if var.is_basic:
            return
        else:
            for dep in var.dependencies.values():
                self.prepareDownTreeDirectUncertainties(dep)
        
    
    def getDependencyUncertainties(self, var,recurse=True):
        root_sources = []
        root_sensitivities = {}
        for dep in var.dependencies.values():
            if not dep.uncertainty.total_uncertainty_calculated:
                if recurse:
                    self.calculateTotalUncertainty(dep)
                else:
                    raise ValueError(f"Cannot calculate uncertainty for variable {var.name}: dependency uncertainty of {dep.name} not calculated. Please enable recursion or do this manually.")
            #If the dependency has no uncertainty sources, we can skip appending to our registries
            if len(dep.uncertainty.all_uncertainty_sources) == 0:
                continue
            #Add all root sources from this dependency to our registry, and all root sensitivities to a new dictionary entry
            root_sources += dep.uncertainty.all_uncertainty_sources
            root_sensitivities[dep.name] = dep.uncertainty.all_uncertainty_sensitivities
        return root_sources, root_sensitivities
    
    
    
    def calculateTotalUncertainty(self, var, recurse=True):
        print(f"Calculating total uncertainty for variable {var.name}")
        if var.uncertainty.total_uncertainty_calculated is True:
            return
        if not var.uncertainty.direct_uncertainties_calculated:
            self._prepareVariableDirectUncertainties(var)
        
        n_values = 1 if isinstance(var.values, (float,int)) else len(var.values)
        if var.is_basic:
            var.uncertainty.all_uncertainty_sources = var.uncertainty.direct_uncertainty_sources
            var.uncertainty.all_uncertainty_sensitivities = np.ones((len(var.uncertainty.all_uncertainty_sources), n_values))
        else:
            root_sources, root_sensitivities = self.getDependencyUncertainties(var, recurse=recurse)
            var.uncertainty.all_uncertainty_sources = var.uncertainty.direct_uncertainty_sources + root_sources
            
            #Update sensitivity coefficients by multiplying old sensitivities with the new sensitivity coefficients
            new_sensitivities = self._getDependencyPartialsValues(var)
            totals = []
            print(f"Working for variable {var.name}, root sensitivities: {[np.shape(a) for a in root_sensitivities.values()]}")
            for dep_name in root_sensitivities.keys():
                totals.append(root_sensitivities[dep_name] * new_sensitivities[dep_name])
            new_sensitivities = self._convertNestedListTo2DArray(totals)
            
            var.uncertainty.all_uncertainty_sensitivities = np.append(np.ones((len(var.uncertainty.direct_uncertainty_sources), n_values)),\
                                                                      new_sensitivities,
                                                                      axis=0)
            print(np.shape(var.uncertainty.all_uncertainty_sensitivities))
        
        if len(var.uncertainty.all_uncertainty_sources) == 0:
            var.uncertainty.total_uncertainty_calculated = True
            var.uncertainty.is_certain = True
            print(f"Finished calculating total uncertainty for variable {var.name}")
            return
        
        all_uncertainty_magnitudes = [u.values for u in var.uncertainty.all_uncertainty_sources]
        all_uncertainty_magnitudes = self._convertNestedListTo2DArray(all_uncertainty_magnitudes)
        
        var.uncertainty.total_uncertainty = np.sqrt(np.sum( (all_uncertainty_magnitudes * var.uncertainty.all_uncertainty_sensitivities)**2 , axis=0) )
        var.uncertainty.total_uncertainty_calculated = True
        print(f"Finished calculating total uncertainty for variable {var.name}")
        return


    def partialAggregation(self, var, harmonization_data):
        """ Performs time-aggregation by temporal subspaces, defined by the harmonization_data TimeHarmonizationData dataclass instance """
        if not var.uncertainty.total_uncertainty_calculated:
            raise ValueError(f"Cannot aggregate uncertainty of variable {var.name} because its uncertainty is not yet evaluated.")
        
        low_index, high_index = harmonization_data.getTotalOffsetSteps()
        #If there the harmonization was only edge pruning and no change in temporal resolution, we simply return the pruned uncertainty
        if harmonization_data.upsample_factor == 1:
            return var.uncertainty.total_uncertainty[low_index:high_index]
        
        #If there was a change in temporal resolution, we must aggregate our uncertainty - while taking correlation into account
        new_length = int((harmonization_data.high_index - harmonization_data.low_index) / harmonization_data.upsample_factor)
        
        #Re-build the all_uncertainty_magnitudes block
        all_uncertainty_magnitudes = [u.values for u in var.uncertainty.all_uncertainty_sources]
        all_uncertainty_magnitudes = self._convertNestedListTo2DArray(all_uncertainty_magnitudes)
        
        #Multiply the magnitudes with their respective sensitivities
        all_magnitudes_sensitivities = all_uncertainty_magnitudes * var.uncertainty.all_uncertainty_sensitivities
        
        
        
        
        test_combi = np.sqrt(np.sum(all_magnitudes_sensitivities**2, axis=0))
        
        print("New method hard aggregation")
        test_mag_sens = all_magnitudes_sensitivities[:, low_index:high_index]
        contribs_per_source = []
        for i in range(len(var.uncertainty.all_uncertainty_sources)):
            contrib = test_mag_sens[i,600]**2 + test_mag_sens[i,601]**2 + test_mag_sens[i,600]*test_mag_sens[i,601]*2
            contribs_per_source.append(contrib)
        contribs_per_source = np.array(contribs_per_source)
        print(np.sqrt(np.sum(contribs_per_source)))
        
        
        
        
        
        #Prune ends, divide into the new segments
        all_magnitudes_sensitivities = all_magnitudes_sensitivities[:, low_index:high_index].reshape((-1, new_length, harmonization_data.upsample_factor))
        #Take root square sum for each bin
        
        test_combi = test_combi[low_index:high_index].reshape((new_length, harmonization_data.upsample_factor))
        ones = np.ones(( harmonization_data.upsample_factor, harmonization_data.upsample_factor))
        total_combi = np.matvec(ones, test_combi)
        total_combi = np.vecdot(test_combi, total_combi)
        total_combi = np.sqrt(total_combi)
        print("Old method results in the following uncertainty")
        print(total_combi * 60 * 2 / 1000)
        print(total_combi[300])
        print()
        
        
        correlation_matrices = np.array([u.getCorrelationMatrix(harmonization_data.upsample_factor) for u in var.uncertainty.all_uncertainty_sources])
        all_magnitudes_sensitivities = np.einsum('stf, sfg, stg -> st', all_magnitudes_sensitivities, correlation_matrices, all_magnitudes_sensitivities)
        all_magnitudes_sensitivities = np.sqrt(all_magnitudes_sensitivities)
        
        #all_magnitudes_sensitivities = np.sqrt(np.sum(all_magnitudes_sensitivities**2, axis=2))
        
        #Now we must identify the correction factor for each uncertainty source independently
        #Correction factor is based on aggregation rule and correlation of each source (aggregation rule is that of the variable the uncertainty acts on)
        correction_factors = np.ones(len(var.uncertainty.all_uncertainty_sources))
        
        
        
        ###!!!
        #aggregation_rules = np.array([source.parent_variable.aggregation_rule for source in var.uncertainty.all_uncertainty_sources])
        #correction_factors[np.where(aggregation_rules=="average")] = 1 / harmonization_data.upsample_factor
        if var.aggregation_rule == "average":
            correction_factors *= 1 / harmonization_data.upsample_factor
        ###!!! Ensure that correction factors are applied in the correct way! Once for each variance in the covariance!
        
        
        ###!!! NOTE: aggregation with correlation is implemented wrong rn, to fix tomorrow!
        
        #correlations = np.array([source.correlation for source in var.uncertainty.all_uncertainty_sources])
        #correction_factors[np.where(correlations == 0)] *= 1/np.sqrt(harmonization_data.upsample_factor)
        
        all_magnitudes_sensitivities = correction_factors[:, np.newaxis] * all_magnitudes_sensitivities
        
        return np.sqrt(np.sum( (all_magnitudes_sensitivities)**2 , axis=0) )
        
        
        
    
    
    
        
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
        
        correlations = {}
        #Loop through dependencies and retrieve correlations
        for dep_name in var.uncertainty.dependency_uncertainty_names:
            dep = var.dependencies[dep_name]
            #Check if uncertainty calculation is performed, optionally auto calculate this
            if not dep.uncertainty.is_calculated:
                if auto_calculate or force_recalculation:
                    self.calculateUncertainty(var, recurse=True)
                else:
                    raise ValueError(f"Cannot retrieve correlation of dependency {dep.name} of variable {var.name} - uncertainty not calculated yet. Please enable auto-calculation or call calculation manually.")
            
            #If no uncertainty: simply return 0 as correlation. This will not interfere with the actual uncertainty correlation: contribution of this dependency will be 0.
            if dep.uncertainty.is_certain: 
            #if dep.uncertainty.total_uncertainty == 0: ###!!!
                correlations[dep_name] = np.array([0])
                continue
            #Now we retrieve or recursively calculate the correlation of this dependency
            if (dep.uncertainty.correlation is None and recurse) or force_recalculation:
                cor = self.calculateCorrelation(dep, auto_calculate=auto_calculate, recurse=recurse, force_recalculation=force_recalculation)
            else:
                cor = dep.uncertainty.correlation
            correlations[dep_name] = cor
        return correlations
    
                
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
    
    
    def _harmonizeDependencyUncertainties(self, var, dependency_uncertainties, time_engine=None): ###!!! perhaps move to time engine
        """  """

        return 
    
    def _harmonizeDependencyCorrelations(self, var, dependency_correlations): ###!!! perhaps move to time engine
        
        return
    
    
    def _convertNestedListTo2DArray(self, lst):
        """ Converts a nested list of various sizes to a 2D array block - extends scalars to the length of the rest of the array """
        max_length = max((np.size(v) if np.ndim(v)>0 else 1) for v in lst)
        lst = [np.full(max_length, v) if (np.isscalar(v) or (isinstance(v, np.ndarray) and v.size == 1))  else v for v in lst]
        return np.vstack(lst).astype(float)
      
    
    
    
            
    def _correlationCalculationRequirementsHelper(self, var, auto_calculate, recurse, force_recalculation):
        """ check if necessary elements are already calculated, and optionally auto-calculate these """
        return
    
   
    def calculateCorrelation(self, var, auto_calculate=False, recurse=True, force_recalculation=False):
        """ Calculates the correlation of the uncertainty for each timestep """
        return
    

    

    def timeAggregateTotalUncertainty(self, var, auto_calculate=False):
        return 
    
    
    
        


