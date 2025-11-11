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


    def _initializeUncertaintyPropagation(self, var):
        if not var.uncertainty.direct_uncertainties_calculated:
            self._prepareVariableDirectUncertainties(var)
        if var.uncertainty.is_certain:
            return [], [], [], [], []
        if var.uncertainty.total_uncertainty_calculated and var.uncertainty.root_weighted_uncertainties is not None:
            return (var.uncertainty.root_sources,
                    var.uncertainty.root_weighted_uncertainties,
                    var.uncertainty.root_total_upsample_factors,
                    var.uncertainty.root_local_upsample_factors,
                    var.uncertainty.root_propagation_paths)
        return None


    def _rootWeightedUncertaintyCalculator(self, var, dep_name, new_sensitivities, dep_weighted_uncertainties, total_upsample_factors, local_upsample_factors):
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
                #First we prune the weighted uncertainties using the original total upsample factor
                dep_weighted_uncertainties[i] = dep_weighted_uncertainties[i][int(harmonization_data.low_index * total_upsample_factors[i]) : \
                                                                              int(harmonization_data.high_index * total_upsample_factors[i]) ]
                #now we update the total upsample factor
                local_upsample_factors[i] += [harmonization_data.upsample_factor]
                total_upsample_factors[i] *= harmonization_data.upsample_factor
                #Here we extend the new sensitivities by copying each entry 'total_upsample_factors[i]' times
                shaped_new_sensitivities = np.kron(new_sensitivities[dep_name], np.ones(int(total_upsample_factors[i])))
                #now the two arrays are both of the same shape and we can multiply the two arrays directly
                dep_weighted_uncertainties[i] = dep_weighted_uncertainties[i] * shaped_new_sensitivities
        else:
            for i in range(len(dep_weighted_uncertainties)):
                #Here we extend the new sensitivities by copying each entry 'total_upsample_factors[i]' times
                shaped_new_sensitivities = np.kron(new_sensitivities[dep_name], np.ones(total_upsample_factors[i]))
                #now the two arrays are both of the same shape and we can multiply the two arrays directly
                
                
                #print()
                #Current problem: the sensitivities are of length 1 when we are dealing with timesums, but we want the partials from the expression inside the timesum
                #print(f"SHape new sensitivities {new_sensitivities[dep_name]}")
                #print(f"Variable {var.name}, dependency {dep_name}, total factor {total_upsample_factors[i]}")
                #print(f"Shape current: {np.shape(dep_weighted_uncertainties[i])}, shape sensitivities {np.shape(shaped_new_sensitivities)}")
                dep_weighted_uncertainties[i] = dep_weighted_uncertainties[i] * shaped_new_sensitivities
                #if var.name == 'my_test':
                #    print(dep_weighted_uncertainties)
        return dep_weighted_uncertainties, total_upsample_factors, local_upsample_factors

    
    def _handleDirectUncertaintyData(self, var):
        #Initialize the length of our timeseries
        n_values = 1 if isinstance(var.values, (float,int)) else len(var.values)
        
        #Initialize the relevant objects using the direct uncertainty sources of this variable
        all_sources, all_weighted_uncertainties, all_total_upsample_factors, all_local_upsample_factors, all_propagation_paths = [], [], [], [], []
        for source in var.uncertainty.direct_uncertainty_sources:
            all_sources                 += [source]
            all_weighted_uncertainties  += [source.values * np.ones(n_values)]
            all_total_upsample_factors  += [1]
            all_local_upsample_factors  += [[1]]
            all_propagation_paths       += [[var]]
        return all_sources, all_weighted_uncertainties, all_total_upsample_factors, all_local_upsample_factors, all_propagation_paths

    def getWeightedRootUncertainties(self, var):
        """ Get all weighted root uncertainties of a variable. Specifically, for all uncertainty sources downtree,
            it returns an array of the weighted uncertainties with respect to the present variable, in the original temporal resolution of the uncertainty.
            Thus, if the temporal timestep of this variable is 10 times greater than that of a root uncertainty,
            each 10 consecutive entries in the weighted root uncertainty will be multiplied by the same top-level sensitivity.
            Returns a list of uncertainty source objects, a list of their corresponding weighted uncertainties, 
            and the total temporal resolution difference factor between the source and the present variable """
        
        #Initialization helper checks whether direct uncertainties are already calculated and short-circuits the function in trivial cases (already calculated, no uncertainty)
        self._initializeUncertaintyPropagation(var)
        
        #Initialize containers
        all_sources, all_weighted_uncertainties, all_total_upsample_factors, all_local_upsample_factors, all_propagation_paths = [], [], [], [], []
        
        #First we retrieve any downtree uncertainties
        if not var.is_basic:
            partial_derivatives = self._getDependencyPartialsValues(var)
            for dep_name in var.dependency_names:
                #Recursively retrieve uncertainty data from the dependencies
                dep_sources, dep_weighted_uncertainties, dep_total_upsample_factors, dep_local_upsample_factors, dep_propagation_paths = \
                    self.getWeightedRootUncertainties(var.dependencies[dep_name])
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
            self._handleDirectUncertaintyData(var)
        all_sources                += dir_sources
        all_weighted_uncertainties += dir_weighted_uncertainties
        all_total_upsample_factors += dir_total_upsample_factors
        all_local_upsample_factors += dir_local_upsample_factors
        all_propagation_paths      += dir_propagation_paths
        
        #In case the variable is a timesum we are at a destructive node in our equation tree.
        #and we must pass the timesummed root uncertainties here.
        if var.is_timesum:
            all_weighted_uncertainties = self.timeSumWeightedRootUncertainties(all_sources, all_weighted_uncertainties,
                                                                               aggregation_rule=var.aggregation_rule)
            all_total_upsample_factors = [1 for _ in all_total_upsample_factors]
        
        #Pass on the package
        return all_sources, all_weighted_uncertainties, all_total_upsample_factors, all_local_upsample_factors, all_propagation_paths
        

    
    
    def timeSumWeightedRootUncertainties(self, sources, weighted_uncertainties, aggregation_rule):
        new_weighted_uncertainties = []
        for i, source in enumerate(sources):
            corr_matrix = source.getCorrelationMatrix(len(weighted_uncertainties[i]))
            result = np.sqrt(np.vecdot(weighted_uncertainties[i], np.matvec(corr_matrix, weighted_uncertainties[i])))
            new_weighted_uncertainties.append(result)
        return new_weighted_uncertainties


                 
    def aggregateWeightedRootUncertainties(self, sources, weighted_uncertainties, total_upsample_factors, aggregation_rule):
        """ Performs a time-aggregation on a retrieved set of root sources, weighted root uncertainties and upsample factors
            In other words, this function brings the weighted root uncertainties for a variable retrieved by the getWeightedRootUncertainties function
            and aggregates all uncertainty arrays to the timestep of the variable.
            Note that the time aggregation is a destructive procedure in general; time aggregation of time aggregates only makes sense for fully (un)correlated error sources """
        #Note, we cannot make this function fully numpy in general, because the arrays inside weighted_uncertainties can be of different lengths
        new_weighted_uncertainties = []
        for i, source in enumerate(sources):
            factor = total_upsample_factors[i]
            #We take a shortcut if the total upsample factor is 1, then no rebinning has to take place for thsi source
            if factor == 1:
                new_weighted_uncertainties.append(weighted_uncertainties[i])
                continue
            #Else: we rebin using the correlation matrix
            corr_matrix = source.getCorrelationMatrix(factor)
            wu = weighted_uncertainties[i].reshape((-1, factor))
            #Calculate new uncertainties, append to list
            result = np.sqrt(np.vecdot(wu, np.matvec(corr_matrix, wu)))
            new_weighted_uncertainties.append(result)
            
        return new_weighted_uncertainties

    
    def calculateTotalUncertainty(self, var, recurse=True):
        if var.uncertainty.total_uncertainty_calculated is True:
            return
        if not var.uncertainty.direct_uncertainties_calculated:
            self._prepareVariableDirectUncertainties(var)
        
        #Retrieve root sources, weighted uncertainties and upsample factors
        var.uncertainty.root_sources, var.uncertainty.root_weighted_uncertainties, \
            var.uncertainty.root_total_upsample_factors, var.uncertainty.root_local_upsample_factors, \
                var.uncertainty.root_propagation_paths = self.getWeightedRootUncertainties(var)
            
        #Here we aggregate all uncertainties to the temporal resolution of the called variable
        aggregated_weighted_uncertainties = self.aggregateWeightedRootUncertainties(var.uncertainty.root_sources, var.uncertainty.root_weighted_uncertainties, 
                                                                                    var.uncertainty.root_total_upsample_factors, aggregation_rule=var.aggregation_rule)
        aggregated_weighted_uncertainties = np.array(aggregated_weighted_uncertainties, dtype=float)
        
        #If there are no root sources we can break our calculation here
        if len(var.uncertainty.root_sources)==0:
            var.uncertainty.total_uncertainty_calculated = True
            var.uncertainty.is_certain = True
            return
        
        var.uncertainty.aggregated_weighted_uncertainties = aggregated_weighted_uncertainties
        var.uncertainty.total_uncertainty = np.sqrt(np.sum(aggregated_weighted_uncertainties**2, axis=0))
        var.uncertainty.total_uncertainty_calculated = True
        
        
        for i, source in enumerate(var.uncertainty.root_sources):
            path = [a.name for a in var.uncertainty.root_propagation_paths[i]]
            print(f"{source.name} : {path}")
            print(var.uncertainty.root_local_upsample_factors[i])
            print(var.uncertainty.root_total_upsample_factors[i])
        return var.uncertainty.total_uncertainty
    
    
   

    
    
    
    
    
    
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
                
                
      
    

    
    
    
        


