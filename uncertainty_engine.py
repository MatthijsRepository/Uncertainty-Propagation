from my_dataclasses import Variable, UncertaintySource
import numpy as np


class UncertaintyEngine:
    def __init__(self, variables, equation_engine=None, calculation_engine=None):
        self.variables = variables
        self.equation_engine = equation_engine
        self.calculation_engine = calculation_engine
        
     
        
    def _calculateDirectUncertainty(self, var):
        """ Calculates the magnitudes of the direct uncertainties for each source and returns these in an array; or 0 if there is no uncertainty """
        if len(var.uncertainty.direct_uncertainty_sources)==0:
            return 0
        else:
            #names_list = [] ###!!!
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
                    
    
    
    
    def _retrieveDependencyUncertainties_OLD(self, var, recurse):
        """ Retrieve total uncertainties for all dependencies, returns tuple(list of names, array of uncertainties) """
        if var.is_basic:
            return 0         

        #Initialize dictionary of dependencies and associated uncertainties
        dep_uncertainty_names = []
        total_dep_uncertainties = []
    
        #Loop through dependencies and populate uncertainties dictionary
        for dep in var.dependencies.values():
            #Check if dependency has uncertainty sources and whether these are calculated already
            if not dep.uncertainty.is_calculated:
                if recurse:
                    self.calculateUncertainty(dep, recurse)
                else:
                    raise ValueError(f"Uncertainty of dependency {dep.name} of variable {var.name} not yet calculated, please ensure this or enable recursive calculation")
            
            dep_uncertainty_names.append(dep.name)
            total_dep_uncertainties.append(dep.uncertainty.total_uncertainty)
        
        #Convert list of total uncertainties to array
        total_dep_uncertainties = np.stack(total_dep_uncertainties, axis=0)
        print(total_dep_uncertainties)
        return dep_uncertainty_names, total_dep_uncertainties
            
    
    
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

        
    def _getDependencyPartialsValues(self, var, equation_engine=None):
        """ Executes partial derivative executable builder form the equation engine and executes calculation of partial values, 
            returns dictionary of partial derivative values per dependency """
        if equation_engine is None:
            equation_engine = self.equation_engine
        #Populate variable partial executables
        equation_engine.buildPartialDerivativeExecutables(var)
        partials_dict = var.executeAllPartials(absolute_values=True, store_results=False, force_recalculation=False)
        return partials_dict
    
    
    def _convertNestedListTo2DArray(self, lst):
        """ Converts a nested list of various sizes to a 2D array block - extends scalars to the length of the rest of the array """
        max_length = max((np.size(v) if np.ndim(v)>0 else 1) for v in lst)
        lst = [np.full(max_length, v) if np.isscalar(v) else v for v in lst]
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
        
        #Flag variable uncertainty to be calculated
        var.uncertainty.is_calculated = True
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
    
    def calculateCorrelation(self, var):
        """ Calculates the correlation of the uncertainty for each timestep """
        #check for direct contributions, dependency contributions
        #Extract values, combine
        
    
    def splitToSourceContributions(self, var, recalculate=False):
        """ Splits uncertainty contributions down to all root uncertainty sources """
        if not var.uncertainty.is_calculated:
            raise ValueError(f"Cannot split uncertainty contributions for variable {var.name}, please calculate uncertainties first")
        
        #If the root split is already calculated and forced recalculation is not chosen, we simply return previously calculated values
        if var.uncertainty.root_uncertainty_contribution_split is not None and not recalculate:
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
            dep_sources, dep_root_split = self.splitToSourceContributions(var.dependencies[dep_name], recalculate=recalculate)
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
            
    



