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
    
    
    def _convertListToArray(self, lst):
        """ Converts a nested list of various sizes to a 2D array block - extends scalars to the length of the rest of the array """
        max_length = max((np.size(v) if np.ndim(v)>0 else 1) for v in lst)
        
        lst = [np.full(max_length, v) if np.isscalar(v) else v for v in lst]
        return np.stack(lst, axis=0).astype(float)
      

        
    
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
            weighted_dep_uncertainties = self._convertListToArray(weighted_dep_uncertainties)
            
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
        var.uncertainty.direct_uncertainty_contributions = var.uncertainty.direct_uncertainties / np.sum(var.uncertainty.direct_uncertainties, axis=0)
        return
    
    def splitTotalUncertaintyContributions(self, var):
        """ Splits the fractional contributions to the total uncertainty between direct sources and all dependencies """
        if not var.uncertainty.is_calculated:
            raise ValueError(f"Cannot split uncertainty contributions for variable {var.name}, please calculate uncertainties first")
        #First we split to fractional contributions between direct sources and from dependencies
        total_sum = np.sum(var.uncertainty.weighted_dependency_uncertainties, axis=0) + var.uncertainty.total_direct_uncertainty
        
        var.uncertainty.total_direct_uncertainty_contribution = var.uncertainty.total_direct_uncertainty / total_sum
        var.uncertainty.dependency_uncertainty_contributions = var.uncertainty.weighted_dependency_uncertainties / total_sum
        return
    
    def splitToSourceContributions(self, var):
        """ Splits uncertainty contributions down to all root uncertainty sources """
        if not var.uncertainty.is_calculated:
            raise ValueError(f"Cannot split uncertainty contributions for variable {var.name}, please calculate uncertainties first")
        
    



