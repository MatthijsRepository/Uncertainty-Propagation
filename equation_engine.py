from my_dataclasses import Variable
import re
import sympy as sp


""" This engine handles the initialization of the equation tree, 
    verifies its internal consistency and populates missing basic values """
class EquationEngine:
    def __init__(self, variables, update_dependencies=False):
        self.variables = variables          #dict: dictionary of variable names and Variable objects
        self.basic_variables, self.derived_variables = self.splitBasicDerived() #lists of variable names for basic and derived variables
        for name in self.derived_variables:
            self.checkVariableEquationConsistency(self.variables[name], update_dependencies)
        
    def splitBasicDerived(self, variables=None):
        """ Splits a tree into basic and derived variables """
        #If no variables provided, act on own registry
        if variables is None:
            variables = self.variables
        
        basic_variables = []
        derived_variables = []
        #Separate variables into basic and derived variables
        for name, var in variables.items():
            if var.is_basic:
                basic_variables.append(name)
            else:
                derived_variables.append(name)
        return basic_variables, derived_variables

    def equationReader(self, variable):
        """ This function extracts constituent variables from an equation """
        dependency_names = re.findall(r"'(.*?)'", variable.equation)                   
        return list(set(dependency_names))
    
    def checkVariableEquationConsistency(self, variable, update_dependencies=False):
        """ This function verifies whether specified variable dependency matches a given equation, optionally updates dependencies of variable to detected variables """
        dependency_names = self.equationReader(variable)
        if not set(dependency_names).issubset(variable.dependency_names):
            raise ValueError(f"Listed equation variables of variable {variable.name} do not match equation. Please check inputs. \nEquation: {variable.equation} \nDetected variables: {dependency_names} \nListed variables: {variable.dependency_names}")
        #update equation variables to detected variables
        if update_dependencies:
            variable.dependency_names = dependency_names
            
    def updateDependencyNames(self, variables=None):
        """ Updates listed dependency names for a variable set to only include those detected in the listed equation """
        if variables is None:
            variables = self.variables
            derived_variables = self.derived_variables
        else:
            derived_variables = self.splitBasicDerived(variables)[1]
            
        for name in derived_variables:
            var = variables[name]
            dependency_names = self.equationReader(var)
            var.dependency_names = dependency_names
            
    def _checkEquationTreeRecursive(self, variables, variables_to_check, silent):
        """ Internal function that recursively checks equation tree consistency """
        #Recursively navigates down the tree and checks if the equation tree is defined in a consistent way
        for name in variables_to_check:
            var = variables.get(name)
            #Check if variable in variables list
            if var is None:
                raise ValueError(f"Equation tree consistency check failed: variable {name} not included in the variable set. \nVariable set: {variables.keys()}.")
            
            #If variable already checked, pass this variable
            if var.is_basic:
                continue
            if var.is_root_consistent: 
                continue
            #Else: loop through all variables in dependencies
            for dep_name in var.dependency_names:
                dep = variables.get(dep_name)
                #If listed variable not included in variables list then the equation tree is not consistent
                if dep is None:
                    raise ValueError(f"Equation tree of variable {name} not consistent: variable {dep_name} not recognized.")
                #If variable is basic we can continue
                if dep.is_basic:
                    continue
                #If tree section is verified we can continue
                elif dep.is_root_consistent:
                    continue
                #If variable is not basic and not yet verified: check tree consistency
                elif self._checkEquationTreeRecursive(variables, dep.dependency_names, silent):
                    continue
                else:
                    raise ValueError(f"Equation of {name} is not root-consistent: check for {dep_name} failed")
            #If we succesfully looped through all dependencies of this variable, we can flag it as safe, and potentially log this to output
            if not silent:
                print(f"{name} is root-consistent")
            var.is_root_consistent = True
        #If no errors were encountered we can safely return True
        return True

    def checkEquationTreeConsistency(self, variables=None, variables_to_check=None, silent=True):
        """ This function that handles input and executes the equation tree consistency checker """
        #If no variables provided: act on own registry
        if variables is None:
            variables = self.variables
        #If no variables to check are provided: split provided variables
        if variables_to_check is None:
            if variables is self.variables:
                variables_to_check = self.derived_variables
            else:
                variables_to_check = self.splitBasicDerived(variables)
  
        #If variables_to_check is a single variable, we convert it to a list
        if isinstance(variables_to_check, str):
            variables_to_check = [variables_to_check]
        
        return self._checkEquationTreeRecursive(variables, variables_to_check, silent)
    
    def populateVariableDependencies(self, var, variables=None):
        """ Populates the dependencies for a single variable """
        if variables is None:
            variables = self.variables
            
        for dep_name in var.dependency_names:
            var.dependencies.append(variables[dep_name])
            
    def populateEquationTreeDependencies(self, variables=None, derived_variables=None):
        """ Populates dependencies for all dependent variables """
        if variables is None:
            variables = self.variables
            derived_variables = self.derived_variables
        if derived_variables is None:
            derived_variables = self.splitBasicDerived(variables)[1]
        
        for name in derived_variables:
            var = variables[name]
            self.populateVariableDependencies(var, variables)
    
    def buildEquation(self, var, symbol_map=None):
        """ Builds equation executable of a given variable, potentially using a provided sympy symbol map """
        #If Symbol map is not provided, build one from the dependency names of the variable
        if symbol_map is None:
            symbol_map = {name: sp.Symbol(name) for name in var.dependency_names}
        
        #Get symbol list for this variable
        symbols = [symbol_map[name] for name in var.dependency_names]
        #Remove quotes from equation
        clean_eq = var.equation
        for dep in var.dependency_names:
            clean_eq = clean_eq.replace(f"'{dep}'", dep)
        
        #Create sympy equation - this is later also used for differentiation
        var.sympy_equation = sp.sympify(clean_eq, locals=symbol_map)
        
        #Create callable function
        executable_equation = sp.lambdify(symbols, var.sympy_equation)
        def wrapper():
            args = [dep.values for dep in var.dependencies]
            return executable_equation(*args)
        var.callable = wrapper
            
    def buildSymPyEquationTree(self, variables=None):
        """ Goes through all derived variables in a variable set and builds their equations using SymPy """
        if variables is None:
            variables = self.variables
            derived_variables = self.derived_variables
        else:
            derived_variables = self.splitBasicDerived(variables)
        
        symbol_map = {name: sp.Symbol(name) for name in variables.keys()}
        for name in derived_variables:
            var = variables[name]
            self.buildEquation(var, symbol_map)

                

    
                
            
            
                            







