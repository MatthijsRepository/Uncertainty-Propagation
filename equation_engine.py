from my_dataclasses import Variable
from sympy_timesum import timesum
import re
import sympy as sp

from calculation_engine import TIMESUM_TEMP ###!!!


""" This engine handles the initialization of the equation tree, 
    verifies its internal consistency and handles everything related to sympy """
class EquationEngine:
    def __init__(self, variables):
        self.variables = variables          #dict: dictionary of variable names and Variable objects
        self.basic_variables, self.derived_variables = self.splitBasicDerived() #lists of variable names for basic and derived variables
        
        self.populateVariableDependencyNames()
        
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
    
    
    def _equationTimeSumExtracter(self, equation):
        """ Detects, extracts and subsequently removes top-level timesum expressions from equations
        returns list of top-level timesum expressions and a cleaned equation """
        timesums = []
        i=0
        while i < len(equation):
            if equation[i:i+3] == "TS_":
                #Start loop, depth is 1
                depth = 1
                #We jump to the position of the opening parenthesis of the timesum, in the next loop we will skip this position
                i += 3
                start = i
                while depth>0:
                    i+=1
                    if equation[i]=="(":
                        depth +=1
                    if equation[i]==")":
                        depth -= 1
                timesums.append(f"TS_{equation[start:i+1]}")
            i+=1
        #Removing detected timesum statements from equation
        clean_eq = equation
        for ts in timesums:
            #We use re.escape in order to avoid regex errors due to presence of parentheses
            clean_eq = re.sub(re.escape(ts), " ", clean_eq)
        return timesums, clean_eq

    def equationReader(self, variable):
        """ This function extracts top-level constituent variables from an equation """
        dependency_names = []
        
        #first clean the equation such that 'timesum' is replaced by 'TS_'
        variable.equation = re.sub("timesum", "TS_", variable.equation)
        
        #Now we extract only top-level timesum statements, nested timesums are ignored
        timesums, clean_eq = self._equationTimeSumExtracter(variable.equation)
        dependency_names.extend(timesums)

        #Extend dependencies with regular variables left in the equation after top-level timesum statements are removed
        dependency_names.extend(re.findall(r"'(.*?)'", clean_eq))
        return list(set(dependency_names))
    
    
    def populateVariableDependencyNames(self, variables=None):
        """ Updates listed dependency names for a variable set to only include those detected in the listed equation """
        if variables is None:
            variables = self.variables
            derived_variables = self.derived_variables
        elif isinstance(variables, dict):
            derived_variables = self.splitBasicDerived(variables)[1]
        else:
            #In this case, variables is a single variable
            variables.dependency_names = self.equationReader(variables)
            return
        #Loop through the variable dictionary and update all dependencies
        for name in derived_variables:
            var = variables[name]
            var.dependency_names = self.equationReader(var)            
        
    def createTimeSumVariable(self, name):
        """ Creates a new variable that is a timesum of other variables """
        #extract equation: name is TS(equation; options)
        data = name[4:-1].strip().split(",")
        if len(data)==2:
            equation = data[0]
            settings = data[1]
        else:
            equation = data[0]
            settings = None
        var = Variable(name=name, description=f"Timesum of: {equation}, with settings {settings}", \
                       is_basic=False, equation=equation, is_timesum=True, timesum_settings=settings)
        self.populateVariableDependencyNames(var)
        return var
        
    def _checkEquationTreeRecursive(self, variables, variables_to_check, silent, indent=""):
        """ Internal function that recursively checks equation tree consistency """
        #Recursively navigates down the tree and checks if the equation tree is defined in a consistent way
        for name in variables_to_check:
            if not silent: 
                print(f"{indent}Tree checking: {name}")
            var = variables.get(name)
            #Check if variable in variables list
            if var is None:
                #If var is a timesum variable we create it here, otherwise we raise an error
                if name.startswith("TS_("):
                    var = self.createTimeSumVariable(name)
                    variables[name] = var
                else:
                    raise ValueError(f"Equation tree consistency check failed: variable {name} not included in the variable set. \nVariable set: {variables.keys()}.")
            
            #If variable already checked, pass this variable
            if var.is_basic:
                if not silent:
                    print(f"{indent}Variable {name} is basic")
                continue
            if var.is_root_consistent: 
                if not silent:
                    print(f"{indent}Variable {name} is consistent")
                continue
            #Else: check consistency of dependencies
            self._checkEquationTreeRecursive(variables, var.dependency_names, silent=silent, indent=f"   {indent}")

            #If we succesfully looped through all dependencies of this variable, we can flag it as safe, and potentially log this to output
            if not silent:
                print(f"{indent}{name} is root-consistent")
            var.is_root_consistent = True
        #If no errors were encountered we can safely return True
        return True

    def checkEquationTreeConsistency(self, variables=None, derived_variables=None, silent=True):
        """ This function that handles input and executes the equation tree consistency checker """
        #If no variables provided: act on own registry
        if variables is None:
            variables = self.variables
        #If no variables to check are provided: split provided variables
        if derived_variables is None:
            if variables is self.variables:
                derived_variables = self.derived_variables
            else:
                derived_variables = self.splitBasicDerived(variables)[1]
  
        #If variables_to_check is a single variable, we convert it to a list
        if isinstance(derived_variables, str):
            derived_variables = [derived_variables]
        
        #Check equation tree
        self._checkEquationTreeRecursive(variables, derived_variables, silent)
        
        #Update derived variables set to include newly created timesum variables
        if variables is self.variables:
            self.derived_variables = self.splitBasicDerived(variables)[1]
        else:
            return self.splitBasicDerived(variables)[1]
    
    def populateVariableDependencies(self, var, variables=None):
        """ Populates the dependencies for a single variable using the variables in the variables registry """
        if variables is None:
            variables = self.variables
        
        if var.dependency_names is None:
            self.populateVariableDependencyNames(var)
        
        for dep_name in var.dependency_names:
            var.dependencies[dep_name] = variables[dep_name]
            
    def populateEquationTreeDependencies(self, variables=None, derived_variables=None):
        """ Populates dependencies for all dependent variables in the tree """
        if variables is None:
            variables = self.variables
            derived_variables = self.derived_variables
        if derived_variables is None:
            derived_variables = self.splitBasicDerived(variables)[1]
        
        for name in derived_variables:
            var = variables[name]
            self.populateVariableDependencies(var, variables)
    
    
    def _buildSymbolMap_OLD(self, variables):
        """ Helper function that builds a symbol map for a variable dictionary or a single variable dependency list
            we add an alias for non-timesum variables as they are in quotes in the equation
            so but G and 'G' will both refer to a symbol G """
        if isinstance(variables, dict):
            names = variables.keys()
        else:
            names = variables.dependency_names
        symbol_map = {}
        #We need to register quotes around variables, so we add them retroactively to our symbol map for non-timesum variables
        #So: TS_('G') is unchanged, but G will become 'G', which is needed for equation readout
        for name in names:
            symbol_map[name] = sp.Symbol(name)
            # Add quoted alias for non-timesum names
            if not name.startswith("TS_("):
                symbol_map[f"'{name}'"] = symbol_map[name]
                
        #symbol_map = {name: sp.Symbol(name) for name in names}
        return symbol_map
    
    def _cleanEquationForSymPy_OLD(self, equation):
        """ Cleans quotes around variables that are outside of timesum operators, for use in sympy """
        i=0
        cleaned_equation = ""
        while i < len(equation):
            #If we register a timesum: pass this section
            if equation[i:i+3] == "TS_":
                #Skip this part, add it to our equation
                cleaned_equation += "TS_"
                i += 3
                #Start loop, depth is 1
                depth = 1

                while depth>0:
                    cleaned_equation += equation[i]
                    i+=1
                    if equation[i]=="(":
                        depth +=1
                    if equation[i]==")":
                        depth -= 1
            #If we did not encounter a timesum: check if there is a quote and remove it if necessary
            if not equation[i] == "'":
                cleaned_equation += equation[i]
            i+=1
        return cleaned_equation
    
    def _buildSymPySymbolMap(self, variables):
        if isinstance(variables, dict):
            names = variables.keys()
        else:
            names = variables.dependency_names
        
        #cleaned_names = [re.sub("'", "", name) for name in names]               #Remove quotes
        #cleaned_names = [re.sub(",", "", name) for name in cleaned_names]       #Remove commas
        #return {name : sp.Symbol(cleaned_names[i]) for i, name in enumerate(names)}
        return {name: sp.Symbol(self._cleanEquationForSymPy(name)) for name in names}
            
    
    def _cleanEquationForSymPy(self, equation):
        cleaned_equation = ""
        #Parse through equation and replace all parentheses related to timesums by double underscores, but not mathematical ones
        i = 0
        while i < len(equation):
            if equation[i:i+3] =="TS_":
                cleaned_equation += "TS__"
                depth=1
                i+=4
                while depth>0:
                    if equation[i]=="(":
                        depth +=1
                        cleaned_equation += "__"
                    elif equation[i]==")":
                        depth -= 1
                        cleaned_equation += "__"
                    else:
                        cleaned_equation += equation[i]
                    i += 1
            else:
                cleaned_equation += equation[i]
                i+=1
        #Remove all quotes
        cleaned_equation = re.sub("'","", cleaned_equation)
        #Remove all commas
        cleaned_equation = re.sub(",", "", cleaned_equation)
        #Remove all spaces
        cleaned_equation = re.sub(" ", "", cleaned_equation)
        return cleaned_equation
        
    
    
    def buildVariableExecutable(self, var, symbol_map=None):
        """ Builds equation executable of a given variable, potentially using a provided sympy symbol map """
        #If Symbol map is not provided, build one from the dependency names of the variable
        if symbol_map is None:
            symbol_map = self._buildSymPySymbolMap(var)        
        
        
        #Access sympy symbols from map
        symbols = [symbol_map[name] for name in var.dependency_names]
        #Clean equation for sympy
        cleaned_equation = self._cleanEquationForSymPy(var.equation)
        #Create sympy equation - this is later also used for differentiation
        var.sympy_equation = sp.sympify(cleaned_equation, locals=symbol_map)
        #Create callable function
        var.executable = sp.lambdify(symbols, var.sympy_equation, modules=[{"timesum": TIMESUM_TEMP}]) ###!!!
        """ 
        def wrapper():
            #args = [dep.values for dep in var.dependencies]
            args = [var.dependencies[dep_name] for dep_name in var.dependency_names]

            var.values = executable_equation(*args)
            return None
        var.executeEquation = wrapper
        """ 
            
    def buildEquationTreeExecutables(self, variables=None):
        """ Goes through all derived variables in a variable set and builds their equations using SymPy """
        if variables is None:
            variables = self.variables
            derived_variables = self.derived_variables
        else:
            derived_variables = self.splitBasicDerived(variables)
        
        symbol_map = self._buildSymPySymbolMap(variables)
        for name in derived_variables:
            var = variables[name]
            self.buildVariableExecutable(var, symbol_map)
            

                

    
                
            
            
                            







