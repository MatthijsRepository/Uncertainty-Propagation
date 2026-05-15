from engines.my_dataclasses import Variable
import re
import sympy as sp


class EquationEngine:
    """
    This engine handles the building of a functional equation tree from a set of variables.
    
    The equation engine has functionality to read the regex of variable equations, 
    extract variable dependencies and match these with a variable registry. building the equation tree.
    The engine also checks whether the equation tree is well-defined and does not contain circular definitions.
    Using SymPy, the equation engine is capable of converting equations to python executables.
    Also using Sympy, the equation engine can take partial derivatives with repsect to dependencies and create their executables,
    which is required for the evaluation of uncertainty.
    
    Attributes
    ----------
    variables: dict[str, Variable]
        Variable registry, in which the engine looks for variables matching other variable's dependencies.
    """
    def __init__(self, variables):
        self.variables = variables          #dict: dictionary of variable names and Variable objects
        self.basic_variables, self.derived_variables = self.splitBasicDerived() #lists of variable names for basic and derived variables
        #Automatically populate the dependency names of the variables        
        self.populateVariableDependencyNames() ###!!!
        
    def splitBasicDerived(self, variables=None):
        """ 
        Function that partitions a dictionary of variables into lists of basic and derived variables.
        
        Parameters
        ----------
        variables: dict[str, Variable] or None, default=None
            Dictionary of variables to be split. If None, uses internal variable registry.
            
        Returns
        -------
        tuple[ list[Variable], list[Variable] ]
            Two lists of variables: the basic variables and the derived variables, respectively.
        """
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
        """ 
        Helper function that detects, extracts and subsequently removes top-level timesum expressions from equations.
        Returns list of top-level timesum expressions and a cleaned equation.
        
        Nested timesums are included in the returned top-level timesum expressions.
        
        Parameters
        ----------
        equation: str
            String following the custom equation regex used in this package.
        
        Returns
        -------
        tuple[ list[str], str]
            A list containing top-level timesum expressions
            A string containing the input equation with the top-level timesum expressions substituted.
        
        Notes
        -----        
        - In the overall workflow, a timesum expression prompts the creation of an auxiliary Variable with the 
          timesum contents as its equation and the ``is_timesum`` flag set to ``True``. 
        - In case of a nested timesum, the workflow detects the top-level timesum, creates an auxiliary variable for it
          and analyze this new variable's equation. Then it will encounter the nested timesum and repeat this workflow.
        """
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
        """ 
        Function that extracts top-level variable dependencies from a variable's equation.
        Top-level variables are not nested in timesums, or dependencies of dependencies.
        
        Parameters
        ----------
        variable: Variable
            Variable for which to read the equation.
        
        Returns
        -------
        list[str]
            List containing the variable names of the detected dependencies.
        """
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
        """ 
        Populates the dependency names for a set of variables, based on the variables detected in their equations.
        
        Parameters
        ----------
        variables: Variable or dict[str, Variable] or None
            Variable or dictionary of variables for which to populate the dependency names.
            If None, function acts on engine's internal variable registry.
        """
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
        
    def createTimeSumVariable(self, ts_str):
        """ 
        Creates a timesum variable from a timesum string segment.
        
        Parameters
        ----------
        ts_str
            String containing the timesum equation fragment.
            String of the form: TS_(equation, 'aggregation=' 'aggregate' or 'sum', 'rate=' 'true' or 'false')
            Only equation is required, if settings are not passed, aggregation rules are inferred from dependencies.
        
        Returns
        -------
        Variable
            Timesum variable with as equation the equation inside the timesum.
        """
        #extract equation; ts_str is of the form: TS_(equation, options) or TS_(equation)
        data = ts_str[4:-1].strip().split(",")
        if len(data)>1:
            equation = data[0]
            aggregation_rule, is_rate = self._getTimeSumSettingsFromString(data[1:])
        else:
            equation = data[0]
            aggregation_rule = is_rate = None
        var = Variable(name=ts_str, description=f"Timesum of: {equation}", aggregation_rule=aggregation_rule, is_rate=is_rate, \
                       is_basic=False, equation=equation, is_timesum=True) 
        self.populateVariableDependencyNames(var)
        return var
    
    def _getTimeSumSettingsFromString(self, settings):
        """ 
        Helper function to `createTimeSumVariable`, converts the settings in the timesum string to corresponding booleans.
        
        Parameters
        ----------
        settings: str
            String containing the settings, in format: 'aggregation=' 'average' or 'sum', 'rate=' 'true' or 'false'
        
        Returns
        -------
        tuple[str or None, bool or None]
            Detected aggregation rule, or `None` if not specified.
            Detected whether variable is a rate or not, or `None` if not specified.
        """        
        aggregation_rule = None
        is_rate = None
        for setting in settings:
            setting = setting.strip().lower()
            if setting.startswith("agg"):
                aggregation_rule = setting.split("=")[1].strip()
            elif setting.startswith("rate"):
                is_rate = setting.split("=")[1].strip()
                if is_rate == "true":
                    is_rate = True
                else:
                    is_rate = False
        return aggregation_rule, is_rate
    
    def _getAggregationRulesFromDependencies(self, var):
        """ 
        Function that infers aggregation rules from a variables dependencies.
        Recurses depth-first down the tree if the aggregation rules of a variable are not specificied.
        
        Note that the function does not perform a dimensional analysis, but infers aggregation rules using a simple decision rule.
        If any dependency is a rate over time and not a timesum, then the variable will also be a rate over time.
        If any dependency has "sum" as its aggregation rule, then the aggregation rule of this variable will also be "sum".
        If no dependency has "sum" as its aggregation rule, but any dependency has "average", then it will be "average".
        Defaults to `None` otherwise.
        
        Parameters
        ----------
        var: Variable
            Variable for which to retrieve the aggregation rules.
        
        Returns
        -------
        tuple[str or None, bool or None]
            Determined aggregation rule, or `None` if no rules specified downtree.
            Determined whether variable is a rate or not, or `None` if not specified downtree.
        
        Notes
        -----
        - It is recommended to manually specify aggregation rules to avoid any mistakes. The decision rules of this function do not
          reflect the actual dimensional analysis which actually determines the aggregation rule.
        - Populates timesum settings of dependencies, if these are not already defined.
        """
        rules = []
        is_rate = None
        for dep in var.dependencies.values():
            #Check whether the dependency has populated aggregation rules and is_rate booleans, otherwise it retrieves and populates them
            if (dep.aggregation_rule is None or dep.is_rate is None) and not dep.is_basic:
                #Retrieve rules from dependencies
                retr_agg_rule, retr_is_rate = self._getAggregationRulesFromDependencies(dep)
                #Populate fields in the dependency if they are empty
                if dep.aggregation_rule is None:
                    dep.aggregation_rule = retr_agg_rule
                if dep.is_rate is None:
                    dep.is_rate = retr_is_rate
                
            #Append found rules to the list
            rules.append(dep.aggregation_rule)
            #If the dependency is a rate and not a timesum, then the result must be a rate.
            if dep.is_rate and not dep.is_timesum:
                is_rate=True
            
        #Declare the fields in the variable under consideration
        # summing > averaging > None    ###!!! Note: ratio between two extensive quantities is an intensive quantity, 
        ### this function does not handle any complex algebra for this. When in doubt: specify manually
        if "sum" in rules:
            rule = "sum"
        elif "average" in rules:
            rule = "average"
        else:
            rule = None
        return rule, is_rate
        
    def populateEquationTreeTimeSumSettings(self, variables=None):
        """ 
        Ensures all timesum variables have their timesum settings defined. If not passed in the equation, infers from dependencies.
        
        Parameters
        ----------
        variables: dict[str, Variable] or None
            Dictionary of variables on which to act, will only act on timesum variables.
            If `None`, acts on internal variable registry.
        
        Raises
        ------
        ValueError
            If a timesum has no aggregation rules defined and none could be inferred from dependencies.
            Aggregation rule should be manually specified in the timesum equation.
        """
        if variables is None:
            variables = self.variables
        
        for var in variables.values():
            if var.is_timesum:
                rule, is_rate = self._getAggregationRulesFromDependencies(var)
                if var.aggregation_rule is None:
                    var.aggregation_rule = rule
                if var.is_rate is None:
                    var.is_rate = is_rate    
                if (var.aggregation_rule is None) or (var.is_rate is None):
                    raise ValueError(f"Could not infer aggregation rules for timesum {var.name}. Please specify manually in the equation tree.")
        
    def _checkEquationTreeRecursive(self, variables, variables_to_check, silent, indent="", stack=[]):
        """ 
        Recursive function that verifies well-definedness of the equation tree and creates timesum variables if encountered.
        Recursion handler helper function for `checkEquationTreeConsistency`.
        
        For a given set of variables to check, performs a depth-first verification that there are no circular definitions in the equation tree
        and that all variable dependencies exists. Creates a timesum variable if one is encountered, and adds it to variable registry.
        Recursion ends if all variables in `variables_to_check` are root consistent. Recursion backtracks if a variable is basic, or flagged as root-consistent.
        
        Parameters
        ----------
        variables: dict[str, Variable]
            Variable registry which we are checking with. Encountered variables are looked up from (and appended to inc ase of timesums) this registry.
        variables_to_check: list[str]
            Backlog of variables' names to check the downtree equation tree for.
        silent: bool
            Whether to print out the checking process. Can be used for debugging.
        indent: str
            Helper variable providing the print indent.
        stack: list[str]
            Stack of variable names detailing the propagation path from the node calling the recursion. Used to catch circular definitions.
        
        Returns
        -------
        bool: True
            Returns `True` if all variables in `variables_to_check` are root-consistent. Raises an error otherwise.
        
        Raises
        ------
        ValueError
            In case a circular definition is detected in the equation tree.
        KeyError
            If a dependency is encountered that is not included in the given variable set (and is not a timesum).
        """
        #Recursively navigates down the tree and checks if the equation tree is defined in a consistent way
        for name in variables_to_check:
            #Check for circular definitions
            if name in stack:
                raise ValueError(f"Equation tree consistency check failed: tree is circularly defined. Variable {name} appeared twice. Equation stack: {stack + [name]}.")
            
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
                    raise KeyError(f"Equation tree consistency check failed: variable {name} not included in the variable set. \nVariable set: {variables.keys()}.")
            
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
            self._checkEquationTreeRecursive(variables, var.dependency_names, silent=silent, indent=f"   {indent}", stack=stack+[name])

            #If we succesfully looped through all dependencies of this variable, we can flag it as safe, and potentially log this to output
            if not silent:
                print(f"{indent}{name} is root-consistent")
            var.is_root_consistent = True
        #If no errors were encountered we can safely return True
        return True

    def checkEquationTreeConsistency(self, variables=None, derived_variables=None, silent=True):
        """ 
        Checks equation tree consistency of the passed variable registry.
        Ensures input is of the correct format, calls the recursive equation tree check and potentially updates derived variable registry.
        
        Parameters
        ----------
        variables: dict[str, Variable] or None, default=None
            Variable registry used to check consistency against. All variables in the equation tree must be contained in this dictionary.
            If `None`, function will act on equation engine's internal variable registry.
        derived_variables: list[str], str or None, default=None
            List of names of all derived variables of the passed `variables` registry.
            Optional, if not passed the function will build this list automatically.
        silent: bool, default=True
            Boolean on whether to print out the recursive consistency check, for debugging purposes.
        
        Returns
        -------
        tuple[dict[str, Variable], list[str]]
            Updated variable dictionary and updated list of derived variable's names.
            Updated registry also contains encountered timesum variables that were built during the consistency check.
        """
        #If no variables provided: act on own registry
        if variables is None:
            variables = self.variables
        #If no variables to check are provided: split provided variables
        if (derived_variables is None) and (variables is self.variables):
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
        return variables, self.splitBasicDerived(variables)[1]
    
    def populateVariableDependencies(self, var, variables=None):
        """ 
        Populates the dependencies for a single variable by matching detected dependency names to variables in a registry.
        
        Parameters
        ----------
        var: Variable
            Variable for which the dependencies must be populated.
        variables: dict[str, Variable] or None, default=None
            Variable registry to retrieve dependencies from. If none is passed, uses engine's internal variable registry.
        """
        if variables is None:
            variables = self.variables
        
        if var.dependency_names is None:
            self.populateVariableDependencyNames(var)
        
        for dep_name in var.dependency_names:
            var.dependencies[dep_name] = variables[dep_name]
            
    def populateEquationTreeDependencies(self, variables=None, derived_variables=None):
        """ 
        Populates dependencies of all derived variables in the equation tree.
        
        Parameters
        ----------
        variables: dict[str, Variable] or None, default=None
            Variable registry to populate dependencies for and retrieve dependencies from. If None is passed, uses engine's internal variable registry.
        """
        if variables is None:
            variables = self.variables
            derived_variables = self.derived_variables
        if derived_variables is None:
            derived_variables = self.splitBasicDerived(variables)[1]
        
        for name in derived_variables:
            var = variables[name]
            self.populateVariableDependencies(var, variables)
        
    def _buildSymPySymbolMap(self, variables):
        """ 
        Builds a dictionary connecting variable names to their `SymPy.Symbol` object.
        Can act on either a variable dictionary, or on a single variable (in which case it returns a dictionary for its dependencies).
        
        Parameters
        ----------
        variables: dict[str, Variable] or Variable
            Dictionary of variables to build the mapping for. If a single variable, builds the map for this variable's dependencies.
        
        Returns
        dict[str, sp.Symbol]
            Mapping of variable names to their respective sympy symbols.
        """
        if isinstance(variables, dict):
            names = variables.keys()
        else:
            names = variables.dependency_names
        return {name: sp.Symbol(self._cleanEquationForSymPy(name)) for name in names}
            
    def _cleanEquationForSymPy(self, equation):
        """ 
        Cleans equation strings such that they are interpretable as equations by SymPy.
        
        Replaces top-level timesums by dummy syntax that does not interfere with sympy.
        Removes spaces and quotation marks from the string.
        Cleaned equation is stored in `variable.sympy_equation` class attribute.
        
        Parameters
        ----------
        equation: str
            Equation string following the format outlined in the example equationt tree.
        
        Returns
        -------
        str: Equation string cleaned for interpretation by SymPy.
        """
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
                    elif equation[i] in ["+", "-", "*","/","="]: #Math symbols are not allowed inside TS expressions
                        cleaned_equation += "_"
                        i += 1 ; continue
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
        """ 
        Builds the executable of a variable based on its equation string.
        Executable parameters are the variables `var` depends on, in the order they appear in the `variable.dependency_names` list.
        Executable is stored in `variable.executable` class attribute.
        
        Parameters
        ----------
        var: Variable
            Variable for which executable must be built.
        symbol_map: dict[str, sp.Symbol] or None, default=None
            Dictionary used to map dependency names to their sympy symbols.
            If None, function will build a symbol map itself.
        """
        #If Symbol map is not provided, build one from the dependency names of the variable
        if symbol_map is None:
            symbol_map = self._buildSymPySymbolMap(var)        
        
        #Access sympy symbols from map
        var.sympy_symbol_map = {name : symbol_map[name] for name in var.dependency_names}
        symbols = var.sympy_symbol_map.values()
        #Clean equation for sympy
        cleaned_equation = self._cleanEquationForSymPy(var.equation)
        #Create sympy equation - this is later also used for differentiation
        var.sympy_equation = sp.sympify(cleaned_equation, locals=symbol_map)
        
        #Create callable function
        #Note: in case the function is trivial (equation = 'X'), we don't want the equation to return the input variable 
        #Instead we want to return the input variable's values, the wrapper handles this
        #If the function is not trivial, we perform regular executable creation
        if isinstance(var.sympy_equation, sp.Symbol):
            def wrapper(temp_var):
                if isinstance(temp_var, Variable):
                    return temp_var.values
                else:
                    return temp_var #If we have numerical input, we simply return it back
            var.executable = wrapper
        else:
            var.executable = sp.lambdify(symbols, var.sympy_equation)
        
        
    def buildEquationTreeExecutables(self, variables=None):
        """ 
        Builds the equation executables for each dependent variable in an equation tree.
        
        Parameters
        ----------
        variables: dict[str, Variable] or None, default=None
            Variable registry for which to build equation executables.
            If None, function acts on the equation engine's internal variable registry.
        """
        if variables is None:
            variables = self.variables
            derived_variables = self.derived_variables
        else:
            derived_variables = self.splitBasicDerived(variables)
        
        symbol_map = self._buildSymPySymbolMap(variables)
        for name in derived_variables:
            var = variables[name]
            self.buildVariableExecutable(var, symbol_map)
        
    def buildPartialDerivativeExecutables(self, var, force_rebuild=False):
        """ 
        Builds the executable for the partial derivatives of the variable with respect to each dependency.
        
        For each dependency, the partial derivative of the variable's equation is taken by sympy.
        An executable is built for the resulting equation.
        Parameters of all executables are the variables `var` depends on, in the order they appear in the 
        `variable.dependency_names` list, regardless of whether they actually play a role in the equation.
        Executables are stored in a dictionary dict[str, func], with dependency names as keys, in `variable.partial_executables` attribute.
        
        Parameters
        ----------
        var: Variable
            Variable for which to build all partial derivative executables.
        force_rebuild: bool, default=False
            Boolean indicating whether the partial derivatives must be rebuilt if they are found to already exist.
        """
        #If partials are already built and forced rebuilding is not selected, simpyl return immediately
        if var.partial_executables is not None and force_rebuild is False:
            return
        #Else: rebuild partial executables
        var.partial_executables = {}
        
        for dep_name in var.dependency_names:
            partial_eq = sp.diff(var.sympy_equation, var.sympy_symbol_map[dep_name])
            #Create callable function
            executable = sp.lambdify(var.sympy_symbol_map.values(), partial_eq)           
            #Append executable to partial executables dictionary
            var.partial_executables[dep_name] = executable
            
            

                

