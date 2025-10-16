from my_dataclasses import Variable


def TIMESUM_TEMP(var, options=None):
    print(var)
    if isinstance(var, (float,int)):
        return var
    elif isinstance(var.values, (float, int)):
        raise ValueError(f"Timesums of constants are not supported yet (variable {var.name})")
    else:
        return var.TimeSum()


