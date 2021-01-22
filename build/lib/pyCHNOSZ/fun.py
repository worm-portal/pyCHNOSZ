import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    
    CHNOSZ = importr("CHNOSZ")

    
def info(species, state=None, check_it=True, messages=True):
    
    """
    Python wrapper for the info() function in CHNOSZ.
    """
    
    args = {}
    output_is_df = False
    
    if not isinstance(species, list):
        args["species"] = species
        if isinstance(species, int):
            output_is_df = True
    else:
        if isinstance(species[0], int):
            output_is_df = True
            args["species"] = ro.IntVector(species)
        else:
            args["species"] = ro.StrVector(species)
    
    if state != None:
        if not isinstance(state, list):
            args["state"] = state
        else:
            args["state"] = ro.StrVector(state)
    
    args["check.it"] = check_it

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        a = CHNOSZ.info(**args)
    
    if messages:
        for warning in w:
            print(warning.message)
    
    if output_is_df:
        return pandas2ri.ri2py_dataframe(a)
    else:
        return list(a)

    
def subcrt(species, coeff=None, state=None,
           property=["logK","G","H","S","V","Cp"],
           T=None, P=None, grid=None,
           convert=True, exceed_Ttr=False, exceed_rhomin=False,
           logact=None, autobalance=True, IS=None, messages=True):
    
    """
    Python wrapper for the subcrt() function in CHNOSZ.
    """
    
    single_species = False
    
    if not isinstance(species, list):
        args = {'species':species}
    else:
        if isinstance(species[0], int):
            args = {'species':ro.IntVector(species)}
        else:
            args = {'species':ro.StrVector(species)}
        
    if coeff != None:
        args["coeff"] = ro.FloatVector(coeff)
    else:
        single_species = True
        
    if state != None:
        if not isinstance(state, list):
            args["state"] = state
        else:
            args["state"] = ro.StrVector(state)
    
    if not isinstance(property, list): property = [property]
    args["property"] = ro.StrVector(property)
    
    if T != None:
        if not isinstance(T, list):
            args['T'] = T
        else:
            args['T'] = ro.FloatVector(T)

    if P != None:
        if not isinstance(P, list):
            args['P'] = P
        else:
            args['P'] = ro.FloatVector(P)
    
    if grid != None: args['grid'] = grid # grid is either 'T' or 'P'
    
    args['convert'] = convert
    args['exceed.Ttr'] = exceed_Ttr
    args['exceed.rhomin'] = exceed_rhomin
    
    if logact != None: args["logact"] = ro.FloatVector(logact)
    
    args['autobalance'] = autobalance
    
    if IS != None:
        if not isinstance(IS, list):
            args["IS"] = IS
        else:
            args["IS"] = ro.FloatVector(IS)
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        a = CHNOSZ.subcrt(**args)
    
    if messages:
        for warning in w:
            print(warning.message)
    
    if len(a) == 3:
        warn = a[2][0] # subcrt's list includes warnings only if they appear
    else:
        warn = None
    
    if not single_species:
        out_dict = {"reaction":pandas2ri.ri2py_dataframe(a[0]),
                    "out":pandas2ri.ri2py_dataframe(a[1])} # the extra [0] is important
    else:
        out_dict = {"species":pandas2ri.ri2py_dataframe(a[0]), "out":{}}
        
        i=0
        for df in a[1]:
            out_dict["out"][out_dict["species"].name[i]] = pandas2ri.ri2py_dataframe(df)
            i += 1
        
    if warn != None:
        out_dict["warnings"] = warn
    
    return out_dict