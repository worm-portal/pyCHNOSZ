import os
import warnings
from contextlib import contextmanager
from IPython.display import Image, display

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.lib import grdevices
    
    CHNOSZ = importr("CHNOSZ")
    grdev = importr('grDevices')


@contextmanager
def r_inline_plot(width=600, height=600, dpi=100):

    with grdevices.render_to_bytesio(grdevices.png, 
                                     width=width,
                                     height=height, 
                                     res=dpi) as b:

        yield

    data = b.getvalue()
    display(Image(data=data, format='png', embed=True))

    
def convert_to_RVector(value, convert_lists=True):
    
    """
    Convert a value into a list (if it is not already) if convert_lists=True
    and returns an R vector of the appropriate type (bool, int, float, str).
    """
    
    if not isinstance(value, list) and not convert_lists:
        return value
    elif not isinstance(value, list) and convert_lists:
        value = [value]
    else:
        pass
    
    if all(isinstance(x, bool) for x in value):
        return ro.BoolVector(value)
    elif all(isinstance(x, int) for x in value):
        return ro.IntVector(value)
    elif all(isinstance(x, float) or isinstance(x, int) for x in value):
        return ro.FloatVector(value)
    else:
        return ro.StrVector(value)


def equilibrate(aout, balance=None, loga_balance=None, ispecies=None,
                normalize=False, messages=True):
    
    """
    Python wrapper for the equilibrate() function in CHNOSZ.
    """
        
    # stay.normal is an 'unused argument' in equilibrate()?
    
    args = {'aout':aout, 'normalize':normalize}
    
    if balance != None: args['balance'] = balance
    if loga_balance != None: args['loga.balance'] = loga_balance
    if ispecies != None: args['ispecies'] = convert_to_RVector(ispecies)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        a = CHNOSZ.equilibrate(**args)

    if messages:
        for warning in w:
            print(warning.message)
    
    return a
    

def diagram(eout, ptype='auto', alpha=False, normalize=False,
            as_residue=False, balance=None, groups=None,
            xrange=None, mar=None, yline=None, side=[1,2,3,4],
            ylog=True, xlim=None, ylim=None, xlab=None, ylab=None,
            cex=None, cex_names=None, cex_axis=None,
            lty=None, lwd=None, dotted=None,
            spline_method=None, contour_method=None, levels=None,
            col=None, col_names=None, fill=None,
            fill_NA="gray80", limit_water=None,
            names=None, format_names=True, bold=False, italic=False, 
            font=None, family=None, adj=0.5,
            dx=0, dy=0, srt=0, min_area=0,
            main=None, legend_x=None,
            add=False, plot_it=True, tplot=True,
            messages=True):
    
    """
    Python wrapper for the diagram() function in CHNOSZ.
    """
    
    args = {'eout':eout, 'ptype':ptype, 'alpha':alpha, 'normalize':normalize,
            'as.residue':as_residue, 'ylog':ylog, 'fill.NA':fill_NA, 'format.names':format_names,
            'bold':bold, 'italic':italic, 'adj':adj, 'dx':dx, 'dy':dy, 'srt':srt,
            'min.area':min_area, 'plot.it':plot_it, 'tplot':tplot}
    
    if balance != None: args["balance"] = balance
    if groups != None: args["groups"] = groups
    if xrange != None: args["xrange"] = convert_to_RVector(xrange)
    if mar != None: args["mar"] = convert_to_RVector(mar)
    if yline != None: args["yline"] = convert_to_RVector(yline)
    if side != None: args["side"] = convert_to_RVector(side)
    if xlim != None: args["xlim"] = convert_to_RVector(xlim)
    if ylim != None: args["ylim"] = convert_to_RVector(ylim)
    if xlab != None: args["xlab"] = xlab
    if ylab != None: args["ylab"] = ylab
    if cex != None: args["cex"] = convert_to_RVector(cex)
    if cex_names != None: args["cex.names"] = convert_to_RVector(cex_names)
    if cex_axis != None: args["cex.axis"] = convert_to_RVector(cex_axis)
    if lty != None: args["lty"] = convert_to_RVector(lty)
    if lwd != None: args["lwd"] = convert_to_RVector(lwd)
    if dotted != None: args["dotted"] = convert_to_RVector(dotted)
    if spline_method != None: args["spline.method"] = spline_method
    if contour_method != None: args["contour.method"] = contour_method
    if levels != None: args["levels"] = convert_to_RVector(levels)
    if col != None: args["col"] = convert_to_RVector(col)
    if col_names != None: args["col.names"] = convert_to_RVector(col_names)
    if fill != None: args["fill"] = convert_to_RVector(fill)
    if limit_water != None: args["limit.water"] = limit_water
    if names != None: args["names"] = convert_to_RVector(names)
    if font != None: args["font"] = convert_to_RVector(font)
    if family != None: args["family"] = convert_to_RVector(family)
    if main != None: args["main"] = main
    if legend_x != None: args["legend.x"] = legend_x
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        with r_inline_plot(width=1024, height=896, dpi=150):
            if isinstance(add, bool):
                if add: # add='True' does not work with the current pyCHNOSZ framework
                    raise Exception("The argument 'add' must be assigned the output of the previous diagram(s).")
                else:
                    a = CHNOSZ.diagram(**args)
            elif isinstance(add, list):
                args.update({'add':True})
                for to_add in add:
                    CHNOSZ.diagram(**to_add)
                a = CHNOSZ.diagram(**args)
                    
            else:
                CHNOSZ.diagram(**add)
                args.update({'add':True})
                a = CHNOSZ.diagram(**args)
            
    if messages:
        for warning in w:
            print(warning.message)
    
    return a, args


def affinity(property=None, sout=None, exceed_Ttr=False,
             exceed_rhomin=False, return_buffer=False, return_sout=False,
             balance="PBB", iprotein=None, loga_protein=-3, transect=None,
             messages=True, **kwargs):

    """
    Python wrapper for the affinity() function in CHNOSZ.
    """
    
    args = {'exceed_Ttr':exceed_Ttr, 'exceed_rhomin':exceed_rhomin,
            'return_buffer':return_buffer, 'return_sout':return_sout,
            'balance':balance, 'loga_protein':loga_protein}
    
    if property != None: args["property"] = property
    if sout != None: args["sout"] = sout
    if iprotein != None: args["iprotein"] = iprotein
    if transect != None: args["transect"] = transect
    
    for key, value in kwargs.items():
        if isinstance(value, list):
            value = ro.FloatVector(value)
        args.update({key:value})
        
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        a = CHNOSZ.affinity(**args)

    if messages:
        for warning in w:
            print(warning.message)

    return a

            
def species(species=None, state=None, delete=False, add=False, index_return=False, messages=True):
    
    """
    Python wrapper for the species() function in CHNOSZ.
    """
    
    args={}
    
    if species != None:
        args["species"] = convert_to_RVector(species, convert_lists=False)
            
    if state != None:
        args["state"] = convert_to_RVector(state, convert_lists=False)
        
    args["add"] = add
    args["delete"] = delete
    args["index.return"] = index_return
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        a = CHNOSZ.species(**args)
        
    if messages:
        for warning in w:
            print(warning.message)
    
    return pandas2ri.ri2py_dataframe(a)


def basis(species=None, state=None, logact=None, delete=False, messages=True):
    
    """
    Python wrapper for the basis() function in CHNOSZ.
    """
    
    args={}
    
    if species != None:
        args["species"] = convert_to_RVector(species, convert_lists=False)
            
    if state != None:
        args["state"] = convert_to_RVector(state, convert_lists=False)
    
    if logact != None: args["logact"] = convert_to_RVector(logact)
    
    args["delete"] = delete
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        a = CHNOSZ.basis(**args)
        
    if messages:
        for warning in w:
            print(warning.message)
    
    return pandas2ri.ri2py_dataframe(a)


def reset(messages=True):
    
    """
    Python wrapper for the reset() function in CHNOSZ.
    """
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        a = CHNOSZ.reset()
        
    if messages:
        for warning in w:
            print(warning.message)


def add_OBIGT(file, species=None, force=True, messages=True):
    
    """
    Python wrapper for the add.OBIGT() function in CHNOSZ.
    """
    
    args={'file':file}
    
    if species != None:
        if not isinstance(species, list):
            args["species"] = species
        else:
            args["species"] = convert_to_RVector(species)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        a = CHNOSZ.add_OBIGT(**args)
    
    if messages:
        for warning in w:
            print(warning.message)
    
    return list(a)
        
    
def info(species, state=None, check_it=True, messages=True):
    
    """
    Python wrapper for the info() function in CHNOSZ.
    """
    
    args = {}
    output_is_df = False
    
    args["species"] = convert_to_RVector(species, convert_lists=False)
    if not isinstance(species, list): species = [species]
    if all(isinstance(x, int) for x in species):
        output_is_df = True
    
    if state != None:
        args["state"] = convert_to_RVector(state, convert_lists=False)
    
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
    
    args = {'species':convert_to_RVector(species, convert_lists=False)}
        
    if coeff != None:
        args["coeff"] = convert_to_RVector(coeff)
    else:
        single_species = True
        
    if state != None:
        args["state"] = convert_to_RVector(state, convert_lists=False)
    
    args["property"] = convert_to_RVector(property)
    
    if T != None:
        args['T'] = convert_to_RVector(T, convert_lists=False)

    if P != None:
        args['P'] = convert_to_RVector(P, convert_lists=False)
    
    if grid != None: args['grid'] = grid # grid is either 'T' or 'P'
    
    args['convert'] = convert
    args['exceed.Ttr'] = exceed_Ttr
    args['exceed.rhomin'] = exceed_rhomin
    
    if logact != None: args["logact"] = convert_to_RVector(logact)
    
    args['autobalance'] = autobalance
    
    if IS != None:
        args["IS"] = convert_to_RVector(IS, convert_lists=False)
    
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