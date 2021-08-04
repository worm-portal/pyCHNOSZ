import os
import warnings
from contextlib import contextmanager
from IPython.display import Image, display
import pandas as pd

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.lib import grdevices
    
    pandas2ri.activate()
    
    CHNOSZ = importr("CHNOSZ")
    grdev = importr('grDevices')

NumberTypes = (int, float, complex)

@contextmanager
def __r_inline_plot(width=600, height=520, dpi=150):
    
    """
    Display R plots inline.

    Parameters
    ----------
    width, height : numeric, default 600 by 520
        Width and height of the plot.
        
    dpi : numeric, default 150
        Resolution of the plot.
    """
    
    with grdevices.render_to_bytesio(grdevices.png, 
                                     width=width,
                                     height=height, 
                                     res=dpi) as b:

        yield

    data = b.getvalue()
    display(Image(data=data, format='png', embed=True))

    
def _convert_to_RVector(value, force_Rvec=True):
    
    """
    Convert a value or list into an R vector of the appropriate type.
    
    Parameters
    ----------
    value : numeric or str, or list of numeric or str
        Value to be converted.
    
    force_Rvec : bool, default True
        If `value` is not a list, force conversion into a R vector?
        False will return an int, float, or str if value is non-list.
        True will always return an R vector.
    
    Returns
    -------
    int, float, str, an rpy2 R vector
        A value or R vector of an appropriate data type.
    """
    
    if not isinstance(value, list) and not force_Rvec:
        return value
    elif not isinstance(value, list) and force_Rvec:
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
            

def water(property=None, T=298.15, P="Psat", P1=True, messages=True):
    
    """
    Python wrapper for the water() function in CHNOSZ.
    Calculate thermodynamic and electrostatic properties of water.
    
    Parameters
    ----------
    property : str or list of str, optional
        Computational setting or property(s) to calculate. To set a water model,
        use 'SUPCRT92' (default) or 'SUPCRT', 'IAPWS95' or 'IAPWS', or 'DEW'.
        To calculate a property, use 'A', 'G', 'S', 'U', etc. See 
        http://chnosz.net/manual/water.html for the complete catalog of
        properties that can be calculated, their units, and their availability
        in water models.

    T : numeric, default 298.15
        Temperature (K)

    P : numeric, default "Psat"
        Pressure (bar), or Psat for vapor-liquid saturation.

    P1 : bool, default True
        Output pressure of 1 bar below 100 °C instead of calculated values of Psat?

    messages : bool, default True
        Display messages from CHNOSZ?

    Returns
    ----------
    out : float or dict
        Calculated value of desired water property. If `property` is a list,
        returns a dictionary of calculated values.
    """
    
    if property == None:
        property = ro.r("NULL")
    elif isinstance(property, list):
        property = _convert_to_RVector(property, force_Rvec=True)
    else:
        pass
    
    args = {'property':property, 'T':T, 'P':P, 'P1':P1}
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = CHNOSZ.water(**args)

    if messages:
        for warning in w:
            print(warning.message)
    
    if property not in ["SUPCRT92", "SUPCRT", "IAPWS95", "IAPWS", "DEW"]:
        out = list(out)
        if len(out) == 1:
            if property == ro.r("NULL"):
                out = out[0]
            else:
                out = out[0][0]
        else:
            # CHNOSZ produces a single-row dataframe when multiple properties
            # are supplied to `property`. Here, a dictionary is returned.
            out = {prop:o[0] for prop,o in zip(property, out)}
        return out
        
        
def entropy(formula, messages=True):
    
    """
    Python wrapper for the entropy() function in CHNOSZ.
    Calculate the standard molal entropy of elements in a compound.
    
    Parameters
    ----------
    formula : str, int, or list of str or int
        Supply a chemical formula (e.g. "CH4"), a species index, or a list of
        formulas or indices.

    messages : bool, default True
        Display messages from CHNOSZ?
    
    Returns
    ----------
    out : float or list of float
        Standard molal entropy of elements in the formula in cal/(mol K).
        Returns a list if `formula` is a list.
    """
    
    formula_R = _convert_to_RVector(formula, force_Rvec=False)
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = CHNOSZ.entropy(formula_R)

    if messages:
        for warning in w:
            print(warning.message)
    
    out = list(out)
    if not isinstance(formula, list):
        out = out[0]
    
    return out


def mass(formula, messages=True):
    
    """
    Python wrapper for the mass() function in CHNOSZ.
    Calculate the mass of the sum of elements in a formula.
    
    Parameters
    ----------
    formula : str, int, or list of str or int
        Supply a chemical formula (e.g. "CH4"), a species index, or a list of
        formulas or indices.

    messages : bool, default True
        Display messages from CHNOSZ?
    
    Returns
    ----------
    out : float or list of float
        Molar mass of the sum of elements in the formula, in g/mol.
        Returns a list if `formula` is a list.
    """
    
    formula_R = _convert_to_RVector(formula, force_Rvec=False)
    print(formula_R)
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = CHNOSZ.mass(formula_R)

    if messages:
        for warning in w:
            print(warning.message)
    
    out = list(out)
    if not isinstance(formula, list):
        out = out[0]
    
    return out


def zc(formula, messages=True):
    
    """
    Python wrapper for the ZC() function in CHNOSZ.
    Calculate the average oxidation state of carbon (ZC) in a molecule.
    
    Parameters
    ----------
    formula : str or list of str
        Chemical formula(s) of molecules.

    messages : bool, default True
        Display messages from CHNOSZ?
    
    Returns
    ----------
    out : float or list of float
        The average oxidation state of carbon in the formula.
        Returns a list if `formula` is a list.
    """
    
    formula_R = _convert_to_RVector(formula, force_Rvec=False)
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = CHNOSZ.ZC(formula_R)

    if messages:
        for warning in w:
            print(warning.message)
    
    out = list(out)
    if not isinstance(formula, list):
        out = out[0]
    
    return out


def makeup(formula, multiplier=1, sum=False, count_zero=False, messages=True):
    
    """
    Python wrapper for the makeup() function in CHNOSZ.
    Parse a formula into a dictionary of elements and their counts.
    
    Parameters
    ----------
    formula : str or list of str
        Chemical formula or a list of chemical formulas.

    multiplier : numeric, default 1
        Multiplier for the elemental counts in each formula.
    
    sum : bool, default False
        Add together the elemental counts in all formulas?
        Ignored if `formula` is not a list.
    
    count_zero : bool, default False
        Include zero counts for an element if it is not in this formula, but is
        found in another formula in the list? Ensures that the dictionaries
        returned for each formula have the same length.
        Ignored if `formula` is not a list.

    messages : bool, default True
        Display messages from CHNOSZ?
    
    Returns
    ----------
    out : dict
        Dictionary of elements and their counts in the formula(s).
    """
    
    formula_R = _convert_to_RVector(formula, force_Rvec=False)
    
    args = {'formula':formula_R, "multiplier":multiplier,
            "sum":sum, "count.zero":count_zero}
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = CHNOSZ.makeup(**args)
        
    if messages:
        for warning in w:
            print(warning.message)
    
    if isinstance(out, ro.ListVector):
        out_dict = {}
        for i, molecule in enumerate(formula):
            molecule_dict = {}
            if not count_zero:
                out_names = out[i].names
            else:
                out_names = out[i].names[0]

            for ii, elem in enumerate(out_names):
                molecule_dict[elem] = out[i][ii]
            out_dict[molecule] = molecule_dict
            
    else:
        out_dict = {}
        if not sum or not isinstance(formula, list):
            out_names = out.names
        else:
            out_names = out.names[0]
        for i, elem in enumerate(out_names):
            out_dict[elem] = out[i]
    
    return out_dict


def seq2aa(protein, sequence, messages=True):
    
    """
    Python wrapper for the seq2aa() function in CHNOSZ.
    Returns a data frame of amino acid composition corresponding to the provided
    sequence.
    
    Parameters
    ----------
    protein : str
        Protein name with an underscore, e.g., 'LYSC_CHICK'
    
    sequence : str
        Amino acid sequence of the protein.

    messages : bool, default True
        Display messages from CHNOSZ?
    
    Returns
    -------
    Pandas dataframe
        Amino acid composition of protein.
    """
    
    args = {'protein':protein, 'sequence':sequence}

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        pout = CHNOSZ.seq2aa(**args)

    if messages:
        for warning in w:
            print(warning.message)
    
    return pandas2ri.ri2py_dataframe(pout)


def add_protein(aa, messages=True):
    
    """
    Python wrapper for the add.protein() function in CHNOSZ.
    Add proteins to the OBIGT thermodynamic database.
    
    Parameters
    ----------
    aa : Pandas dataframe
        Amino acid composition of protein(s) from `seq2aa`.
        
    messages : bool, default True
        Display messages from CHNOSZ?
    
    Returns
    -------
    list of int
        List of protein indices, iprotein.
    """
    aa = pandas2ri.py2ri(aa)
    args = {'aa':aa}

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        apout = CHNOSZ.add_protein(**args)

    if messages:
        for warning in w:
            print(warning.message)
    
    return list(apout)


def equilibrate(aout, balance=None, loga_balance=None, ispecies=None,
                normalize=False, messages=True):
    
    """
    Python wrapper for the equilibrate() function in CHNOSZ.
    Calculate equilibrium chemical activities of species from the affinities
    of formation of the species at unit activity.
    
    Parameters
    ----------
    aout : rpy2.ListVector
        Output from `affinity`.
    
    balance : str or numeric, optional
        How to balance the transformations.
    
    loga_balance : numeric or list of numeric, optional
        Logarithm of total activity of balanced quantity.
    
    ispecies : numeric, optional
        Which species to include.
    
    normalize : bool, default False
        Normalize the molar formulas of species by the balancing coefficients?
    
    messages : bool, default True
        Display messages from CHNOSZ?
    
    Returns
    -------
    eout : rpy2.ListVector
        Output from `equilibrate`.
    """
        
    # stay.normal is an unused argument in equilibrate()?
    
    args = {'aout':aout, 'normalize':normalize}
    
    if balance != None: args['balance'] = balance
    if loga_balance != None: args['loga.balance'] = loga_balance
    if ispecies != None: args['ispecies'] = _convert_to_RVector(ispecies)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        eout = CHNOSZ.equilibrate(**args)

    if messages:
        for warning in w:
            print(warning.message)
    
    return eout


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
            width=600, height=520, dpi=150,
            messages=True):
    
    """
    Python wrapper for the diagram() function in CHNOSZ.
    Plot equilibrium chemical activity (1-D speciation) or equal-activity
    (2-D predominance) diagrams as a function of chemical activities of
    basis species, temperature and/or pressure.
    
    Parameters
    ----------
    eout : rpy2.ListVector
        Output from `equilibrate` or `affinity`.
    
    ptype : str, default 'auto'
        Type of plot, or name of basis species whose activity to plot.

    alpha : bool or str (balance), default False
        For speciation diagrams, plot degree of formation instead of
        activities?

    normalize : bool, default False
        Divide chemical affinities by balance coefficients (rescale to whole
        formulas)?

    as_residue : bool, default False
        Divide chemical affinities by balance coefficients (no rescaling)?

    balance : str or numeric, optional
        How to balance the transformations.

    groups : rpy2.ListVector of numeric, optional
        Groups of species to consider as a single effective species.

    xrange : numeric, optional
        Range of x-values between which predominance field boundaries are
        plotted.

    mar : numeric, optional
        Margins of plot frame.

    yline : numeric, optional
        Margin line on which to plot the y-axis name.

    side : numeric, optional
        Which sides of plot to draw axes.

    xlim : numeric, optional
        Limits of x-axis.

    ylim : numeric, optional
        Limits of y-axis.

    xlab : str, optional
        Label to use for x-axis.

    ylab : str, optional
        Label to use for y-axis.

    ylog : bool, optional
        Use a logarithmic y-axis (on 1D degree diagrams)?

    cex : numeric, optional
        Character expansion (scaling relative to current).

    cex_names : numeric, optional
        Character expansion factor to be used for names of species on plots.

    cex_axis : numeric, optional
        Character expansion factor for names of axes.

    lty : numeric, optional
        Line types to be used in plots.

    lwd : numeric, optional
        Line width.

    dotted : numeric, optional
        How often to skip plotting points on predominance field boundaries (to
        gain the effect of dotted or dashed boundary lines).

    spline_method : str, optional
        Method used in splinefun.

    contour_method : str, optional
        Labelling method used in contour (use None for no labels).

    levels : numeric, optional
        Levels at which to draw contour lines.

    col : str, optional
        Color of activity lines (1D diagram) or predominance field boundaries
        (2D diagram).

    col_names : str, optional
        Colors for labels of species.

    fill : str, optional
        Colors used to fill predominance fields.

    fill_NA : str, optional
        Color for grid points with NA values.

    limit_water : None or bool, optional
        Set NaN values beyond water stability limits?

    names : str, optional
        Names of species for activity lines or predominance fields.

    format_names : bool, default True
        Apply formatting to chemical formulas?

    bold : bool, default False
        Use bold formatting for names?

    italic : bool, default False
        Use italic formatting for names?

    font : str, optional
        Font type for names (has no effect if format.names is TRUE).

    family : str, optional
        Font family for names.

    adj : numeric or list, default 0.5
        Adjustment for line labels.

    dx : numeric or list, default 0
        X offset for line or field labels.

    dy : numeric or list, default 0
        Y offset for line or field labels.

    srt : numeric, default 0
        Rotation for line labels.

    min_area : numeric, default 0
        Minimum area of fields that should be labeled, expressed as a fraction
        of the total plot area.

    main : str, optional
        A main title for the plot; None means to plot no title.

    legend_x : str, optional
        Description of legend placement passed to legend.

    add : bool, default False
        Add to current plot?

    plot_it : bool, default True
        Make a plot?

    tplot : bool, default True
        Set up plot with thermo.plot.new?
    
    width, height : numeric, default 600 by 520
        Width and height of the plot.
        
    dpi : numeric, default 150
        Resolution of the plot.
    
    messages : bool, default True
        Display messages from CHNOSZ?
    
    Returns
    -------
    a : rpy2.ListVector
        Output from `diagram`.
    args : dict
        Dictionary of arguments supplied to `diagram`.
    """
    
    args = {'eout':eout, 'ptype':ptype, 'alpha':alpha, 'normalize':normalize,
            'as.residue':as_residue, 'ylog':ylog, 'fill.NA':fill_NA, 'format.names':format_names,
            'bold':bold, 'italic':italic, 'adj':adj, 'dx':dx, 'dy':dy, 'srt':srt,
            'min.area':min_area, 'plot.it':plot_it, 'tplot':tplot}
    
    if balance != None: args["balance"] = balance
    if groups != None: args["groups"] = groups
    if xrange != None: args["xrange"] = _convert_to_RVector(xrange)
    if mar != None: args["mar"] = _convert_to_RVector(mar)
    if yline != None: args["yline"] = _convert_to_RVector(yline)
    if side != None: args["side"] = _convert_to_RVector(side)
    if xlim != None: args["xlim"] = _convert_to_RVector(xlim)
    if ylim != None: args["ylim"] = _convert_to_RVector(ylim)
    if xlab != None: args["xlab"] = xlab
    if ylab != None: args["ylab"] = ylab
    if cex != None: args["cex"] = _convert_to_RVector(cex)
    if cex_names != None: args["cex.names"] = _convert_to_RVector(cex_names)
    if cex_axis != None: args["cex.axis"] = _convert_to_RVector(cex_axis)
    if lty != None: args["lty"] = _convert_to_RVector(lty)
    if lwd != None: args["lwd"] = _convert_to_RVector(lwd)
    if dotted != None: args["dotted"] = convert_to_RVector(dotted)
    if spline_method != None: args["spline.method"] = spline_method
    if contour_method != None: args["contour.method"] = contour_method
    if levels != None: args["levels"] = _convert_to_RVector(levels)
    if col != None: args["col"] = _convert_to_RVector(col)
    if col_names != None: args["col.names"] = _convert_to_RVector(col_names)
    if fill != None: args["fill"] = _convert_to_RVector(fill)
    if limit_water != None: args["limit.water"] = limit_water
    if names != None: args["names"] = _convert_to_RVector(names)
    if font != None: args["font"] = _convert_to_RVector(font)
    if family != None: args["family"] = _convert_to_RVector(family)
    if isinstance(adj, list):
        args["adj"] = _convert_to_RVector(adj)
    if isinstance(dx, list):
        args["dx"] = _convert_to_RVector(dx)
    if isinstance(dy, list):
        args["dy"] = _convert_to_RVector(dy)
    if main != None: args["main"] = main
    if legend_x != None: args["legend.x"] = legend_x
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        with __r_inline_plot(width=width, height=height, dpi=dpi):
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
    Calculate the chemical affinities of formation reactions of species.
    
    Parameters
    ----------
    property : str, default 'A'
        The property to be calculated. Default is A, for chemical affinity of
        formation reactions of species of interest.

    sout : list
        Output from `subcrt`.

    exceed_Ttr : bool, default False
        Allow subcrt to compute properties for phases beyond their transition
        temperature?

    exceed_rhomin : bool, default False
        Allow subcrt to compute properties of species in the HKF model below
        0.35 g cm-3?

    return_buffer : bool, default False
        If TRUE, and a buffer has been associated with one or more basis
        species in the system, return the values of the activities of the basis
        species calculated using the buffer. Default is FALSE.

    return_sout : bool, default False
        Return only the values calculated with subcrt?

    balance : str, default 'PBB'
        This argument is used to identify a conserved basis species (or PBB) in
        a chemical activity buffer.

    iprotein : numeric, optional
        Indices of proteins in thermo$protein for which to calculate
        properties.

    loga_protein : numeric, default -3
        Logarithms of activities of proteins identified in iprotein.

    transect : bool, optional
        Force a transect calculation, even for three or fewer values of the
        variables?
    
    messages : bool, default True
        Display messages from CHNOSZ?
        
    **kwargs : dict
        Numeric, zero or more named arguments, used to identify the variables
        of interest in the calculations.
    
    Returns
    -------
    a : rpy2.ListVector
        Output from `affinity`.
    """
    
    args = {'exceed_Ttr':exceed_Ttr, 'exceed_rhomin':exceed_rhomin,
            'return_buffer':return_buffer, 'return_sout':return_sout,
            'balance':balance, 'loga_protein':loga_protein}
    
    if property != None: args["property"] = property
    if sout != None: args["sout"] = sout
    if iprotein != None: args["iprotein"] = _convert_to_RVector(iprotein)
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

            
def species(species=None, state=None, delete=False, add=False,
            index_return=False, messages=True):
    
    """
    Python wrapper for the species() function in CHNOSZ.
    Define the species of interest in a system; modify their physical states
    and logarithms of activities.
    
    Parameters
    ----------
    species : str, int, or list of str or int
        Names or formulas of species to add to the species definition; int,
        rownumbers of species in OBIGT to modify or delete.

    state : str or list of str, optional
        physical states; numeric, logarithms of activities or fugacities.

    delete : bool, default False
        Delete the species identified by numeric values of species (or all
        species if that argument is missing)?

    add : bool, default False
        Delete a previous species definition instead of adding to it?

    index_return : bool, default False
        return the affected rownumbers of species in the OBIGT database instead
        of the normal output of `species`?
    
    messages : bool, default True
        Display messages from CHNOSZ?
        
    Returns
    ----------
    pd.DataFrame
        Pandas dataframe containing a stoichometric matrix of reactions to form
        species from basis species defined by `basis`.
    """
    
    args={}
    
    if species != None:
        args["species"] = _convert_to_RVector(species, force_Rvec=False)
            
    if state != None:
        args["state"] = _convert_to_RVector(state, force_Rvec=False)
        
    args["add"] = add
    args["delete"] = delete
    args["index.return"] = index_return
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        sout = CHNOSZ.species(**args)
        
    if messages:
        for warning in w:
            print(warning.message)
    
    return pandas2ri.ri2py_dataframe(sout)


def basis(species=None, state=None, logact=None, delete=False, messages=True):
    
    """
    Python wrapper for the basis() function in CHNOSZ.
    Define the basis species of a chemical system.
    
    Parameters
    ----------
    species : str, int, or list of str or int
        Names or formulas of species, or numeric, indices of species.

    state : str or list of str, optional
        Physical states or names of buffers.

    logact : numeric or list of numeric, optional
        Logarithms of activities or fugacities.

    delete : bool, default False
        Delete the current basis species definition?

    messages : bool, default True
        Display messages from CHNOSZ?
        
    Returns
    ----------
    pd.DataFrame
        Pandas dataframe containing a stoichometric matrix of basis species'
        elemental composition.
    """
    
    args={}
    
    if species != None:
        args["species"] = _convert_to_RVector(species, force_Rvec=False)
            
    if state != None:
        args["state"] = _convert_to_RVector(state, force_Rvec=False)
    
    if logact != None: args["logact"] = _convert_to_RVector(logact)
    
    args["delete"] = delete
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        bout = CHNOSZ.basis(**args)
        
    if messages:
        for warning in w:
            print(warning.message)
    
    return pandas2ri.ri2py_dataframe(bout)


def reset(messages=True):
    
    """
    Python wrapper for the reset() function in CHNOSZ.
    Reset all of the data used in CHNOSZ to default values. This includes the
    computational settings, thermodynamic database, and system settings
    (chemical species).
    
    Parameters
    ----------
    messages : bool, default True
        Print messages from CHNOSZ?
    
    """
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        CHNOSZ.reset()
        
    if messages:
        for warning in w:
            print(warning.message)


def add_OBIGT(file, species=None, force=True, messages=True):
    
    """
    Python wrapper for the add.OBIGT() function in CHNOSZ.
    Add or overwrite species in the OBIGT thermodynamic database by supplying
    a comma separated value (csv) file with custom data.
    
    Parameters
    ----------
    file : str
        Path to a file.

    species : str, optional
        Names of species to load from file.

    force : bool, default True
        Force replacement of already existing species?
    
    messages : bool, default True
        Print messages from CHNOSZ?
    
    Returns
    ----------
    a : list of int
        A list of OBIGT database indices (ispecies) that have been added or
        modified.
    """
    
    args={'file':file}
    
    if species != None:
        if not isinstance(species, list):
            args["species"] = species
        else:
            args["species"] = _convert_to_RVector(species)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ispecies = CHNOSZ.add_OBIGT(**args)
    
    if messages:
        for warning in w:
            print(warning.message)
    
    return list(ispecies)


def mod_OBIGT(*args, messages=True, **kwargs):
    
    """
    Python wrapper for the mod.OBIGT() function in CHNOSZ.
    Modify species in the OBIGT thermodynamic database. Optionally, supply a
    Pandas dataframe containing custom data.
    
    Parameters
    ----------
    *args : str, numeric, or Pandas dataframe
        Species names or thermodynamic database index numbers to modify. If a
        Pandas dataframe, database entries to add or modify.
    
    *kwargs : str or numeric
        Properties of species to modify in the thermodynamic database.
    
    messages : bool, default True
        Print messages from CHNOSZ?
    
    Returns
    ----------
    list of int
        A list of OBIGT database indices (ispecies) that have been added or
        modified.
    """

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        if isinstance(args[0], pd.DataFrame):
            arg_list = list(args)
            arg_list[0] = pandas2ri.py2ri(arg_list[0])
            args = tuple(arg_list)
        else:
            pass
        
        ispecies = CHNOSZ.mod_OBIGT(*args, **kwargs)
    
    if messages:
        for warning in w:
            print(warning.message)
    
    return list(ispecies)
    
    
def info(species, state=None, check_it=True, messages=True):
    
    """
    Python wrapper for the info() function in CHNOSZ.
    Search for species by name or formula, retrieve their thermodynamic
    properties and parameters, and add proteins to the thermodynamic database.
    
    Parameters
    ----------
    species : str, int, or list of str or int
        Name or formula of species, or numeric, rownumber of species in the
        OBIGT database.

    state : str or list of str, optional
        State(s) of species.

    check_it : bool, default True
        Check GHS and EOS parameters for self-consistency?
    
    messages : bool, default True
        Print messages from CHNOSZ?
    
    Returns
    ----------
    pd.DataFrame or list of int
        Returns a Pandas dataframe if supplied a list of OBIGT database indices
        (ispecies), or a list of ispecies if given a list of species names or
        formulas.
    """
    
    args = {}
    output_is_df = False
    
    args["species"] = _convert_to_RVector(species, force_Rvec=False)
    if not isinstance(species, list): species = [species]
    if all(isinstance(x, int) for x in species):
        output_is_df = True
    
    if state != None:
        args["state"] = _convert_to_RVector(state, force_Rvec=False)
    
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
           property=["logK", "G", "H", "S", "V", "Cp"],
           T=None, P=None, grid=None,
           convert=True, exceed_Ttr=False, exceed_rhomin=False,
           logact=None, autobalance=True, IS=None, messages=True,
           show=True):
    
    """
    Python wrapper for the subcrt() function in CHNOSZ.
    Calculate the standard molal thermodynamic properties of one or more
    species or a reaction between species as a function of temperature and
    pressure.
    
    Parameters
    ----------
    species : str, int, or list of str or int
        Name or formula of species, or numeric, rownumber of species in the
        OBIGT database.

    coeff : numeric or list of numeric, optional
        Reaction coefficients on species.

    state : str or list of str, optional
        State(s) of species.

    property : str or list of str, optional
        Property(s) to calculate.

    T : numeric or list of numeric, optional
        Temperature(s) of the calculation.

    P : numeric, list of numeric, or str if 'Psat', default 'Psat'
        Pressure(s) of the calculation.

    grid : str, default None
        Type of PxT grid to produce (None, the default, means no gridding).

    exceed_Ttr : bool, default False
        Calculate Gibbs energies of mineral phases and other species beyond
        their transition temperatures?

    exceed_rhomin : bool, default False
        Return properties of species in the HKF model below 0.35 g cm-3?

    logact : numeric or list of numeric, optional
        Logarithms of activities of species in reaction.

    convert : bool, default True
        Are input and output units of T and P those of the user (True) (see
        T_units), or are they Kelvin and bar (False)?

    autobalance : bool, default True
        Attempt to automatically balance reaction with basis species?

    IS : numeric or list of numeric, optional
        Ionic strength(s) at which to calculate adjusted molal properties,
        mol kg^-1.
    
    messages : bool, default True
        Print messages from CHNOSZ?
        
    show : bool, default True
        Display CHNOSZ tables?
    
    Returns
    ----------
    out : object of class SubcrtOutput
        An object that stores the output of `subcrt`.
    """
    
    single_species = False
    
    args = {'species': _convert_to_RVector(species, force_Rvec=False)}
        
    if coeff != None:
        args["coeff"] = _convert_to_RVector(coeff)
    else:
        single_species = True
        
    if state != None:
        args["state"] = _convert_to_RVector(state, force_Rvec=False)
    
    args["property"] = _convert_to_RVector(property)
    
    if T != None:
        args['T'] = _convert_to_RVector(T, force_Rvec=False)

    if P != None:
        args['P'] = _convert_to_RVector(P, force_Rvec=False)
    
    if grid != None: args['grid'] = grid # grid is either 'T' or 'P'
    
    args['convert'] = convert
    args['exceed.Ttr'] = exceed_Ttr
    args['exceed.rhomin'] = exceed_rhomin
    
    if logact != None: args["logact"] = _convert_to_RVector(logact)
    
    args['autobalance'] = autobalance
    
    if IS != None:
        args["IS"] = _convert_to_RVector(IS, force_Rvec=False)
    
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
    
    out = SubcrtOutput(out_dict)
    
    if show:
        for table in out.__dict__.keys():
            if not isinstance(out[table], dict):
                display(out[table])
            else:
                for subtable in out[table].keys():
                    print("\n"+subtable+":") # species name
                    display(out[table][subtable])
    
    return out


class SubcrtOutput(object):
    """
    Stores the output of a `subcrt` calculation.
    
    Attributes
    ----------
    reaction : pd.DataFrame
        Pandas dataframe summary of the reaction. Only present if a reaction is
        specified.
        
    species : pd.Dataframe
        Pandas dataframe summary of species. Only present if a reaction is not
        specified.
        
    out : pd.Dataframe
        Pandas dataframe of `subcrt` output.
        
    warnings : pd.Dataframe
        Pandas dataframe of calculation warnings. Only present if warnings were
        generated.
        
    """
    def __init__(self, args): 
        for k in args:
            setattr(self, k, args[k])

    def __getitem__(self, item):
         return getattr(self, item)


class thermo(object):
    
    """
    Python wrapper for the thermo() object in CHNOSZ.
    See the original CHNOSZ documentation for in-depth descriptions of each
    attribute: https://chnosz.net/manual/thermo.html
    
    Attributes
    ----------
    OBIGT : pd.DataFrame
        A thermodynamic database of standard molal thermodynamic properties and
        equations of state parameters of species.
    
    basis : pd.DataFrame
        Initially `None`, reserved for a dataframe written by basis upon
        definition of the basis species.
    
    buffer : pd.DataFrame
        Contains definitions of buffers of chemical activity.
    
    element : pd.DataFrame
        Containins the thermodynamic properties of elements taken from Cox et
        al., 1989, Wagman et al., 1982, and (for Am, Pu, Np, Cm) Thoenen et al.,
        2014.
    
    groups : pd.DataFrame
        A dataframe with 22 columns for the amino acid sidechain, backbone and
        protein backbone groups ([Ala]..[Tyr],[AABB],[UPBB]) whose rows
        correspond to the elements C, H, N, O, S. It is used to quickly
        calculate the chemical formulas of proteins that are selected using the
        iprotein argument in `affinity`.
    
    opar : dict
        Stores parameters of the last plot generated in pyCHNOSZ. If a plot has
        not yet been generated, `opar` is `None`.
        
    opt : dict
        Dictionary of computational settings.

    protein : pd.DataFrame
        Amino acid compositions of selected proteins.
        
    refs : pd.DataFrame
        References for thermodynamic data.
        
    stoich : pd.DataFrame
        A precalculated stoichiometric matrix for the default database.
        
    species : pd.DataFrame
        Initially `None`, reserved for a dataframe generated by species to
        define the species of interest.
        
    """
    
    def __init__(self, messages=True, **kwargs):
        
        args = {}
        
        for key, value in kwargs.items():
            if isinstance(value, list) or isinstance(value, str) or isinstance(value, NumberTypes):
                value = _convert_to_RVector(value, force_Rvec=False)
            elif isinstance(value, pd.DataFrame):
                value = pandas2ri.py2ri_dataframe(value)
            else:
                pass
            args.update({key:value})
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            t = CHNOSZ.thermo(**args)

        if messages:
            for warning in w:
                print(warning.message)
        
        for i, name in enumerate(t.names):
            if isinstance(t[i], ro.DataFrame) or isinstance(t[i], ro.Matrix):
                attr = pandas2ri.ri2py_dataframe(t[i])
            elif isinstance(t[i], ro.ListVector):
                attr = {}
                for ii, subname in enumerate(t[i].names):
                    attr[subname] = list(t[i][ii])
            elif isinstance(t[i], ro.FloatVector) or isinstance(t[i], ro.IntVector) or isinstance(t[i], ro.StrVector):
                attr = list(t[i])
                if len(t[i]) == 1:
                    attr = attr[0]
            elif t[i] == ro.r("NULL"):
                attr = None
            else:
                attr = t[i]
            
            
            setattr(self, name, attr)
            
    def __getitem__(self, item):
         return getattr(self, item)