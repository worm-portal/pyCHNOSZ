Module pyCHNOSZ.fun
===================

Functions
---------

    
`add_OBIGT(file, species=None, force=True, messages=True)`
:   Python wrapper for the add.OBIGT() function in CHNOSZ.
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

    
`add_protein(aa, messages=True)`
:   Python wrapper for the add.protein() function in CHNOSZ.
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

    
`affinity(property=None, sout=None, exceed_Ttr=False, exceed_rhomin=False, return_buffer=False, return_sout=False, balance='PBB', iprotein=None, loga_protein=-3, transect=None, messages=True, **kwargs)`
:   Python wrapper for the affinity() function in CHNOSZ.
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
        
    **kwargs : numeric, zero or more named arguments
        Used to identify the variables of interest in the calculations. E.g.,
        pH, T, P, Eh, ...
    
    Returns
    -------
    a : rpy2.ListVector
        Output from `affinity`.

    
`animation(basis_args={}, species_args={}, affinity_args={}, equilibrate_args=None, diagram_args={}, anim_var='T', anim_range=[0, 350, 8], save_as='newanimationframe', save_format='png', height=300, width=400, save_scale=1, messages=False)`
:   Produce an animated interactive affinity, activity, or predominance diagram.
    
    Parameters
    ----------
    basis_args : dict
        Dictionary of options for defining basis species (see `basis`) in the
        animated diagram.
        Example: basis_args={'species':['CO2', 'O2', 'H2O', 'H+']}
    
    species_args : dict
        Dictionary of options for defining species (see `species`) in the
        animated diagram, or a list of dicts.
        Example 1: species_args={'species':['CO2', 'HCO3-', 'CO3-2']}
        Example 2: species_args=[
                {'species':['CO2', 'HCO3-', 'CO3-2'], 'state':[-4]},
                {'species':['graphite'], state:[0], 'add':True}]
    
    affinity_args : dict
        Dictionary of options for defining the affinity calculation (see
        `affinity`).
        Example: affinity_args={"pH":[2, 12, 100]}
        Example: affinity_args={"pH":[2, 12, 100], "P":[2000, 4000, 100]}
    
    equilibrate_args : dict or None, default None
        Dictionary of options for defining equilibration calculation
        (see `equilibrate`). If None, plots output from `affinity`.
        Example: equilibrate_args={"balance":1}
    
    diagram_args : dict
        Dictionary of options for diagramming (see `diagram`). Diagram option
        `interactive` is set to True.
        Example: diagram_args={"alpha":True}
    
    anim_var : str, default "T"
        Variable that changes with each frame of animation.
    
    anim_range : list of numeric, default [0, 350, 8]
        The first two numbers in the list are the starting and ending
        values for `anim_var`. The third number in the list is the desired
        number of animation frames.
    
    messages : bool, default True
        Display messages from CHNOSZ?
    
    Returns
    -------
    An interactive animated plot.

    
`basis(species=None, state=None, logact=None, delete=False, messages=True)`
:   Python wrapper for the basis() function in CHNOSZ.
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

    
`diagram(eout, ptype='auto', alpha=False, normalize=False, as_residue=False, balance=None, groups=None, xrange=None, mar=None, yline=None, side=[1, 2, 3, 4], ylog=True, xlim=None, ylim=None, xlab=None, ylab=None, cex=None, cex_names=None, cex_axis=None, lty=None, lwd=None, dotted=None, spline_method=None, contour_method=None, levels=None, col=None, col_names=None, fill=None, fill_NA='gray80', limit_water=None, names=None, format_names=True, bold=False, italic=False, font=None, family=None, adj=0.5, dx=0, dy=0, srt=0, min_area=0, main=None, legend_x=None, add=False, plot_it=True, tplot=True, annotation=None, annotation_coords=[0, 0], width=600, height=520, dpi=150, messages=True, interactive=False, save_as=None, save_format=None, save_scale=1, fig_out=False)`
:   Python wrapper for the diagram() function in CHNOSZ.
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
    
    annotation : str, optional
        Annotation to add to the plot. Interactive plots only (`interactive`
        is set to True).
    
    annotation_coords : list of numeric, default [0, 0], optional
        Coordinates of annotation, where 0,0 is bottom left and 1,1 is top
        right. Interactive plots only (`interactive` is set to True).
    
    width, height : numeric, default 600 by 520
        Width and height of the plot.
        
    dpi : numeric, default 150
        Resolution of the plot.
    
    messages : bool, default True
        Display messages from CHNOSZ?
    
    interactive : bool, default False
        Display an interactive plot?
    
    save_as : str, optional
        For interactive plots only (`interactive`=True). Provide a filename to
        save this figure. Filetype of saved figure is determined by
        `save_format`.
        Note: interactive plots can be saved by clicking the 'Download plot'
        button in the plot's toolbar.
    
    save_format : str, default "png"
        For interactive plots only (`interactive`=True). Desired format of saved
        or downloaded figure. Can be 'png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf',
        'eps', 'json', or 'html'. If 'html', an interactive plot will be saved.
        Only 'png', 'svg', 'jpeg', and 'webp' can be downloaded with the
        'download as' button in the toolbar of an interactive plot.
    
    save_scale : numeric, default 1
        For interactive plots only (`interactive`=True). Multiply
        title/legend/axis/canvas sizes by this factor when saving the figure.
        
    fig_out : bool, default False
        Function output is a plotly figure? Ignored if `interactive` is False.
        
    Returns
    -------
    a : rpy2.ListVector
        Output from `diagram`.
    args : dict
        Dictionary of arguments supplied to `diagram`.

    
`diagram_interactive(data, title=None, borders=0, names=None, annotation=None, annotation_coords=[0, 0], balance=None, xlab=None, ylab=None, colormap='viridis', width=600, height=520, alpha=False, messages=True, plot_it=True, save_as=None, save_format=None, save_scale=1)`
:   Produce an interactive diagram.
    
    Parameters
    ----------
    data : rpy2.ListVector
        Output from `equilibrate` or `affinity`.
    
    title : str, optional
        Title of the plot.
    
    borders : float, default 0
        If set to a value greater than 0, shows lines that indicate borders
        between regions in predominance diagrams. Value indicates thickness
        of the border, in pixels.
    
    names : str, optional
        Names of species for activity lines or predominance fields.
    
    annotation : str, optional
        Annotation to add to the plot.
    
    annotation_coords : list of numeric, default [0, 0], optional
        Coordinates of annotation, where 0,0 is bottom left and 1,1 is top
        right.
        
    balance : str or numeric, optional
        How to balance the transformations.
    
    xlab, ylab : str
        Custom x and y axes labels.
    
    width, height : numeric, default 600 by 520
        Width and height of the plot.
        
    alpha : bool or str (balance), default False
        For speciation diagrams, plot degree of formation instead of
        activities?
    
    messages : bool, default True
        Display messages from CHNOSZ?
    
    plot_it : bool, default True
        Show the plot?
    
    save_as : str, optional
        For interactive plots only (`interactive`=True). Provide a filename to
        save this figure. Filetype of saved figure is determined by
        `save_format`.
        Note: interactive plots can be saved by clicking the 'Download plot'
        button in the plot's toolbar.
    
    save_format : str, default "png"
        For interactive plots only (`interactive`=True). Desired format of saved
        or downloaded figure. Can be 'png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf',
        'eps', 'json', or 'html'. If 'html', an interactive plot will be saved.
        Only 'png', 'svg', 'jpeg', and 'webp' can be downloaded with the
        'download as' button in the toolbar of an interactive plot.
    
    save_scale : numeric, default 1
        For interactive plots only (`interactive`=True). Multiply
        title/legend/axis/canvas sizes by this factor when saving the figure.
    
    Returns
    -------
    An interactive plot.

    
`entropy(formula, messages=True)`
:   Python wrapper for the entropy() function in CHNOSZ.
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

    
`equilibrate(aout, balance=None, loga_balance=None, ispecies=None, normalize=False, messages=True)`
:   Python wrapper for the equilibrate() function in CHNOSZ.
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

    
`html_chemname_format(name)`
:   Format a chemical formula to display subscripts and superscripts in HTML
    (e.g., Plotly plots)
    Example, "CH3COO-" becomes "CH<sub>3</sub>COO<sup>-</sup>"
    
    Parameters
    ----------
    name : str
        A chemical formula.
    
    Returns
    -------
    A formatted chemical formula string.

    
`info(species, state=None, check_it=True, messages=True)`
:   Python wrapper for the info() function in CHNOSZ.
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

    
`makeup(formula, multiplier=1, sum=False, count_zero=False, messages=True)`
:   Python wrapper for the makeup() function in CHNOSZ.
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

    
`mass(formula, messages=True)`
:   Python wrapper for the mass() function in CHNOSZ.
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

    
`mod_OBIGT(*args, messages=True, **kwargs)`
:   Python wrapper for the mod.OBIGT() function in CHNOSZ.
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

    
`ratlab(top='K+', bottom='H+', molality=False)`
:   Python wrapper for the ratlab() function in CHNOSZ.
    Produces a expression for the activity ratio between the ions in the top and
    bottom arguments. The default is a ratio with H+, i.e.
    (activity of the ion) / [(activity of H+) ^ (charge of the ion)]
    
    Parameters
    ----------
    top : str, default "K+"
        The ion in the numerator of the ratio.
    
    bottom : str, default "H+"
        The ion in the denominator of the ratio.
    
    molality : bool, default False
        Use molality (m) instead of activity (a) for aqueous species?
    
    Returns
    -------
    A formatted string representing the activity ratio.

    
`reset(db='OBIGT', messages=True)`
:   Python wrapper for the reset() function in CHNOSZ.
    Reset all of the data used in CHNOSZ to default values. This includes the
    computational settings, thermodynamic database, and system settings
    (chemical species).
    
    Parameters
    ----------
    db : str, default "OBIGT"
        Accepts "WORM" or "OBIGT". Which thermodynamic database should be used?
        If "WORM", the database from
        https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv
        will be loaded. Otherwise, the default database for CHNOSZ, called
        OBIGT, will be loaded.
    
    messages : bool, default True
        Print messages from CHNOSZ?

    
`retrieve(elements=None, ligands=None, state=None, T=None, P='Psat', add_charge=True, hide_groups=True, must_have=None, messages=True)`
:   Python wrapper for the retrieve() function in CHNOSZ.
    Retrieve species in the database containing one or more chemical elements.
    
    Parameters
    ----------
    elements : str, list of str, or tuple of str
        Elements in a chemical system. If `elements` is a string, retrieve
        species containing that element.
        
        E.g., `retrieve("Au")` will return all species containing Au.
        
        If `elements` is a list, retrieve species that have all of the elements
        in the list.
        
        E.g., `retrieve(["Au", "Cl"])` will return all species that have both
        Au and Cl.
        
        If `elements` is a tuple, retrieve species relevant to the system,
        including charged species.
        
        E.g., `retrieve(("Au", "Cl"))` will return species that have Au
        and/or Cl, including charged species, but no other elements.
    
    ligands : str, list of str, or tuple of str, optional
        Elements present in any ligands.
        
    state : str, list of str, or tuple of str, optional
        Filter the result on these state(s).
        
    T : float, optional
        Temperature where DeltaG0 of species must not be NA
    
    P : float or "Psat", default "Psat"
        Pressure where DeltaG0 of species must not be NA
    
    add_charge : bool, default True
        Add charge to the system?
    
    hide_groups : bool, default True
        Exclude groups from the result?
    
    must_have : str or list of str, optional
        Retrieved species must have the element(s).
    
    messages : bool, default True
        Display messages from CHNOSZ?
    
    Returns
    ----------
    out: list or dict
        A list of the OBIGT indices of retrieved chemical species

    
`seq2aa(protein, sequence, messages=True)`
:   Python wrapper for the seq2aa() function in CHNOSZ.
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

    
`solubility(iaq=None, in_terms_of=None, dissociate=False, find_IS=False, messages=True, **kwargs)`
:   Python wrapper for the solubility() function in CHNOSZ.
    Calculate chemical activities of aqueous species in equilibrium with a mineral or gas.
    
    Parameters
    ----------
    iaq : str, int, or list of str or int
        Name(s) of aqueous species (if str). If int, the index of the aqueous
        species in the OBIGT database.
    
    in_terms_of : str, optional
        Express the total solubility in terms of moles of this species.
    
    dissociate : bool, default False
        Does the mineral undergo a dissociation reaction?
        
    find_IS : bool, default False
        Find the equilibrium ionic strength by iteration?
    
    messages : bool, default True
        Display messages from CHNOSZ?
        
    **kwargs : named arguments
        Arguments for `affinity` or `mosaic` (i.e. plotting variables).
    
    Returns
    -------
    s : rpy2.ListVector
        Output from `solubility`.

    
`species(species=None, state=None, delete=False, add=False, index_return=False, messages=True)`
:   Python wrapper for the species() function in CHNOSZ.
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

    
`subcrt(species, coeff=None, state=None, property=['logK', 'G', 'H', 'S', 'V', 'Cp'], T=None, P=None, grid=None, convert=True, exceed_Ttr=False, exceed_rhomin=False, logact=None, autobalance=True, IS=None, messages=True, show=True)`
:   Python wrapper for the subcrt() function in CHNOSZ.
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

    
`syslab(system=['K2O', 'Al2O3', 'SiO2', 'H2O'], dash='-')`
:   Python wrapper for the syslab() function in CHNOSZ.
    Formats the given thermodynamic components and adds intervening en dashes.
    
    Parameters
    ----------
    system : list of str, default ["K2O", "Al2O3", "SiO2", "H2O"]
        Thermodynamic components.
    
    dash : str, default "-"
        Character to use for dash between components.
    
    Returns
    -------
    A formatted string representing the thermodynamic system.

    
`unicurve(logK, species, coeff, state, pressures=1, temperatures=25, IS=0, minT=0.1, maxT=100, minP=1, maxP=500, tol=None, solve='T', width=600, height=520, dpi=90, plot_it=True, messages=True, show=True)`
:   Solve for temperatures or pressures of equilibration for a given logK
    value and produce a plot.
    
    Parameters
    ----------
    logK : numeric
        Logarithm (base 10) of an equilibrium constant.
    
    species : str, int, or list of str or int
        Name or formula of species, or numeric, rownumber of species in the
        OBIGT database.
    
    coeff : numeric or list of numeric, optional
        Reaction coefficients on species.
        
    state : str or list of str, optional
        State(s) of species.
    
    pressures, temperatures : numeric or list of numeric
        Pressures (if solving for temperature) or temperatures (if solving for
        pressures) to resolve with `logK`. Pressures are in bars and
        temperatures are in degrees Celcius.
    
    IS : numeric or list of numeric, optional
        Ionic strength(s) at which to calculate adjusted molal properties,
        mol kg^-1.
    
    minT, maxT : numeric
        Minimum and maximum temperatures (degrees Celcius) to search within for
        a solution for the given pressure(s) and logK. Ignored if solving for
        pressure.
    
    minP, maxP : numeric
        Minimum and maximum pressures (bars) to search within for
        a solution for the given temperatures(s) and logK. Ignored if solving
        for temperature.
    
    tol : float
        Tolerance for finding a temperature or pressure that converges on the
        given logK. Will attempt to find a solution that satisfies logK plus or
        minus `tol`. By default, tol is equal to 1/(10^(n+2)) where n is the
        number of decimal places in logK, with a maximum default tolerance of
        1e-5.
    
    solve : "T" or "P"
        Solve for temperature or pressure?
    
    width, height : numeric, default 600 by 520
        Width and height of the plot.
        
    dpi : numeric, default 90
        Resolution of the plot.
    
    plot_it : bool, default True
        Show the plot?
    
    messages : bool, default True
        Display messages from CHNOSZ?
    
    show : bool, default True
        Display CHNOSZ tables?
    
    Returns
    -------
    out : object of class SubcrtOutput
        An object that stores the output of `subcrt` along the univariant curve.

    
`univariant_TP(logK, species, coeff, state, Trange, Prange, IS=0, xlim=None, ylim=None, line_type='markers+lines', tol=None, title=None, res=10, width=500, height=400, save_as=None, save_format=None, save_scale=1, show=False, messages=False, plot_it=True)`
:   Solve for temperatures or pressures of equilibration for a given logK
    value and produce a plot with temperature and pressure axes.
    
    Parameters
    ----------
    logK : numeric
        Logarithm (base 10) of an equilibrium constant.
    
    species : str, int, or list of str or int
        Name or formula of species, or numeric, rownumber of species in the
        OBIGT database.
    
    coeff : numeric or list of numeric, optional
        Reaction coefficients on species.
    
    state : str or list of str, optional
        State(s) of species.
    
    Trange : list of two numeric
        List containing the minimum and maximum temperature, in degrees C, to
        solve for the specified logK value. Does not necessarily correspond to
        temperature axis range in the resulting plot.
    
    Prange : list of two numeric
        List containing the minimum and maximum pressure, in bars, to
        solve for the specified logK value. Does not necessarily correspond to
        pressure axis range in the resulting plot.
    
    IS : numeric or list of numeric, optional
        Ionic strength(s) at which to calculate adjusted molal properties,
        mol kg^-1.
    
    xlim, ylim : list of numeric, optional
        Define the range of the x and y axes, respectively. For example,
        `xlim=[0, 400]` will set the x axis range from 0 to 400°C.
    
    line_type : str, default "markers+lines"
        What type of line should be plotted? Options are "markers+lines",
        "markers", and "lines".
    
    tol : float
        Tolerance for finding a temperature or pressure that converges on the
        given logK. Will attempt to find a solution that satisfies logK plus or
        minus `tol`. By default, tol is equal to 1/(10^(n+2)) where n is the
        number of decimal places in logK, with a maximum default tolerance of
        1e-5.
    
    width, height : numeric, default 500 by 400
        Width and height of the plot.
    
    save_as : str, optional
        For interactive plots only (`interactive`=True). Provide a filename to
        save this figure. Filetype of saved figure is determined by
        `save_format`.
        Note: interactive plots can be saved by clicking the 'Download plot'
        button in the plot's toolbar.
    
    save_format : str, default "png"
        For interactive plots only (`interactive`=True). Desired format of saved
        or downloaded figure. Can be 'png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf',
        'eps', 'json', or 'html'. If 'html', an interactive plot will be saved.
        Only 'png', 'svg', 'jpeg', and 'webp' can be downloaded with the
        'download as' button in the toolbar of an interactive plot.
    
    save_scale : numeric, default 1
        For interactive plots only (`interactive`=True). Multiply
        title/legend/axis/canvas sizes by this factor when saving the figure.
    
    messages : bool, default True
        Display messages from CHNOSZ?
    
    show : bool, default True
        Display CHNOSZ tables?
    
    plot_it : bool, default True
        Show the plot?
    
    Returns
    -------
    out : object of class SubcrtOutput
        An object that stores the output of `subcrt` along the univariant curve.

    
`water(property=None, T=298.15, P='Psat', P1=True, messages=True)`
:   Python wrapper for the water() function in CHNOSZ.
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

    
`zc(formula, messages=True)`
:   Python wrapper for the ZC() function in CHNOSZ.
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

Classes
-------

`R_output()`
:   

    ### Methods

    `capture_r_output(self)`
    :   Capture and create a list of R console messages

`SubcrtOutput(args)`
:   Stores the output of a `subcrt` calculation.
    
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

`thermo(db='OBIGT', messages=True, **kwargs)`
:   Python wrapper for the thermo() object in CHNOSZ.
    See the original CHNOSZ documentation for in-depth descriptions of each
    attribute: https://chnosz.net/manual/thermo.html
    
    Parameters
    ----------
    db : str, default "OBIGT"
        Accepts "WORM" or "OBIGT". Which thermodynamic database should be used?
        If "WORM", the database from
        https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv
        will be loaded. Otherwise, the default database for CHNOSZ, called
        OBIGT, will be loaded.
    
    messages : bool, default True
        Print messages from CHNOSZ?
    
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