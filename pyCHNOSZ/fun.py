import os
from contextlib import contextmanager
from IPython.display import Image, display
import pandas as pd
import numpy as np
import re
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import copy
import statistics
import decimal
import chemparse
import copy

from urllib.request import urlopen
from io import StringIO
import time

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)
    import rpy2.robjects as ro
    import rpy2.rinterface_lib
    from rpy2.robjects.packages import importr
    from rpy2.robjects.lib import grdevices

from WORMutils import chemlabel, can_connect_to, import_package_file
from wormutils_r import R_output, pd_to_r_df, r_df_to_pd, rpy2float

import matplotlib.pyplot as plt

CHNOSZ = importr("CHNOSZ")
grdev = importr('grDevices')

NumberTypes = (int, float, complex)


@contextmanager
def __r_inline_plot(width=600, height=520, dpi=150, plot_it=True):
    
    """
    Display R plots inline.

    Parameters
    ----------
    width, height : numeric, default 600 by 520
        Width and height of the plot.
        
    dpi : numeric, default 150
        Resolution of the plot.

    plot_it : bool, default True
        Render the plot?
    """
    
    with grdevices.render_to_bytesio(grdevices.png, 
                                     width=width,
                                     height=height, 
                                     res=dpi) as b:

        yield

    data = b.getvalue()
    if plot_it:
        display(Image(data=data, format='png', embed=True))


def __flatten_list(_2d_list):
    flat_list = []
    # Iterate through the outer list
    for element in _2d_list:
        if type(element) is list:
            # If the element is of type list, iterate through the sublist
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list


def __save_figure(fig, save_as, save_format, save_scale, plot_width, plot_height, ppi):

    if isinstance(save_format, str) and save_format not in ['png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'eps', 'json', 'html']:
        raise Exception("{}".format(save_format)+" is an unrecognized "
                        "save format. Supported formats include 'png', "
                        "'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'eps', "
                        "'json', or 'html'")

    if isinstance(save_format, str) and save_as != None:
        if not isinstance(save_as, str):
            save_as = "newplot"
        if save_format=="html":
            fig.write_html(save_as+".html")
            print("Saved figure as {}".format(save_as)+".html")
            save_format = 'png'
        elif save_format in ['pdf', 'eps', 'json']:
            pio.write_image(fig, save_as+"."+save_format, format=save_format, scale=save_scale,
                            width=plot_width*ppi, height=plot_height*ppi)
            print("Saved figure as {}".format(save_as)+"."+save_format)
            save_format = "png"
        else:
            pio.write_image(fig, save_as+"."+save_format, format=save_format, scale=save_scale,
                            width=plot_width*ppi, height=plot_height*ppi)
            print("Saved figure as {}".format(save_as)+"."+save_format)
    else:
        save_format = "png"

    return save_as, save_format


def __seq(start, end, by=None, length_out=None):
    
    """
    Mimic the seq() function in base R.
    """
    
    len_provided = True if (length_out is not None) else False
    by_provided = True if (by is not None) else False
    if (not by_provided) & (not len_provided):
        raise ValueError('At least by or n_points must be provided')
    width = end - start
    eps = pow(10.0, -14)
    if by_provided:
        if (abs(by) < eps):
            raise ValueError('by must be non-zero.')
        absby = abs(by)
        if absby - width < eps: 
            length_out = int(width / absby)
        else: 
            # by is too great, we assume by is actually length_out
            length_out = int(by)
            by = width / (by - 1)
    else:
        length_out = int(length_out)
        by = width / (length_out - 1) 
    out = [float(start)]*length_out
    for i in range(1, length_out):
        out[i] += by * i
    if abs(start + by * length_out - end) < eps:
        out.append(end)
    return out


def convert(value, units, T=298.15, P=1, pH=7, logaH2O=0, messages=True):
    """
    Python wrapper for the convet() function in CHNOSZ.
    Convert values between units.

    Parameters
    ----------
    units : str
        Name of units to set or convert to/from.

    value : float
        Value(s) to be converted.
    
    T : float, default 298.15
        Temperature (Kelvin), used in G⁠-⁠logK, pe-Eh⁠ and ⁠logfO2⁠-E0⁠ conversions.
    
    P : float, default 1
        Pressure (bar), used in logfO2⁠-E0⁠ conversions.
    
    pH : float, default 7
        pH, used in logfO2⁠-⁠E0⁠ conversions.
    
    logaH2O : float, default 0
        Logarithm of activity of water, used in ⁠logfO2⁠-E0⁠ conversions

    messages : bool, default True
        Display messages from CHNOSZ?

    Returns
    -------
    out : float
        Converted value
    
    """

    if isinstance(value, list):
        raise Exception("The convert() function in pyCHNOSZ does not currently "
                "support lists of values, nor the output from solubility().")
    
    args = {"value":value, "units":units, "T":T, "P":P, "pH":pH, "logaH2O":logaH2O}

    capture = R_output()
    capture.capture_r_output()

    out = CHNOSZ.convert(**args)

    if messages:
        capture.print_captured_r_output()

    return out[0]


def add_saturation_lines(asat, d, line_color=None, line_type=None, line_width=1,
                         unique_styles=True, showlegend=True, messages=True):
    """
    Add saturation lines to an interactive plot produced by `diagram`.

    Parameters
    ----------
    asat : output from `affinity`
        Output from `affinity`.

    d : plotly figure object
        Interactive plot from `diagram`.
    
    line_color : str or list of str, optional
        Color of saturation lines. Can be a color name, hex code, or rgb (e.g.,
        "cornflowerblue", "#6495ED", or "rgb(100, 149, 237)")
        or a list of color names, hex codes, or rgb values to cycle through. Note
        that line types defined by line_type are cycled through first before
        moving on to a new color in a list provided to line_color.

    line_type : str or list of str, optional
        The type (style) of saturation lines. Can be 'dot', 'dash', 'dashdot',
        'longdash', 'longdashdot', or 'solid', or a list of line types to cycle
        through. Note that line types are cycled through first before
        moving on to a new color in a list provided to `line_color`.

    line_width : float, default 1
        The width of the saturation lines.
    
    messages : bool, default True
        Print messages produced by this function?

    Returns
    -------
    dout : plotly figure object
        A diagram with saturation lines.
    
    """
    
    dout = copy.deepcopy(d)

    x_var = asat.rx2("args").names[0]
    y_var = asat.rx2("args").names[1]

    x_range = asat.rx2("args").rx2(str(x_var))
    y_range = asat.rx2("args").rx2(str(y_var))

    if len(x_range) == 2:
        res = 256
    else:
        res = int(x_range[2])
    if len(y_range) == 2:
        res = 256
    else:
        res = int(y_range[2])
    
    x_vals = np.linspace(start=x_range[0], stop=x_range[1], num=res)
    y_vals = np.linspace(start=y_range[0], stop=y_range[1], num=res)

    sp_idx = asat.rx2("values").names # species indices as strings
    sp_names = list(info([int(name) for name in sp_idx])["name"])

    if isinstance(line_type, str):
        line_styles = [line_type]
    elif isinstance(line_type, list):
        line_styles = line_type
    else:
        line_styles = ['dot', 'dash', 'dashdot', 'longdash', 'longdashdot']

    if not isinstance(line_color, str) and not isinstance(line_color, list):
        line_color = ["#000000", "#6495ED", "#ED95B4", "#9AD559", "#BF6C8A"]
    
    if isinstance(line_color, str):
        line_color = [line_color]


    sp_idx_approved = []
    p_storage = []
    for i,isp in enumerate(sp_idx):
        m = r_df_to_pd(asat.rx2("values").rx2(str(isp))).T
        cs = plt.contour(x_vals, y_vals, m, [0]) # contour for affinity=0
        plt.close()
        try:
            p = cs.get_paths()[0]
            p_storage.append(p)
            sp_idx_approved.append(isp)
        except:
            if messages:
                print("'"+sp_names[i]+"' could not be included because affinity=0 for "
                      "its formation reaction from basis species could not be found "
                      "within the bounds of the plot.")

    if len(sp_idx_approved) == 0:
        if messages:
            print("Saturation lines for species of interest could not be plotted because "
                  "affinity=0 for their formation reactions from basis species could not be "
                  "found within the bounds of the plot.")
        return d
    
    if unique_styles:
        if len(sp_idx_approved) > len(line_color)*len(line_styles):
                raise Exception("There are more mineral saturation lines than can "
                        "be supported by line styles and colors. Set line_type "
                        "and line_color to lists with more styles and colors.")
    
    sp_names_approved = list(info([int(name) for name in sp_idx_approved])["name"])
    
    line_style_counter = 0
    line_color_counter = 0
    for i,isp in enumerate(sp_idx_approved):

        v = p_storage[i].vertices
        x = v[:,0]
        y = v[:,1]

        dout.add_scatter(x=x, y=y, mode='lines', name=sp_names_approved[i])
        dout.update_traces(line_color=line_color[line_color_counter],
                           selector=dict(name=sp_names_approved[i]), showlegend=showlegend,
                           line=dict(dash=line_styles[line_style_counter],
                                     width=line_width))
        if unique_styles:
            line_style_counter += 1
    
            if line_style_counter+1 > len(line_styles):
                line_style_counter = 0
                line_color_counter += 1

    return dout


def solubility(iaq=None, in_terms_of=None, dissociate=False, find_IS=False,
               messages=True, **kwargs):

    """
    Python wrapper for the solubility() function in CHNOSZ.
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
    """
    
    args = {'dissociate':dissociate, 'find_IS':find_IS}
    
    if in_terms_of != None: args["in_terms_of"] = in_terms_of
    args["iaq"] = _convert_to_RVector(iaq, force_Rvec=True)
    
    for key, value in kwargs.items():
        if isinstance(value, list):
            value = ro.FloatVector(value)
        args.update({key:value})

    capture = R_output()
    capture.capture_r_output()

    s = CHNOSZ.solubility(**args)

    if messages:
        capture.print_captured_r_output()

    return s
    

def retrieve(elements=None, ligands=None, state=None, T=None, P="Psat",
             add_charge=True, hide_groups=True, must_have=None, messages=True):
    """
    Python wrapper for the retrieve() function in CHNOSZ.
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
    """

    if elements == None:
        elements = ro.r("NULL")
    elif isinstance(elements, list) or isinstance(elements, str):
        elements = _convert_to_RVector(elements, force_Rvec=True)
    elif isinstance(elements, tuple):
        elements = ro.ListVector({chr(ord('`')+i):l for i,l in zip(range(1, len(elements)+1), elements)})
    else:
        pass

    if ligands == None:
        ligands = ro.r("NULL")
    elif isinstance(ligands, list) or isinstance(ligands, str):
        ligands = _convert_to_RVector(ligands, force_Rvec=True)
    elif isinstance(ligands, tuple):
        ligands = ro.ListVector({chr(ord('`')+i):l for i,l in zip(range(1, len(ligands)+1), ligands)})
    else:
        pass

    if state == None:
        state = ro.r("NULL")
    elif isinstance(state, list) or isinstance(state, str):
        state = _convert_to_RVector(state, force_Rvec=True)
    elif isinstance(state, tuple):
        state = ro.ListVector({chr(ord('`')+i):l for i,l in zip(range(1, len(state)+1), state)})
    else:
        pass

    if T == None:
        T = ro.r("NULL")
    if P == None:
        P = ro.r("NULL")
    
    args = {'elements':elements, 'ligands':ligands, 'state':state, 'T':T, 'P':P,
            'add_charge':add_charge, 'hide_groups':hide_groups}
    
    capture = R_output()
    capture.capture_r_output()

    out = CHNOSZ.retrieve(**args)

    if messages:
        capture.print_captured_r_output()

    out = list(out)
    
    if isinstance(must_have, str):
        must_have = [must_have]
        
    if isinstance(must_have, list):
        df = info(out, messages=False)
        formulas = df["formula"].apply(chemparse.parse_formula)
        keep_ind = []
        for i,formula in enumerate(formulas):
            if all(elem in formula.keys()  for elem in must_have):
                keep_ind.append(i)
        keep_ind = list(set(keep_ind))
        out = [out[i] for i in keep_ind]

    return out

    
    
def animation(basis_args={}, species_args={}, affinity_args={},
              equilibrate_args=None, diagram_args={},
              anim_var="T", anim_range=[0, 350, 8], xlab=None, ylab=None,
              save_as="newanimationframe", save_format="png", height=300,
              width=400, save_scale=1,
              messages=False):
    
    """
    Produce an animated interactive affinity, activity, or predominance diagram.
    
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

    xlab, ylab : str, optional
        Custom names for the X and Y axes.
    
    messages : bool, default True
        Display messages from CHNOSZ?
    
    Returns
    -------
    An interactive animated plot.
    """
    
    # cap number of frames in animation. Remove limitation after more testing.
    if isinstance(anim_range, list):
        if len(anim_range) == 3:
            if anim_range[2] > 30:
                raise Exception("anim_range is limited to 30 frames.")
        else:
            raise Exception("anim_range must be a list with three values: starting "
                            "value of anim_var, stopping value, and number of "
                            "frames in the animation")
    else:
        raise Exception("anim_range must be a list with three values: starting "
                        "value of anim_var, stopping value, and number of "
                        "frames in the animation")

    if isinstance(basis_args, dict):
        if "species" not in basis_args.keys():
            raise Exception("basis_args needs to contain a list of species for 'species'. "
                            "Example: basis_args={'species':['CO2', 'O2', 'H2O', 'H+']}")
    else:
        raise Exception("basis_args needs to be a Python dictionary with a key "
                        "called 'species' (additional keys are optional). "
                        "Example: basis_args={'species':['CO2', 'O2', 'H2O', 'H+']}")

    
    basis_out = basis(**basis_args)
    basis_sp = list(basis_out.index)
    basis_state = list(basis_out["state"])
    
    if isinstance(species_args, dict):
        if "species" not in species_args.keys():
            raise Exception("species_args needs to contain a list of species for 'species'. "
                            "Example: species_args={'species':['CO2', 'HCO3-', 'CO3-2']}")
        species_args_list = [species_args]
    elif isinstance(species_args, list):
        species_args_list = species_args
        for species_args in species_args_list:
            if "species" not in species_args.keys():
                raise Exception("species_args needs to contain a list of species for 'species'. "
                                "Example: species_args={'species':['CO2', 'HCO3-', 'CO3-2']}")
    else:
        raise Exception("species_args needs to be either a Python dictionary with a key "
                        "called 'species' (additional keys are optional). "
                        "Example: species_args={'species':['CO2', 'HCO3-', 'CO3-2']}"
                        "or else species_args needs to be a list of Python dictionaries."
                        "Example: species_args=[{'species':['CO2', 'HCO3-', 'CO3-2'], 'state':[-4]},"
                        "{'species':['graphite'], state:[0], 'add':True}]")

    # There may be multiple arguments passed to species, especially in cases
    # where add=True. Loop through all the arguments to apply them.
    for species_args in species_args_list:
        if "logact" in species_args.keys():
            mod_species_logact = copy.copy(species_args['logact'])
            del species_args['logact']
        else:
            mod_species_logact = []
    
        species_out = species(**species_args)
        
        if len(mod_species_logact)>0:
            for i in range(0, len(mod_species_logact)) :
                species_out = species(species_args["species"][i], mod_species_logact[i])

    sp = list(species_out["name"])

    if isinstance(sp[0], (int, np.integer)):
        sp = [info(s, messages=False)["name"].values[0] for s in sp]

    dfs = []
    dmaps = []
    dmaps_names = []
    
    if len(anim_range) == 2:
        anim_res = 8
        anim_range = anim_range + [anim_res]
    elif len(anim_range) == 3:
        anim_res = anim_range[2]
        anim_range = [anim_range[0], anim_range[1]]
    
    zvals = __seq(anim_range[0], anim_range[1], length_out=anim_res)
        
    if "messages" not in affinity_args.keys():
        affinity_args["messages"] = messages
    if "messages" not in diagram_args.keys():
        diagram_args["messages"] = messages
    if "plot_it" not in diagram_args.keys():
        diagram_args["plot_it"] = False
    diagram_args["interactive"] = True
    if "format_names" not in diagram_args.keys():
        format_names=True
        format_x_names=True
        format_y_names=True
    
    for z in zvals:

        if anim_var in basis_out.index:
            basis_out = basis(anim_var, z)
        elif anim_var in list(species_out["name"]):
            species_out = species(anim_var, -z)
        elif anim_var == "pH":
            basis_out = basis("H+", -z)
        else:
            affinity_args[anim_var] = z
        
        aeout = affinity(**affinity_args)
        
        if equilibrate_args != None:
            equilibrate_args["aout"] = aeout
            if "messages" not in equilibrate_args.keys():
                equilibrate_args["messages"] = messages
            aeout = equilibrate(**equilibrate_args)
        
        aeout_args = aeout.rx2("args")
        xvar = aeout_args.names[0]
        xrange = list(aeout_args[0])

        res_default = 256 # default affinity resolution
        if len(xrange) == 3:
            xres = int(xrange[2])
        else:
            xres = res_default
        
        diagram_args["eout"] = aeout
        
        df = diagram(**diagram_args)
        
        df[anim_var] = z
        dfs.append(df)

        if 'pred' not in df.columns:
            # affinity/activity plot
            is_predom_plot = False
            
        else:
            # predominance plot
            is_predom_plot = True
            yvar = aeout_args.names[1]
            yrange = list(aeout_args[1])
            if len(yrange) == 3:
                yres = int(yrange[2])
            else:
                yres = res_default

            data = np.array(df.pred)
            shape = (xres, yres)
            dmap = data.reshape(shape)
            dmaps.append(dmap)

            data = np.array(df.prednames)
            shape = (xres, yres)
            dmap_names = data.reshape(shape)
            dmaps_names.append(dmap_names)
        
    xvals = __seq(xrange[0], xrange[1], length_out=xres)


    unit_dict = {"P":"bar", "T":"°C", "pH":"", "Eh":"volts", "IS":"mol/kg"}
    
    if any([anim_var in basis_out.index, anim_var in list(species_out["name"])]) and anim_var not in unit_dict.keys():
        unit_dict[anim_var] = "logact "+anim_var

    for i,s in enumerate(basis_sp):
        if basis_state[i] in ["aq", "liq", "cr"]:
            if format_names:
                unit_dict[s] = "log <i>a</i><sub>{}</sub>".format(chemlabel(s))
            else:
                unit_dict[s] = "log <i>a</i><sub>{}</sub>".format(s)
        else:
            if format_names:
                unit_dict[s] = "log <i>f</i><sub>{}</sub>".format(chemlabel(s))
            else:
                unit_dict[s] = "log <i>f</i><sub>{}</sub>".format(s)

    xlab = xvar+", "+unit_dict[xvar]
    
    if xvar in basis_sp:
        xlab = unit_dict[xvar]
    if xvar == "pH":
        xlab = "pH"
    
    if is_predom_plot:
        ylab = yvar+", "+unit_dict[yvar]
        if yvar in basis_sp:
            ylab = unit_dict[yvar]
        if yvar == "pH":
            yvar = "pH"

    
    if not is_predom_plot:

        if 'loga.equil' not in aeout.names:
            yvar = "A/(2.303RT)"
        else:
            yvar = "log a"
        if "alpha" in diagram_args.keys():
            if diagram_args["alpha"]:
                yvar = "alpha"
        
        df_c = pd.concat(dfs)

        if "fill" in diagram_args.keys():
            if isinstance(diagram_args["fill"], list):
                colormap = {key:col for key,col in zip(list(dict.fromkeys(df_c["variable"])), diagram_args["fill"])}
            else:
                colormap = diagram_args["fill"]

            # with color mapping
            fig = px.line(df_c, x=xvar, y="value", color='variable', template="simple_white",
                          width=500,  height=400, animation_frame=anim_var,
                          color_discrete_map = colormap,
                          labels=dict(value=yvar, x=xvar),
                         )
        else:
            # without color mapping
            fig = px.line(df_c, x=xvar, y="value", color='variable', template="simple_white",
                          width=500,  height=400, animation_frame=anim_var,
                          labels=dict(value=yvar, x=xvar),
                         )
        
        if "annotation" in diagram_args.keys():
            if "annotation_coords" not in diagram_args.keys():
                diagram_args["annotation_coords"] = [0, 0]
            fig.add_annotation(x=diagram_args["annotation_coords"][0],
                               y=diagram_args["annotation_coords"][1],
                               xref="paper",
                               yref="paper",
                               align='left',
                               text=diagram_args["annotation"],
                               bgcolor="rgba(255, 255, 255, 0.5)",
                               showarrow=False)
        
        if 'main' in diagram_args.keys():
            fig.update_layout(title={'text':diagram_args["main"], 'x':0.5, 'xanchor':'center'})

        if isinstance(xlab, str):
            fig.update_layout(xaxis_title=xlab)
        if isinstance(ylab, str):
            fig.update_layout(yaxis_title=ylab)
        
        if 'fill' in diagram_args.keys():
            if isinstance(diagram_args["fill"], list):
                for i,v in enumerate(diagram_args["fill"]):
                    fig['data'][i]['line']['color']=v
        
        fig.update_layout(legend_title=None)

        config = {'displaylogo': False,
                  'modeBarButtonsToRemove': ['resetScale2d', 'toggleSpikelines'],
                  'toImageButtonOptions': {
                                             'format': save_format, # one of png, svg, jpeg, webp
                                             'filename': save_as,
                                             'height': height,
                                             'width': width,
                                             'scale': save_scale,
                                          },
                 }

        fig.show(config=config)
        return
    
    else:
        yvals = __seq(yrange[0], yrange[1], length_out=yres)


    
    frames = []
    slider_steps = []
    annotations = []
    cst_data = []
    heatmaps = []
    
    # i is a frame in the animation
    for i in range(0, len(zvals)):

        annotations_i = []
        for s in sp:
            if s in set(dfs[i]["prednames"]):
                # if an annotation should appear, create one for this frame
                df_s = dfs[i].loc[dfs[i]["prednames"]==s,]
                namex = df_s[xvar].mean()
                namey = df_s[yvar].mean()
                a = go.layout.Annotation(
                    x=namex,
                    y=namey,
                    xref="x",
                    yref="y",
                    text=chemlabel(s),
                    bgcolor="rgba(255, 255, 255, 0.5)",
                    showarrow=False,
                    )
            else:
                # if an annotation shouldn't appear, make an invisible annotation
                # (workaround for a plotly bug where annotations won't clear in an animation)
                namex = statistics.mean(xvals)
                namey = statistics.mean(yvals)
                a = go.layout.Annotation(
                    x=namex,
                    y=namey,
                    xref="x",
                    yref="y",
                    text="",
                    bgcolor="rgba(255, 255, 255, 0)",
                    showarrow=False,
                    )
            annotations_i.append(a)
            
        # allows adding a custom annotation; append to frame
        if "annotation" in diagram_args.keys():
            if "annotation_coords" not in diagram_args.keys():
                diagram_args["annotation_coords"] = [0, 0]
            custom_annotation = go.layout.Annotation(
                    x=diagram_args["annotation_coords"][0],
                    y=diagram_args["annotation_coords"][1],
                    xref="paper",
                    yref="paper",
                    align='left',
                    text=diagram_args["annotation"],
                    bgcolor="rgba(255, 255, 255, 0.5)",
                    showarrow=False,
                    )
            annotations_i.append(custom_annotation)
    
        annotations.append(annotations_i)

        if 'ylab' in diagram_args.keys():
            ylab = diagram_args["ylab"]
            hover_ylab = ylab+': %{y} '
        else:
            ylab = chemlabel(ylab)
            hover_ylab = yvar+': %{y} '+unit_dict[yvar]

        if 'xlab' in diagram_args.keys():
            xlab = diagram_args["xlab"]
            hover_xlab = xlab+': %{x} '
        else:
            xlab = chemlabel(xlab)
            hover_xlab = xvar+': %{x} '+unit_dict[xvar]
        
        heatmaps_i = go.Heatmap(z=dmaps[i], x=xvals, y=yvals, zmin=0, zmax=len(sp)-1,
                                customdata=dmaps_names[i],
                                hovertemplate=hover_xlab+'<br>'+hover_ylab+'<br>Region: %{customdata}<extra></extra>')

        heatmaps.append(heatmaps_i)

        frame = go.Frame(data=[heatmaps_i],
                         name=str(i),
                         layout=go.Layout(annotations=annotations_i))

        frames.append(frame)

        slider_step = dict(
            method='animate',
            label=zvals[i],
            value=i,
            args=[
                [i],
                dict(
                    frame=dict(duration=300, redraw=True),
                    mode='immediate',
                    transition=dict(duration=0)
                )
            ]
        )

        slider_steps.append(slider_step)

    fig = go.Figure(
        data = heatmaps[0],
        layout=go.Layout(
    #         title="Frame 0",
            title_x=0.5,
            width=500, height=500,
            annotations=annotations[0],
            sliders=[dict(
                active=0,
                yanchor='top',
                xanchor='left',
                currentvalue=dict(
                    font=dict(size=12),
                    prefix='{}: '.format(anim_var),
                    suffix=' '+unit_dict[anim_var],
                    visible=True,
                    xanchor='right'
                ),
                transition=dict(duration=0, easing='cubic-in-out'),
                pad=dict(b=10, t=50),
                len=0.9,
                x=0.1,
                y=0,
                steps=slider_steps
            )],
            updatemenus=[dict(
                type="buttons",
                buttons=[dict(label="Play",
                              method="animate",
                              args=[None, {"fromcurrent":True}]),
                        dict(label="Pause",
                             method="animate",
                             args=[[None],
                                   {"frame": {"duration": 0, "redraw": True},
                                    "mode": "immediate",
                                    "transition": {"duration": 0}}],
                             )],
                direction="left",
                pad={"r": 10, "t": 87},
                showactive=False,
                x=0.1,
                xanchor="right",
                y=0,
                yanchor="top",
            )]
        ),
        frames=frames

    )


    if 'fill' in diagram_args.keys():
        if isinstance(diagram_args["fill"], list):
            colorscale_temp = []
            for i,v in enumerate(diagram_args["fill"]):
                colorscale_temp.append([i, v])
            colorscale = colorscale_temp
        elif isinstance(diagram_args["fill"], str):
            colorscale = diagram_args["fill"]
    else:
        colorscale = "viridis"
    
    fig.update_traces(dict(showscale=False,
                           colorscale=colorscale),
                      selector={'type':'heatmap'})
    
    fig.update_layout(
        xaxis_title=xlab,
        yaxis_title=ylab,
        xaxis={"range":[list(dfs[0][xvar])[0], list(dfs[0][xvar])[-1]]},
        yaxis={"range":[list(dfs[0][yvar])[0], list(dfs[0][yvar])[-1]]},
        margin={"t": 60, "r":60},
    )

    if 'main' in diagram_args.keys():
        fig.update_layout(title={'text':diagram_args['main'], 'x':0.5, 'xanchor':'center'})

    config = {'displaylogo': False,
              'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d',
                                         'autoScale2d', 'toggleSpikelines',
                                         'hoverClosestCartesian', 'hoverCompareCartesian'],
              'toImageButtonOptions': {
                                       'format': save_format, # one of png, svg, jpeg, webp
                                       'filename': save_as,
                                       'height': height,
                                       'width': width,
                                       'scale': save_scale,
                                        },
             }
    

    fig.show(config=config)

    
def diagram_interactive(data, title=None, borders=0, names=None, format_names=True,
                        annotation=None, annotation_coords=[0, 0],
                        balance=None, xlab=None, ylab=None, colormap="viridis",
                        width=600, height=520, alpha=False, messages=True,
                        plot_it=True, save_as=None, save_format=None, save_scale=1):
    
    """
    Produce an interactive diagram.
    
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

    format_names : bool, default True
        Apply formatting to chemical formulas?
    
    annotation : str, optional
        Annotation to add to the plot.
    
    annotation_coords : list of numeric, default [0, 0], optional
        Coordinates of annotation, where 0,0 is bottom left and 1,1 is top
        right.
        
    balance : str or numeric, optional
        How to balance the transformations.
    
    xlab, ylab : str
        Custom x and y axes labels.

    colormap : str
        Name of the 
    
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
    """

    if format_names:
        format_x_name = True
        format_y_name = True
    else:
        format_x_name = False
        format_y_name = False
    
    basis_sp = list(data.rx2("basis").rownames)
    basis_state = list(data.rx2("basis").rx2("state"))
    xyvars = list(data.rx2("vars"))
    xyvals = list(data.rx2("vals"))

    if 'loga.equil' not in data.names:
        calc_type = "a"
    else:
        calc_type = "e"

    try:
        balance = float(balance)
    except:
        pass

    if balance is None or balance == "":
        a, args = diagram(data, messages=False, plot_it=False)
        balance = list(a.rx2("n.balance"))

    if calc_type=="a":
        # handling output of affinity()
        out_vals = data.rx2["values"]
        out_units = "A/(2.303RT)"
    else:
        # handling output of equilibrate()
        out_vals = data.rx2("loga.equil")
        out_units = "log a"

    nsp = len(out_vals)

    if isinstance(out_vals[0], rpy2.robjects.vectors.FloatMatrix):
        # if there are two variables in affinity (predominance plot), data.rx2["values"][0] will
        # return a rpy2.robjects.vectors.FloatMatrix representing all species.
        # Turn this floatmatrix into a list of lists, then turn into a dataframe.
        lsp = out_vals[0].nrow*out_vals[0].ncol
        flat_out_vals = np.concatenate(tuple(arr.transpose() for arr in out_vals), axis=0).reshape(nsp, lsp).tolist()
        df = pd.DataFrame(flat_out_vals)
    elif isinstance(out_vals[0], rpy2.robjects.vectors.FloatArray):
        # if there is only a single variable in affinity (activity plot), data.rx2["values"][0] will
        # return a rpy2.robjects.vectors.FloatArray representing only the first species,
        # so loop through all species to create a list of lists, then turn into a dataframe.
        flat_out_vals = []
        for v in out_vals:
            flat_out_vals.append(v)
        df = pd.DataFrame(flat_out_vals)

    if calc_type=="a":
        # handling output of affinity()
        if isinstance(balance, str):
            #df["n.balance"] = list(species().rx2(balance))
            df["n.balance"] = list(species()[balance])
        else:
            df["n.balance"] = balance
        # divide values by balance
        df = df.apply(lambda row: row/row["n.balance"], axis=1)
        df = df.drop(["n.balance"], axis=1)

        sp = list(info([int(val) for val in list(out_vals.names)], messages=False)["name"])
    else:
        # handling output of equilibrate()
        sp = list(info([int(val) for val in list(data.rx2("values").names)], messages=False)["name"])

    if isinstance(names, list) and len(names)==len(sp):
        sp = names

    df.index = sp
    df = df.transpose()

    if alpha and len(xyvars) == 1:
        df = df.apply(lambda x: 10**x)
        df = df[sp].div(df[sp].sum(axis=1), axis=0)

    xvar = xyvars[0]
    xvals = [float(val) for val in xyvals[0]]

    if len(xyvars) == 1:
        df[xvar] = xvals
        if not isinstance(ylab, str):
            if alpha:
                ylab = "alpha"
            else:
                ylab = out_units # A/(2.303RT)
                format_y_name = False
            if format_y_name:
                ylab = chemlabel(ylab)
        df = pd.melt(df, id_vars=xyvars, value_vars=sp)

    elif len(xyvars)==2:
        # predominance plot
        yvar = xyvars[1]
        yvals = [float(val) for val in xyvals[1]]

        df["pred"] = df.idxmax(axis=1,skipna=True)
        df["prednames"] = df["pred"]

        xvals_full = xvals*len(yvals)

        yvals_full = __flatten_list([[y]*len(xvals) for y in yvals])
        df[xvar] = xvals_full
        df[yvar] = yvals_full

    unit_dict = {"P":"bar", "T":"°C", "pH":"", "Eh":"volts", "IS":"mol/kg"}

    for i,s in enumerate(basis_sp):
        if basis_state[i] in ["aq", "liq", "cr"]:
            if format_names:
                unit_dict[s] = "log <i>a</i><sub>{}</sub>".format(chemlabel(s))
            else:
                unit_dict[s] = "log <i>a</i><sub>{}</sub>".format(s)
        else:
            if format_names:
                unit_dict[s] = "log <i>f</i><sub>{}</sub>".format(chemlabel(s))
            else:
                unit_dict[s] = "log <i>f</i><sub>{}</sub>".format(s)

    if not isinstance(xlab, str):
        xlab = xvar+", "+unit_dict[xvar]
        if xvar == "pH":
            xlab = "pH"
        if xvar in basis_sp:
            xlab = unit_dict[xvar]
        if format_x_name:
            xlab = chemlabel(xlab)

    if len(xyvars) == 1:

        if format_names:
            df["variable"] = df["variable"].apply(chemlabel)

        fig = px.line(df, x=xvar, y="value", color='variable', template="simple_white",
                      width=width,  height=height,
                      labels=dict(value=ylab, x=xlab), render_mode='svg',
                     )

        if isinstance(colormap, list):
            for i,v in enumerate(colormap):
                fig['data'][i]['line']['color']=v

        fig.update_layout(xaxis_title=xlab,
                          yaxis_title=ylab,
                          legend_title=None)

        if isinstance(title, str):
            fig.update_layout(title={'text':title, 'x':0.5, 'xanchor':'center'})

        if isinstance(annotation, str):
            fig.add_annotation(
                x=annotation_coords[0],
                y=annotation_coords[1],
                text=annotation,
                showarrow=False,
                xref="paper",
                yref="paper",
                align='left',
                bgcolor="rgba(255, 255, 255, 0.5)")

        save_as, save_format = __save_figure(fig, save_as, save_format, save_scale,
                                             plot_width=width, plot_height=height, ppi=1)

        config = {'displaylogo': False,
                  'modeBarButtonsToRemove': ['resetScale2d', 'toggleSpikelines'],
                  'toImageButtonOptions': {
                                             'format': save_format, # one of png, svg, jpeg, webp
                                             'filename': save_as,
                                             'height': height,
                                             'width': width,
                                             'scale': save_scale,
                                          },
                 }

    if len(xyvars) == 2:
        mappings = {'pred': {s:lab for s,lab in zip(sp,range(0,len(sp)))}}
        with pd.option_context('future.no_silent_downcasting', True):
            df.replace(mappings, inplace=True)

        data = np.array(df.pred)
        shape = (len(xvals), len(yvals))
        dmap = data.reshape(shape)

        data = np.array(df.prednames)
        shape = (len(xvals), len(yvals))
        dmap_names = data.reshape(shape)

        if not isinstance(ylab, str):
            ylab = yvar+", "+unit_dict[yvar]
            if yvar in basis_sp:
                ylab = unit_dict[yvar]

            if yvar == "pH":
                ylab = "pH"

            if format_y_name:
                ylab = chemlabel(ylab)

        fig = px.imshow(dmap, width=width, height=height, aspect="auto",
                        labels=dict(x=xlab, y=ylab, color="region"),
                        x=xvals, y=yvals, template="simple_white",
                       )

        fig.update(data=[{'customdata': dmap_names,
            'hovertemplate': xlab+': %{x}<br>'+ylab+': %{y}<br>Region: %{customdata}<extra></extra>'}])

        if colormap == 'none':
            colormap = [[0, 'white'], [1, 'white']]
        elif isinstance(colormap, list):
            colmap_temp = []
            for i,v in enumerate(colormap):
                colmap_temp.append([i,v])
            colormap = colmap_temp

        fig.update_traces(dict(showscale=False, 
                               coloraxis=None, 
                               colorscale=colormap),
                          selector={'type':'heatmap'})

        fig.update_yaxes(autorange=True)

        if isinstance(title, str):
            fig.update_layout(title={'text':title, 'x':0.5, 'xanchor':'center'})

        for s in sp:
            if s in set(df["prednames"]):
                df_s = df.loc[df["prednames"]==s,]
                namex = df_s[xvar].mean()
                namey = df_s[yvar].mean()

                if format_names:
                    annot_text = chemlabel(s)
                else:
                    annot_text = str(s)

                fig.add_annotation(x=namex, y=namey,
                                   text=annot_text,
                                   bgcolor="rgba(255, 255, 255, 0.5)",
                                   showarrow=False)

        if isinstance(annotation, str):
            fig.add_annotation(
                x=annotation_coords[0],
                y=annotation_coords[1],
                text=annotation,
                showarrow=False,
                xref="paper",
                yref="paper",
                align='left',
                bgcolor="rgba(255, 255, 255, 0.5)")

        if borders > 0:

            unique_x_vals = list(dict.fromkeys(df[xvar]))
            unique_y_vals = list(dict.fromkeys(df[yvar]))

            def mov_mean(numbers=[], window_size=2):
                i = 0
                moving_averages = []
                while i < len(numbers) - window_size + 1:
                    this_window = numbers[i : i + window_size]

                    window_average = sum(this_window) / window_size
                    moving_averages.append(window_average)
                    i += 1
                return moving_averages

            x_mov_mean = mov_mean(unique_x_vals)
            y_mov_mean = mov_mean(unique_y_vals)

            x_plot_min = x_mov_mean[0] - (x_mov_mean[1] - x_mov_mean[0])
            y_plot_min = y_mov_mean[0] - (y_mov_mean[1] - y_mov_mean[0])

            x_plot_max = x_mov_mean[-1] + (x_mov_mean[1] - x_mov_mean[0])
            y_plot_max = y_mov_mean[-1] + (y_mov_mean[1] - y_mov_mean[0])

            x_vals_border = [x_plot_min] + x_mov_mean + [x_plot_max]
            y_vals_border = [y_plot_min] + y_mov_mean + [y_plot_max]

            data = np.array(df.pred)
            shape = (len(xvals), len(yvals))
            dmap = data.reshape(shape)

            def find_line(dmap, row_index):
                return [i for i in range(0, len(dmap[row_index])-1) if dmap[row_index][i] != dmap[row_index][i+1]]

            nrows, ncols = dmap.shape
            vlines = []
            for row_i in range(0, nrows):
                vlines.append(find_line(dmap, row_i))

            dmap_transposed = dmap.transpose()
            nrows, ncols = dmap_transposed.shape
            hlines = []
            for row_i in range(0, nrows):
                hlines.append(find_line(dmap_transposed, row_i))
            y_coord_list_vertical = []
            x_coord_list_vertical = []
            for i,row in enumerate(vlines):
                for line in row:
                    x_coord_list_vertical += [x_vals_border[line+1], x_vals_border[line+1], np.nan]
                    y_coord_list_vertical += [y_vals_border[i], y_vals_border[i+1], np.nan]

            y_coord_list_horizontal = []
            x_coord_list_horizontal = []
            for i,col in enumerate(hlines):
                for line in col:
                    y_coord_list_horizontal += [y_vals_border[line+1], y_vals_border[line+1], np.nan]
                    x_coord_list_horizontal += [x_vals_border[i], x_vals_border[i+1], np.nan]

            fig.add_trace(
                go.Scatter(
              mode= 'lines',
              x= x_coord_list_horizontal,
              y= y_coord_list_horizontal,
              line= {
                "width": borders,
                "color": 'black'
              },
              hoverinfo= 'skip',
              showlegend=False,
                )
            )
            fig.add_trace(
                go.Scatter(
              mode= 'lines',
              x= x_coord_list_vertical,
              y= y_coord_list_vertical,
              line= {
                "width": borders,
                "color": 'black'
              },
              hoverinfo= 'skip',
              showlegend=False,
                )
            )

            fig.update_yaxes(range=[min(yvals), max(yvals)], autorange=False, mirror=True)
            fig.update_xaxes(range=[min(xvals), max(xvals)], autorange=False, mirror=True)


    save_as, save_format = __save_figure(fig, save_as, save_format, save_scale,
                                         plot_width=width, plot_height=height, ppi=1)

    config = {'displaylogo': False,
              'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d',
                                         'autoScale2d', 'resetScale2d', 'toggleSpikelines',
                                         'hoverClosestCartesian', 'hoverCompareCartesian'],
              'toImageButtonOptions': {
                                         'format': save_format, # one of png, svg, jpeg, webp
                                         'filename': save_as,
                                         'height': height,
                                         'width': width,
                                         'scale': save_scale,
                                      },
             }

    if plot_it:
        fig.show(config=config)

    return df, fig


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
    elif all(isinstance(x, (int, np.integer)) for x in value):
        return ro.IntVector(value)
    elif all(isinstance(x, (int, np.integer, float)) for x in value):
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
    
    args = rpy2float(args)
    
    capture = R_output()
    capture.capture_r_output()
    
    out = CHNOSZ.water(**args)

    if messages:
        capture.print_captured_r_output()
    
    if property not in ["SUPCRT92", "SUPCRT", "IAPWS95", "IAPWS", "DEW"]:

        return r_df_to_pd(out)
        
        
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
    
    capture = R_output()
    capture.capture_r_output()
        
    out = CHNOSZ.entropy(formula_R)

    if messages:
        capture.print_captured_r_output()
    
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
    
    capture = R_output()
    capture.capture_r_output()
    
    out = CHNOSZ.mass(formula_R)

    if messages:
        capture.print_captured_r_output()
    
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
    
    capture = R_output()
    capture.capture_r_output()
    
    out = CHNOSZ.ZC(formula_R)

    if messages:
        capture.print_captured_r_output()
    
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
    
    capture = R_output()
    capture.capture_r_output()
        
    out = CHNOSZ.makeup(**args)

    return out
    
    if messages:
        capture.print_captured_r_output()
    
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

    capture = R_output()
    capture.capture_r_output()
    
    pout = CHNOSZ.seq2aa(**args)

    if messages:
        capture.print_captured_r_output()
    
    return r_df_to_pd(pout)


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
    aa = pd_to_r_df(aa)
    args = {'aa':aa}

    capture = R_output()
    capture.capture_r_output()
    
    apout = CHNOSZ.add_protein(**args)

    if messages:
        capture.print_captured_r_output()
    
    return list(apout)


def equilibrate(aout, balance=None, loga_balance=None, ispecies=None,
                normalize=False, messages=True, dash=False):
    
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

    dash : bool, default False
        Is this function being used in a Dash web application? Required for rpy2
        package compatibility.
    
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

    capture = R_output()
    capture.capture_r_output()

    eout = CHNOSZ.equilibrate(**args)
        
    if messages:
        capture.print_captured_r_output()
    
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
            annotation=None, annotation_coords=[0,0],
            width=600, height=520, dpi=150,
            messages=True, interactive=False,
            save_as=None, save_format=None,
            save_scale=1, fig_out=False):
    
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
    """
    
    if lwd == None:
        borders = 0
    else:
        borders = lwd
    
    if interactive:
            
        df, fig = diagram_interactive(
                                 data=eout,
                                 title=main,
                                 borders=borders,
                                 names=names,
                                 format_names=format_names,
                                 annotation=annotation,
                                 annotation_coords=annotation_coords,
                                 balance=balance,
                                 xlab=xlab, ylab=ylab,
                                 colormap=fill,
                                 width=width, height=height,
                                 alpha=alpha, plot_it=plot_it,
                                 save_as=save_as, save_format=save_format,
                                 save_scale=save_scale, messages=messages)
        if fig_out:
            return df, fig
        else:
            return df
    
    
    args = {'eout':eout, 'type':ptype, 'alpha':alpha, 'normalize':normalize,
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
    
    capture = R_output()
    capture.capture_r_output()
        
    with __r_inline_plot(width=width, height=height, dpi=dpi, plot_it=plot_it):
        if isinstance(add, bool):
            if add: # add='True' does not work with the current pyCHNOSZ framework
                raise Exception("The argument 'add' must be assigned the output of the previous diagram(s).")
            else:
                d = CHNOSZ.diagram(**args)
        elif isinstance(add, list):
            args.update({'add':True})
            for to_add in add:
                CHNOSZ.diagram(**to_add)
            d = CHNOSZ.diagram(**args)

        else:
            CHNOSZ.diagram(**add)
            args.update({'add':True})
            d = CHNOSZ.diagram(**args)
            
    if messages:
        capture.print_captured_r_output()
    
    return d, args


def affinity(property=None, sout=None, exceed_Ttr=False,
             exceed_rhomin=False, return_buffer=False, return_sout=False,
             balance="PBB", iprotein=None, loga_protein=-3, transect=None,
             messages=True, dash=False, **kwargs):

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

    dash : bool, default False
        Is this function being used in a Dash web application? Required for
        rpy2 compatibility.
        
    **kwargs : numeric, zero or more named arguments
        Used to identify the variables of interest in the calculations. E.g.,
        pH, T, P, Eh, ...
    
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
        
    capture = R_output()
    capture.capture_r_output()

    a = CHNOSZ.affinity(**args)

    if messages:
        capture.print_captured_r_output()

    return a

            
def species(species=None, state=None, delete=False, add=False,
            index_return=False, messages=True, dash=False):
    
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

    dash : bool, default False
        Is this function being used in a Dash web application? Required for rpy2
        compatibility.
        
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
    
    capture = R_output()
    capture.capture_r_output()

    sout = CHNOSZ.species(**args)
    sout = r_df_to_pd(sout)
        
    if messages:
        capture.print_captured_r_output()
    
    return sout


def basis(species=None, state=None, logact=None, delete=False,
          messages=True, dash=False):
    
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

    dash : bool, default False
        Is this function being used in a Dash web application? Required for rpy2
        compatibility.
        
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
    
    capture = R_output()
    capture.capture_r_output()

    bout = CHNOSZ.basis(**args)
    bout = r_df_to_pd(bout)
        
    if messages:
        capture.print_captured_r_output()
    
    return bout


def reset(db="OBIGT", messages=True):
    
    """
    Python wrapper for the reset() function in CHNOSZ.
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
    
    """
    
    if db == "WORM":
        url = "https://raw.githubusercontent.com/worm-portal/WORM-db/master/wrm_data.csv"
        url_ref = "https://raw.githubusercontent.com/worm-portal/WORM-db/master/references.csv"

        
        
        # Load WORM database by detault
        if can_connect_to(url):
            # Download from URL and decode as UTF-8 text.
            with urlopen(url) as webpage:
                content = webpage.read().decode()
            WORM_DB = pd.read_csv(StringIO(content), sep=",")

            obigt_r = thermo()["OBIGT"]
            
            obigt_pd = r_df_to_pd(obigt_r)
            obigt_pd = _assign_OBIGT_pandas_dtypes(obigt_pd)
            
            thermo_trunc = obigt_pd[obigt_pd["name"].isin(["water", "H+", "e-"])]
            thermo(**{"OBIGT":thermo_trunc})
            
            if can_connect_to(url_ref):
                with urlopen(url_ref) as webpage:
                    content = webpage.read().decode()
                WORM_refs = pd.read_csv(StringIO(content), sep=",")
                thermo(**{"refs":WORM_refs})

            _ = add_OBIGT(WORM_DB, messages=False)

            if messages:
                naq = WORM_DB[WORM_DB["state"] == "aq"].shape[0]
                ntot = WORM_DB.shape[0]
                print("The WORM thermodynamic database has been loaded: {0} aqueous, {1} total species".format(str(naq), str(ntot)))
        else:
            db = "OBIGT"
            if messages:
                print("Could not reach the WORM database repository. The OBIGT database has been loaded instead.")
            
    if db == "OBIGT":
        capture = R_output()
        capture.capture_r_output()

        CHNOSZ.reset()

        if messages:
            capture.print_captured_r_output()


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

    add_obigt_builtins = ["SUPCRT92", "AS04", "DEW", "OldAA",
                          "AkDi", "SLOP98", "Berman_cr", "biotic_aq",
                          "H2O_aq", "inorganic_aq", "inorganic_cr",
                          "inorganic_gas", "organic_aq", "organic_cr",
                          "organic_gas", "organic_liq"]

    if isinstance(file, str):
        
        if file in add_obigt_builtins:
            capture = R_output()
            capture.capture_r_output()

            args={'file':file}

            if species != None:
                if not isinstance(species, list):
                    args["species"] = species
                else:
                    args["species"] = _convert_to_RVector(species)

            ispecies = CHNOSZ.add_OBIGT(**args)

            if messages:
                capture.print_captured_r_output()

            return list(ispecies)
        
        elif ".csv" not in file[-4:] or ".CSV" in file[-4:]:
            raise Exception("File must be in .csv format")
        else:
            df = pd.read_csv(file, keep_default_na=False, na_values=['NA']) # keep_default_na=False keeps NA in data file instead of converting them to NaN. NaNs cause errors when converting to an R dataframe with rpy2 3.4.5.
            if df.shape[0] == 0:
                raise Exception("file is empty")
    elif isinstance(file, pd.DataFrame):
        df = file
        if df.shape[0] == 0:
            raise Exception("file is empty")
    else:
        raise Exception("file parameter must be either the name of a .csv file "
                        "to import, or a Pandas dataframe of thermodynamic data.")

    OBIGT_cols = ['name', 'abbrv', 'formula', 'state',
                  'ref1', 'ref2', 'date', 'E_units',
                  'G', 'H', 'S', 'Cp', 'V',
                  'a1.a', 'a2.b', 'a3.c', 'a4.d',
                  'c1.e', 'c2.f', 'omega.lambda', 'z.T']

    if all(col in df.columns for col in OBIGT_cols):
        
        df_mod_OBIGT = copy.deepcopy(df[OBIGT_cols])

        if "formula_ox" in df.columns:
            thermo(formula_ox=copy.deepcopy(df[["name", "formula_ox"]]))
        
        if species != None:
            if isinstance(species, list):
                df_mod_OBIGT = df_mod_OBIGT[df_mod_OBIGT["name"].isin(species)]
                
            else:
                raise Exception("species must be a list of names of species to load from file.")

        if not force:
            t = thermo()
            df_mod_OBIGT = df_mod_OBIGT[~df_mod_OBIGT["name"].isin(t.OBIGT.name)]
            if df_mod_OBIGT.shape[0] == 0:
                raise Exception("No species to add while force=False")
                
    else:
        msg = "The file must contain all of the columns: {}".format(str(OBIGT_cols))
        raise Exception(msg)
    
    df_mod_OBIGT = _assign_OBIGT_pandas_dtypes(df_mod_OBIGT)
    
    return mod_OBIGT(df_mod_OBIGT, messages=messages)


def _assign_OBIGT_pandas_dtypes(pd_df):
    str_cols = ['name', 'abbrv', 'formula', 'state', 'ref1',
                'ref2', 'date', 'E_units']
    
    numeric_cols = ['G', 'H', 'S', 'Cp', 'V',
                  'a1.a', 'a2.b', 'a3.c', 'a4.d',
                  'c1.e', 'c2.f', 'omega.lambda', 'z.T']
    
    pd_df[str_cols] = pd_df[str_cols].astype(str)
    pd_df[numeric_cols] = pd_df[numeric_cols].apply(pd.to_numeric)
    
    return pd_df


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
    capture = R_output()
    capture.capture_r_output()
        
    if isinstance(args[0], pd.DataFrame):
        
        for col in args[0].columns:
            if sum(args[0][col].isna()) > 0 and args[0][col].dtype == "O":
                if messages:
                    print("mod_OBIGT Warning: The column " + str(col) + " "
                          "has dtype 'object' but contains NA values that "
                          "can be interpreted as 'float'. Removing NA values.")
                args[0][col] = args[0][col].fillna('')
        
        
        arg_list = list(args)
        arg_list[0] = pd_to_r_df(arg_list[0])
        args = tuple(arg_list)
    else:
        pass
    
    ispecies = CHNOSZ.mod_OBIGT(*args, **kwargs)
    
    if messages:
        capture.print_captured_r_output()
    
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
    if all(isinstance(x, (int, np.integer)) for x in species):
        output_is_df = True
    
    if state != None:
        args["state"] = _convert_to_RVector(state, force_Rvec=False)
    
    args["check.it"] = check_it

    capture = R_output()
    capture.capture_r_output()
    
    a = CHNOSZ.info(**args)
    
    if messages:
        capture.print_captured_r_output()
    
    if output_is_df:
        return r_df_to_pd(a)
    else:
        return list(a)

    
def subcrt(species, coeff=None, state=None,
           property=["logK", "G", "H", "S", "V", "Cp"],
           T=None, P=None, grid=None,
           convert=True, exceed_Ttr=True, exceed_rhomin=False,
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

    exceed_Ttr : bool, default True
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
    
    capture = R_output()
    capture.capture_r_output()
        
    a = CHNOSZ.subcrt(**args)

    if messages:
        capture.print_captured_r_output()
    
    if "warnings" in a.names:
        warn = a.rx2("warnings")[0] # subcrt's list includes warnings only if they appear
    else:
        warn = None
    
    if "polymorphs" in a.names:
        poly = a.rx2("polymorphs") # subcrt's list includes polymorphs only if they appear
    else:
        poly = None
    
    if not single_species:
        out_dict = {"reaction":r_df_to_pd(a[0]),
                    "out":r_df_to_pd(a[1])} # the extra [0] is important
    else:
        out_dict = {"species":r_df_to_pd(a[0]), "out":{}}
        
        i=0
        for df in a[1]:
            out_dict["out"][out_dict["species"].name[i]] = r_df_to_pd(df)
            i += 1
        
    if isinstance(warn, str):
        out_dict["warnings"] = warn
    if isinstance(poly, pd.DataFrame):
        out_dict["polymorphs"] = poly
    
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


class thermo:
    
    """
    Python wrapper for the thermo() object in CHNOSZ.
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
        
    """
    
    def __init__(self, db="OBIGT", messages=True, **kwargs):
        
        if isinstance(db, str):
            if db == "WORM":
                reset("WORM", messages=messages)
        
        args = {}
        for key, value in kwargs.items():
            if isinstance(value, list) or isinstance(value, str) or isinstance(value, NumberTypes):
                value = _convert_to_RVector(value, force_Rvec=False)
            elif isinstance(value, pd.DataFrame):
                value = pd_to_r_df(value)
            else:
                pass
            args.update({key:value})
        
        capture = R_output()
        capture.capture_r_output()
            
        t = CHNOSZ.thermo(**args)

        if messages:
            capture.print_captured_r_output()
        
        for i, name in enumerate(t.names):
            if isinstance(t[i], ro.DataFrame) or isinstance(t[i], ro.Matrix):
                attr = r_df_to_pd(t[i])
            elif isinstance(t[i], ro.ListVector):
                attr = {}
                for ii, subname in enumerate(t[i].names):
                    attr[subname] = list(t[i][ii])
            elif isinstance(t[i], ro.FloatVector) or isinstance(t[i], ro.IntVector) or isinstance(t[i], ro.StrVector):
                attr = list(t[i])
                if len(t[i]) == 1:
                    attr = attr[0]
            elif isinstance(t[i], rpy2.rinterface_lib.sexp.NULLType):
                attr = None
            else:
                attr = t[i]
            
            
            setattr(self, name, attr)
            
    def __getitem__(self, item):
         return getattr(self, item)

        
def unicurve(logK, species, coeff, state, pressures=1, temperatures=25, IS=0,
             minT=0.1, maxT=100, minP=1, maxP=500, tol=None,
             solve="T", width=600, height=520, dpi=90, plot_it=True,
             messages=True, show=True):
    
    """
    Solve for temperatures or pressures of equilibration for a given logK
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
    
    """

    if tol == None:
        d = decimal.Decimal(str(logK))
        n_decimals = abs(d.as_tuple().exponent)
        tol = float("0."+"".join(n_decimals*["0"])+"01")
        if tol > 0.00001:
            tol = 0.00001
                
    species = _convert_to_RVector(species, force_Rvec=False)
    state = _convert_to_RVector(state, force_Rvec=False)
    coeff = _convert_to_RVector(coeff, force_Rvec=False)
    pressures = _convert_to_RVector(pressures, force_Rvec=False)
    temperatures = _convert_to_RVector(temperatures, force_Rvec=False)

    capture = R_output()
    capture.capture_r_output()

    r_univariant = import_package_file(__name__, "univariant.r")

    ro.r(r_univariant)
    
    if solve=="T":
        a = ro.r.uc_solveT(logK=logK,
                           species=species,
                           coeff=coeff,
                           state=state,
                           pressures=pressures,
                           IS=IS,
                           minT=minT,
                           maxT=maxT,
                           tol=tol)
        if plot_it:
            with __r_inline_plot(width=width, height=height, dpi=dpi, plot_it=plot_it):
                ro.r.create_output_plot_T(logK=logK,
                                          species=species,
                                          coeff=coeff,
                                          state=state,
                                          pressures=pressures,
                                          minT=minT,
                                          maxT=maxT)
    elif solve=="P":
        a = ro.r.uc_solveP(logK=logK,
                           species=species,
                           state=state,
                           coeff=coeff,
                           temperatures=temperatures,
                           IS=IS,
                           minP=minP,
                           maxP=maxP,
                           tol=tol)
        if plot_it:
            with __r_inline_plot(width=width, height=height, dpi=dpi, plot_it=plot_it):
                ro.r.create_output_plot_P(logK=logK,
                                          species=species,
                                          state=state,
                                          coeff=coeff,
                                          temperatures=temperatures,
                                          minP=minP,
                                          maxP=maxP)
    if messages:
        capture.print_captured_r_output()

    if len(a) == 3:
        warn = a[2][0] # subcrt's list includes warnings only if they appear
    else:
        warn = None
    
    out_dict = {"reaction":r_df_to_pd(a[0]),
                "out":r_df_to_pd(a[1])} # the extra [0] is important
        
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


def univariant_TP(logK, species, coeff, state, Trange, Prange, IS=0,
                  xlim=None, ylim=None, line_type="markers+lines",
                  tol=None, title=None, res=10, width=500, height=400,
                  save_as=None, save_format=None, save_scale=1,
                  show=False, messages=False, plot_it=True):

    """
    Solve for temperatures or pressures of equilibration for a given logK
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
    
    """
    
    if not isinstance(logK, list):
        logK = [logK]
    
    fig = go.Figure()
    
    output = []
    
    for this_logK in logK:
        
        if tol == None:
            d = decimal.Decimal(str(this_logK))
            n_decimals = abs(d.as_tuple().exponent)
            tol = float("0."+"".join(n_decimals*["0"])+"01")
            if tol > 0.00001:
                tol = 0.00001
        
        out = unicurve(solve="T",
                       logK=this_logK,
                       species=species,
                       state=state,
                       coeff=coeff,
                       pressures=list(np.linspace(start=Prange[0], stop=Prange[1], num=res)),
                       minT=Trange[0],
                       maxT=Trange[1],
                       IS=IS,
                       tol=tol,
                       show=show,
                       plot_it=False, messages=messages)
        
        if not out["out"]["T"].isnull().all():
            fig.add_trace(go.Scatter(
                x=out["out"]["T"],
                y=out["out"]["P"],
                mode=line_type,
                name="logK="+str(this_logK),
                text = ["logK="+str(this_logK) for i in range(0, len(out["out"]["T"]))],
                hovertemplate = '%{text}<br>T, °C=%{x:.2f}<br>P, bar=%{y:.2f}<extra></extra>',
            ))
        else:
            print("Could not find any T or P values in this range that correspond to a logK value of {}".format(this_logK))
        output.append(out)
    
    if title == None:
        react_grid = output[0]["reaction"]
        react_grid["name"] = [name  if name != "water" else "H2O" for name in react_grid["name"]] # replace any "water" with "H2O" in the written reaction
        
        reactants = " + ".join([(str(-int(react_grid["coeff"].iloc[i]) if isinstance(react_grid["coeff"].iloc[i], (int, np.integer)) else -react_grid["coeff"].iloc[i])+" " if -react_grid["coeff"].iloc[i] != 1 else "") + chemlabel(react_grid["name"].iloc[i]) for i in range(0, len(react_grid["name"])) if react_grid["coeff"].iloc[i] < 0])
        products = " + ".join([(str(int(react_grid["coeff"].iloc[i]) if isinstance(react_grid["coeff"].iloc[i], (int, np.integer)) else react_grid["coeff"].iloc[i])+" " if react_grid["coeff"].iloc[i] != 1 else "") + chemlabel(react_grid["name"].iloc[i]) for i in range(0, len(react_grid["name"])) if react_grid["coeff"].iloc[i] > 0])
        
        title = reactants + " = " + products
    
    fig.update_layout(template="simple_white",
                      title=str(title),
                      xaxis_title="T, °C",
                      yaxis_title="P, bar",
                      width=width,
                      height=height,
                      hoverlabel=dict(bgcolor="white"),
    )
    
    if isinstance(xlim, list):
        fig.update_xaxes(range = xlim)
    if isinstance(ylim, list):
        fig.update_yaxes(range = ylim)
    
    save_as, save_format = __save_figure(fig, save_as, save_format, save_scale,
                                         plot_width=width, plot_height=height, ppi=1)
    
    config = {'displaylogo': False,
              'modeBarButtonsToRemove': ['resetScale2d', 'toggleSpikelines'],
              'toImageButtonOptions': {
                                         'format': save_format, # one of png, svg, jpeg, webp
                                         'filename': save_as,
                                         'height': height,
                                         'width': width,
                                         'scale': save_scale,
                                      },
             }
    
    if plot_it:
        fig.show(config=config)
    
    return output


def syslab(system=["K2O", "Al2O3", "SiO2", "H2O"], dash="-"):

    """
    Python wrapper for the syslab() function in CHNOSZ.
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
    """
    
    return dash.join([chemlabel(sp)for sp in system])


def ratlab(top="K+", bottom="H+", molality=False):
    
    """
    Python wrapper for the ratlab() function in CHNOSZ.
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
    """
    
    top_formula = chemparse.parse_formula(top)
    if "+" in top_formula.keys():
        top_charge = top_formula["+"]
    elif "-" in top_formula.keys():
        top_charge = -top_formula["-"]
    else:
        raise Exception("Cannot create an ion ratio involving one or more neutral species.")
    
    bottom_formula = chemparse.parse_formula(bottom)
    if "+" in bottom_formula.keys():
        bottom_charge = bottom_formula["+"]
    elif "-" in bottom_formula.keys():
        bottom_charge = -bottom_formula["-"]
    else:
        raise Exception("Cannot create an ion ratio involving one or more neutral species.")
    
    if top_charge.is_integer():
        top_charge = int(top_charge)

    if bottom_charge.is_integer():
        bottom_charge = int(bottom_charge)
    
    if top_charge != 1:
        top_charge = "<sup>"+str(top_charge)+"</sup>"
    else:
        top_charge = ""
    if bottom_charge != 1:
        bottom_charge = "<sup>"+str(bottom_charge)+"</sup>"
    else:
        bottom_charge = ""
    
    if molality:
        sym = "m"
    else:
        sym = "a"
        
    return "log("+sym+bottom_charge+"<sub>"+chemlabel(top)+"</sub>/"+sym+top_charge+"<sub>"+chemlabel(bottom)+"</sub>)"


def get_formula_ox(name):
    """
    Get quantities of elements and their oxidation states in a chemical compound
    of interest. This function only works when the WORM thermodynamic database
    is loaded. For example, an input of "magnetite" would return the following:
    `{'Fe+3': 2.0, 'Fe+2': 1.0, 'O-2': 4.0}`.

    Parameters
    ----------
    name : str or int
        The name or database index of the chemical species of interest. Example:
        `"magnetite"` or `738`.

    Returns
    -------
    out : dict
        A dictionary where each key represents an element in a specific
        oxidation state, and its value is the number of that element in the
        chemical species' formula.
    
    """
    
    if not isinstance(name, str) and not isinstance(name, int):
        raise Exception("Must provide input as a string (chemical species name) or an integer (chemical species index).")

    # convert ispecies to name
    if isinstance(name, int):
        name = info(name, messages=False).name.iloc[0]
    
    if "formula_ox" in dir(thermo()):
        df = thermo().formula_ox

        if name not in list(df["name"]):
            raise Exception("The species " + str(name) + " was not found in the loaded thermodynamic database.")

        try:
            df[df["name"]==name]["formula_ox"].iloc[0].split()
        except:
            raise Exception("The species " + str(name) + " does not have "
                    "elemental oxidation states given in the 'formula_ox' "
                    "column of the loaded thermodynamic database.")
        
        split_list = df[df["name"]==name]["formula_ox"].iloc[0].split()
        split_list_clean = [s.replace(" ", "") for s in split_list]

        try:
            elem_ox_names = [re.findall(r"^(?:\d+|)([A-Z].*$)", s)[0] for s in split_list_clean]
        except:
            elem_ox_names = []

        elem_ox_list = []
        for s in split_list:
            coeff = re.findall(r"(\d+)[A-Z]", s)
            if len(coeff) == 0:
                coeff = 1
            else:
                coeff = float(coeff[0])
            elem_ox_list.append(coeff)

        return {key:val for key,val in zip(elem_ox_names, elem_ox_list)}


def get_n_element_ox(names, element_ox, binary=False):
    """
    Get the number of an element of a chosen oxidation state in the formula of a
    list of chemical species. This function only works when the WORM
    thermodynamic database is loaded. Example:
    
    If binary is False, returns a list containing the number of the chosen
    element and oxidation state in the chemical species. For example, how many
    ferrous irons are in the formulae of hematite, fayalite, and magnetite,
    respectively?
    ```
    get_n_element_ox(names=["hematite", "fayalite", "magnetite"],
                     element_ox="Fe+2",
                     binary=False)
    ```
    will return the following list representing the number of ferrous irons
    in the formulas of hematite, fayalite, and magnetite, respectively:
    ```
    [0, 2.0, 1.0]
    ```
    If binary is True, returns a list of whether or not ferrous iron is in their
    formulas:
    ```
    [False, True, True]
    ```

    Parameters
    ----------
    names : str or int, or a list of str or int
        The name or database index of a chemical species, or a list of
        names or indices. Example: `["hematite", "fayalite", "magnetite"]` or
        `[788, 782, 798]`.

    element_ox : str
        An element with a specific oxidation state. For example: `"Fe+2"` for
        ferrous iron.

    binary : bool, default False
        Should the output list show True/False for presence or absence of the
        element defined by `element_ox`? By default, this parameter is set to
        False so the output list shows quantities of the element instead.

    Returns
    -------
    out_list : list of float
        A list containing quantities of the chosen element oxidation state in
        the formulas of the chemical species (if `binary=False`) or whether the
        chosen element oxidation state is present in the formulae (if `binary=
        True`).
    
    """

    if not isinstance(names, list):
        names = list(names)
    
    n_list = []
    for name in names:
        d = get_formula_ox(name)
        n_list.append(d.get(element_ox, 0))

    if binary:
        out_list = [True if n!=0 else False for n in n_list]
    else:
        out_list = n_list
    
    return out_list


def water_lines(d, line_color=None, line_type=None, line_width=1):
    messages = False
    add_saturation_lines(a, d)