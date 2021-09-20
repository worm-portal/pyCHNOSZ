import os
import warnings
from contextlib import contextmanager
from IPython.display import Image, display
import pandas as pd
import numpy as np
import re
import plotly.express as px
import plotly.graph_objects as go
import copy
import statistics
import pkg_resources

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
def __r_inline_plot(width=600, height=520, dpi=150, plot_it=True):
    
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


def html_chemname_format(name):
    p = re.compile(r'(?P<sp>\+|-\d+?$)')
    name = p.sub(r'<sup>\g<sp></sup>', name)
    
    name_no_charge = re.match(r'(?:(?!<|$).)*', name).group(0)
    mapping = {"0": "<sub>0</sub>", "1": "<sub>1</sub>", "2": "<sub>2</sub>", "3": "<sub>3</sub>", "4": "<sub>4</sub>", 
           "5": "<sub>5</sub>", "6": "<sub>6</sub>", "7": "<sub>7</sub>", "8": "<sub>8</sub>", "9": "<sub>9</sub>",
           ".":"<sub>.</sub>"}
    name_no_charge_formatted = "".join([mapping.get(x) or x for x in list(name_no_charge)])
    name = re.sub(name_no_charge, name_no_charge_formatted, name)

    return(name)


def __seq(start, end, by=None, length_out=None):
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


def animation(basis_args={}, species_args={}, affinity_args={},
              equilibrate_args=None, diagram_args={},
              anim_var="T", anim_range=[0, 350, 8],
              messages=False):
    
    """
    Produce an animated interactive affinity, activity, or predominance diagram.
    
    Parameters
    ----------
    basis_args : dict
        Dictionary of options for defining basis species (see `basis`) in the
        animated diagram.
        Example: basis_args={'species':['CO2', 'O2', 'H2O', 'H+']}

    basis_args : dict
        Dictionary of options for defining species (see `species`) in the
        animated diagram.
        Example: species_args={'species':['CO2', 'HCO3-', 'CO3-2']}

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
            raise Exception("basis_args needs a list of basis species for 'species'. "
                            "Example: basis_args={'species':['CO2', 'O2', 'H2O', 'H+']}")
    else:
        raise Exception("basis_args needs to be a Python dictionary with a key "
                        "called 'species' (additional keys are optional). "
                        "Example: basis_args={'species':['CO2', 'O2', 'H2O', 'H+']}")
    
    if isinstance(species_args, dict):
        if "species" not in species_args.keys():
            raise Exception("species_args needs a list of species for 'species'. "
                            "Example: species_args={'species':['CO2', 'HCO3-', 'CO3-2']}")
    else:
        raise Exception("species_args needs to be a Python dictionary with a key "
                        "called 'species' (additional keys are optional). "
                        "Example: species_args={'species':['CO2', 'HCO3-', 'CO3-2']}")
    
    basis_sp = basis_args["species"]
    sp = species_args["species"]
    basis(**basis_args)
    species(**species_args)

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
    
    for z in zvals:

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
    
    if not is_predom_plot:

        if 'loga.equil' not in aeout.names:
            yvar = "A/(2.303RT)"
        else:
            yvar = "log a"
        if "alpha" in diagram_args.keys():
            if diagram_args["alpha"]:
                yvar = "alpha"
        
        df_c = pd.concat(dfs)

        fig = px.line(df_c, x=xvar, y="value", color='variable', template="simple_white",
                      width=500,  height=400, animation_frame=anim_var,
                      labels=dict(value=yvar, x=xvar),
                     )

        config = {'displaylogo': False,
                  'modeBarButtonsToRemove': ['resetScale2d', 'toggleSpikelines']}

        fig.show(config=config)
        return
    
    else:
        yvals = __seq(yrange[0], yrange[1], length_out=yres)

    unit_dict = {"P":"bar", "T":"°C", "pH":"", "Eh":"volts", "IS":"mol/kg"}

    for s in basis_sp:
        unit_dict[s] = "logact "+s

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
                    text=html_chemname_format(s),
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
        annotations.append(annotations_i)

        heatmaps_i = go.Heatmap(z=dmaps[i], x=xvals, y=yvals, zmin=0, zmax=len(sp)-1, # is this zmax valid? Double check.
                                customdata=dmaps_names[i],
                                hovertemplate=xvar+': %{x} '+unit_dict[xvar]+'<br>'+yvar+': %{y} '+unit_dict[yvar]+'<br>Region: %{customdata}<extra></extra>')

        heatmaps.append(heatmaps_i)
        
        frame = go.Frame(data=[heatmaps_i],
                         name=str(i),
                         layout=go.Layout(annotations=annotations_i,
                                          coloraxis={'colorscale':["green", "blue", "red", "yellow", "orange", "brown"]}))

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

    fig.update_traces(dict(showscale=False,
                           colorscale='viridis'),
                      selector={'type':'heatmap'})



    fig.update_layout(
        xaxis_title=html_chemname_format(xlab),
        yaxis_title=html_chemname_format(ylab),
        margin={"t": 40, "r":60},
#         updatemenus = [
#         {
#             "buttons": [
#                 {
#                     "args": [None, {"frame": {"duration": 500, "redraw": False},
#                                     "fromcurrent": True, "transition": {"duration": 300,
#                                                                         "easing": "quadratic-in-out"}}],
#                     "label": "Play",
#                     "method": "animate"
#                 },
#                 {
#                     "args": [[None], {"frame": {"duration": 0, "redraw": False},
#                                       "mode": "immediate",
#                                       "transition": {"duration": 0}}],
#                     "label": "Pause",
#                     "method": "animate"
#                 }
#             ],
#             "direction": "left",
#             "pad": {"r": 10, "t": 87},
#             "showactive": False,
#             "type": "buttons",
#             "x": 0.1,
#             "xanchor": "right",
#             "y": 0,
#             "yanchor": "top"
#         }
#     ])
    )

    config = {'displaylogo': False,
              'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d',
                                         'autoScale2d', 'resetScale2d', 'toggleSpikelines',
                                         'hoverClosestCartesian', 'hoverCompareCartesian']}

    fig.show(config=config)



def diagram_interactive(data, title=None,
                        annotation=None, annotation_coords=[0, 0],
                        balance=None, xlab=None, ylab=None, colormap="viridis",
                        width=600, height=520, alpha=False, messages=True,
                        plot_it=True):
    
    """
    Produce an interactive Plotly plot.
    
    Parameters
    ----------
    data : rpy2.ListVector
        Output from `equilibrate` or `affinity`.
    
    title : str, optional
        Title of the plot.
    
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
    
    Returns
    -------
    An interactive plot.
    """
    
    basis_sp = data.rx2("basis").rownames

    xyvars = list(data.rx2("vars"))
    xyvals = list(data.rx2("vals"))

    if 'loga.equil' not in data.names:
        calc_type = "a"
    else:
        calc_type = "e"

    data = equilibrate(data, balance=balance, messages=messages)

    if calc_type=="a":
        out_vals = data.rx2("values")
        out_units = "A/(2.303RT)"
        df = pandas2ri.ri2py_dataframe(out_vals)
        df["n.balance"] = list(data.rx2("n.balance"))
        
        # divide values by balance
        df = df.apply(lambda row: row/row["n.balance"], axis=1)
        df = df.drop(["n.balance"], axis=1)
        sp = list(info([int(val) for val in list(out_vals.names)], messages=False)["name"])
        
    elif calc_type=="e":
        out_vals = data.rx2("loga.equil")
        out_units = "log a"
        df = pandas2ri.ri2py_dataframe(out_vals)
        sp = list(info([int(val) for val in list(data.rx2("values").names)], messages=False)["name"])
    
    df.index = sp
    df = df.transpose()

    if alpha and len(xyvars) == 1:
        df = df.applymap(lambda x: 10**x)
        df = df[sp].div(df[sp].sum(axis=1), axis=0)
        
    xvar = xyvars[0]
    xvals = [float(val) for val in xyvals[0]]

    if len(xyvars) == 1:
        df[xvar] = xvals
        if not isinstance(ylab, str):
            if alpha:
                ylab = "alpha"
            else:
                ylab = out_units
            ylab = html_chemname_format(ylab)
        df = pd.melt(df, id_vars=xyvars, value_vars=sp)

    elif len(xyvars)==2:
        # predominance plot
        yvar = xyvars[1]
        yvals = [float(val) for val in xyvals[1]]
        df["pred"] = df.idxmax(axis=1)
        df["prednames"] = df["pred"]

        xvals_full = xvals*len(yvals)
        yvals_full = __flatten_list([[y]*len(xvals) for y in yvals])
        df[xvar] = xvals_full
        df[yvar] = yvals_full

    unit_dict = {"P":"bar", "T":"°C", "pH":"", "Eh":"volts", "IS":"mol/kg"}

    for s in basis_sp:
        unit_dict[s] = "logact "+s
    
    if not isinstance(xlab, str):
        xlab = xvar+", "+unit_dict[xvar]
        if xvar == "pH":
            xlab = "pH"
        if xvar in basis_sp:
            xlab = unit_dict[xvar]
        xlab = html_chemname_format(xlab)

    if len(xyvars) == 1:
        
        fig = px.line(df, x=xvar, y="value", color='variable', template="simple_white",
                      width=width,  height=height,
                      labels=dict(value=ylab, x=xlab),
                     )
        
        fig.update_layout(xaxis_title=xlab,
                          yaxis_title=ylab,
                          )
        
        if isinstance(title, str):
            fig.update_layout(title={'text':title, 'x':0.5, 'xanchor':'center'})
    
        config = {'displaylogo': False,
                  'modeBarButtonsToRemove': ['resetScale2d', 'toggleSpikelines']}

    if len(xyvars) == 2:
        mappings = {'pred': {s:lab for s,lab in zip(sp,range(0,len(sp)))}}
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
            ylab = html_chemname_format(ylab)
        
        fig = px.imshow(dmap, width=width, height=height, aspect="auto",
                        labels=dict(x=xlab, y=ylab, color="region"),
                        x=xvals, y=yvals, template="simple_white",
                       )

        fig.update(data=[{'customdata': dmap_names,
            'hovertemplate': xlab+': %{x}<br>'+ylab+': %{y}<br>Region: %{customdata}<extra></extra>'}])

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
                fig.add_annotation(x=namex, y=namey,
                                   text=html_chemname_format(s),
                                   bgcolor="rgba(255, 255, 255, 0.5)",
                                   showarrow=False)

        config = {'displaylogo': False,
                  'modeBarButtonsToRemove': ['zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d',
                                             'autoScale2d', 'resetScale2d', 'toggleSpikelines',
                                             'hoverClosestCartesian', 'hoverCompareCartesian']}
        
    if isinstance(annotation, str):
        fig.add_annotation(
            x=annotation_coords[0],
            y=annotation_coords[1],
            text=annotation,
            showarrow=False,
            xref="paper",
            yref="paper",
            bgcolor="rgba(255, 255, 255, 0.5)")
        
    if plot_it:
        fig.show(config=config)
    
    return df


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
            messages=True, interactive=False,
            annotation=None, annotation_coords=[0,0]):
    
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
    
    interactive : bool, default False
        Experimental! Display an interactive plot?

    annotation : str, optional
        Annotation to add to the plot. Interactive plots only (`interactive`
        is set to True).
    
    annotation_coords : list of numeric, default [0, 0], optional
        Coordinates of annotation, where 0,0 is bottom left and 1,1 is top
        right. Interactive plots only (`interactive` is set to True).
    
    Returns
    -------
    a : rpy2.ListVector
        Output from `diagram`.
    args : dict
        Dictionary of arguments supplied to `diagram`.
    """
    
    if interactive:
        df = diagram_interactive(data=eout, title=main, annotation=annotation,
                                 annotation_coords=annotation_coords,
                                 balance=balance,
                                 xlab=xlab, ylab=ylab,
                                 colormap=fill,
                                 width=width, height=height,
                                 alpha=alpha, plot_it=plot_it,
                                 messages=messages)
        return df
    
    
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
        with __r_inline_plot(width=width, height=height, dpi=dpi, plot_it=plot_it):
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

        
def unicurve(logK, species, phase, stoich, pressures=1, temperatures=25,
             minT=0.1, maxT=100, minP=1, maxP=500,
             solve="T", width=600, height=520, dpi=90, plot_it=True,
             messages=True, show=True):
    
    species = _convert_to_RVector(species, force_Rvec=False)
    phase = _convert_to_RVector(phase, force_Rvec=False)
    stoich = _convert_to_RVector(stoich, force_Rvec=False)
    pressures = _convert_to_RVector(pressures, force_Rvec=False)
    temperatures = _convert_to_RVector(temperatures, force_Rvec=False)
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r_univariant = pkg_resources.resource_string(
            __name__, 'univariant.r').decode("utf-8")
        ro.r(r_univariant)
        if solve=="T":
            a = ro.r.uc_solveT(logK=logK,
                               species=species,
                               phase=phase,
                               stoich=stoich,
                               pressures=pressures,
                               minT=minT,
                               maxT=maxT)
            with __r_inline_plot(width=width, height=height, dpi=dpi, plot_it=plot_it):
                ro.r.create_output_plot_T(logK=logK,
                                          species=species,
                                          phase=phase,
                                          stoich=stoich,
                                          pressures=pressures,
                                          minT=minT,
                                          maxT=maxT)
        elif solve=="P":
            a = ro.r.uc_solveP(logK=logK,
                               species=species,
                               phase=phase,
                               stoich=stoich,
                               temperatures=temperatures,
                               minP=minP,
                               maxP=maxP)
            with __r_inline_plot(width=width, height=height, dpi=dpi, plot_it=plot_it):
                ro.r.create_output_plot_P(logK=logK,
                                          species=species,
                                          phase=phase,
                                          stoich=stoich,
                                          temperatures=temperatures,
                                          minP=minP,
                                          maxP=maxP)
    if messages:
        for warning in w:
            print(warning.message)

    if len(a) == 3:
        warn = a[2][0] # subcrt's list includes warnings only if they appear
    else:
        warn = None
    
    out_dict = {"reaction":pandas2ri.ri2py_dataframe(a[0]),
                "out":pandas2ri.ri2py_dataframe(a[1])} # the extra [0] is important
        
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