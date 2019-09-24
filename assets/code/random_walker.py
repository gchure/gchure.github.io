"""
The Random Walker 
--------------------------------------------------------------------------------
Author: Griffin Chure
Last Modified: September 23, 2019
License: MIT

Description
--------------------------------------------------------------------------------
This file generates a standalone HTML file of a random walker with 100,000 steps
in two dimensions. The left plot shows the overview of the entire random walk
while the right plot shows the current position and the last 5000 steps of the
walker

"""
import numpy as np
import bokeh.io
import bokeh.plotting
import bokeh.layouts
import bokeh.palettes
from bokeh.models import *
bokeh.io.output_file('random_walk.html')
colors = {'black':'#444147', 'purple': '#7E59A2', 'orange':'#E39943'} 

# Define the initial number of steps
n_steps = int(1E5)

# Assemble the three necessary data sources
source = ColumnDataSource({'x':[], 'y':[]})
current_position = ColumnDataSource({'x':[], 'y':[]})
display_source = ColumnDataSource({'x':[], 'y':[]})

# Define the interactive components
run_button = Button(label='click to generate a random walk')
step_slider = Slider(title='number of steps', start=10, end=int(n_steps), step=1,
                    value=10, bar_color=colors['purple'])

# Define the javascript for each component
generate_walk = """

// Instantiate position vectors
var xs = [0]
var ys = [0]

// Loop through each step, starting with the second.
for (var i = 1; i < n_steps; i++) { 

    // Determine a random angle to step towards
    theta = Math.random() * 2 * Math.PI;

    // Compute and store the new xy positions
    xs.push(xs[i - 1] + Math.cos(theta)); 
    ys.push(ys[i - 1] + Math.sin(theta)); }


// Update the data source
source.data['x'] = xs;
source.data['y'] = ys;
source.change.emit()

// Reset the step slider to the default value
stepSlider.value = 10;
"""

display_steps = """
// Determine where to start slicing the inset data
if (stepSlider.value < 5000) {
    var init = 0;
}
else { 
    var init = stepSlider.value - 5000;
}

// Assign the inset data for display
display_source.data['x'] = source.data['x'].slice(init, stepSlider.value);
display_source.data['y'] = source.data['y'].slice(init, stepSlider.value);

// Highlight the current position -- Can do this with an indexfilter as well
current_position.data['x'] = display_source.data['x'].slice(-1);
current_position.data['y'] = display_source.data['y'].slice(-1);

// Update the data sources
display_source.change.emit()
current_position.change.emit()
"""


# Define and assignthe callbacks
args={'stepSlider':step_slider, 'source':source,
'display_source':display_source, 'current_position':current_position,
'n_steps':n_steps}
generate_cb = CustomJS(args=args, code=generate_walk + display_steps)
display_cb = CustomJS(args=args, code=display_steps)
run_button.js_on_click(generate_cb)
step_slider.js_on_change('value', display_cb)

# Define the axes
overview_ax = bokeh.plotting.figure(width=300, height=300, match_aspect=True,
                x_axis_label = 'x position', y_axis_label='y position')
inset_ax = bokeh.plotting.figure(width=300, height=300, match_aspect=True,
                x_axis_label='x position', y_axis_label = 'y position')

# Populate teh axes
overview_ax.circle('x', 'y', source=current_position, color=colors['orange'],
level='overlay', size=5)
overview_ax.square('x', 'y', source=current_position, fill_color='grey', fill_alpha=0.25,
line_color=colors['orange'], level='overlay', size=30)
overview_ax.line('x', 'y', source=source, color=colors['black'], line_width=0.5,
alpha=0.5)

inset_ax.line('x', 'y', source=display_source, color=colors['purple'], line_width=0.5, alpha=1)
inset_ax.circle('x', 'y', source=current_position, color=colors['orange'],
level='overlay', size=5)


# Define the layout
row = bokeh.layouts.row(overview_ax, inset_ax)
lay = bokeh.layouts.column(run_button, step_slider, row)


# Define the theme and save
theme_json =  {
    'attrs' : {
        'Figure' : {
            'background_fill_color': '#EEEEEE',
        },
        'Axis': {
            'axis_line_color': 'slategray',
            'major_tick_line_color': None,
            'minor_tick_line_color': None,
        },
        'Legend': {
            'border_line_color': 'slategray',
            'background_fill_color': '#EEEEEE',
            'border_line_width': 0.75,
            'background_fill_alpha': 0.75,
        },
        'Grid': {
            'grid_line_color': '#FFFFFF',
            'grid_line_width': 0.75,
        },
        'Text': {
            'text_font_style': 'italic',
            'text_font': 'Arial', 
            'text_font_size':10,
        },
        'Title': {
            'text_color': '#3c3c3c',
            'align': 'left',
            'text_font': 'Arial',
            'text_font_style': 'italic',
            'offset': 2,
         }
    }
}
       
theme = bokeh.themes.Theme(json=theme_json)
bokeh.io.curdoc().theme = theme
bokeh.io.save(lay)